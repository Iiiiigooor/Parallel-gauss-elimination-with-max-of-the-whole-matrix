#include <iostream>
#include <chrono>
#include <pthread.h>
#include "matrix.h"
#include "gaussAll.h"
#include "functions.h"

double residualNorm = 0.0;
double normError = 0.0;

// data for calculateResidualNorm
struct ResidualNormData {
    const std::vector<std::vector<double>>& A;
    const std::vector<double>& x;
    const std::vector<double>& b;
    int n;
    int num_threads;
};

// data for calculateNormError
struct NormErrorData {
    const std::vector<double>& x;
    int n;
    int num_threads;
};

// Функция потока для вычисления нормы невязки
void* calculateResidualNormThread(void* arg) {
    ResidualNormData* data = static_cast<ResidualNormData*>(arg);
    residualNorm = calculateResidualNorm(data->A, data->x, data->b, data->n, data->num_threads);
    return nullptr;
}

void* calculateNormErrorThread(void* arg) {
    NormErrorData* data = static_cast<NormErrorData*>(arg);
    normError = calculateNormError(data->x, data->n, data->num_threads);
    return nullptr;
}

int main(int argc, char* argv[]) {
    if (argc < 5 || argc > 6) {
        std::cerr << "Usage: " << argv[0] << " n p m k filename" << std::endl;
        return EXIT_FAILURE;
    }

    int n = std::stoi(argv[1]);       // matrix dimension
    int p = std::stoi(argv[2]);       // amount of threads
    int m = std::stoi(argv[3]);       // amount of output values
    int k = std::stoi(argv[4]);       // formula number

    if (m > n) {
        std::cerr << "Error:Number of values to be output is greater than the matrix dimension " << std::endl;
        return 1;
    }

    if (p <= 0) {
        std::cerr << "Error: Incorrect amount of threads " << std::endl;
        return 2;
    }

    int err = 0;
    std::vector<std::vector<double>> A;
    std::vector<double> b, x;
    A.resize(n, std::vector<double>(n));
    b.resize(n);

    if (k == 0) {
        std::string filename = argv[4];
        err = readMatrixFromFile(filename, A, n);
        if (!err) {
            std::cerr << "Error: Can't open file " << filename << std::endl;
            return err;
        }
    }
    else {
        initializeMatrix(A, k, n, p);
    }

    std::cout << "Initial matrix A:" << std::endl;
    printMatrix(A, m);

    // constructing b
    for (int i = 0; i < n; i++) {
        double sum_value = 0.0;
        for (int k = 0; (2 * k + 1) < n; k++) {
            sum_value += A[i][2 * k + 1];
        }
        b[i] = sum_value;
    }

    std::cout << "The right side b:" << std::endl;
    printVector(b, m);

    auto start = std::chrono::high_resolution_clock::now();

    err = gaussianElimination(A, b, x, n, p);
    if (!err) {
        std::cerr << "Error: Matrix is singular " << std::endl;
        return err;
    }

    // Остановка таймера
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Решение x:" << std::endl;
    printVector(x, m);

    // data for threads
    ResidualNormData residualNormData = {A, x, b, n, std::max((p - 2) / 2, 1)};
    NormErrorData normErrorData = {x, n, std::max((p - 2) / 2, 1)};

    pthread_t residualNormThread, normErrorThread;
    pthread_create(&residualNormThread, nullptr, calculateResidualNormThread, &residualNormData);
    pthread_create(&normErrorThread, nullptr, calculateNormErrorThread, &normErrorData);


    pthread_join(residualNormThread, nullptr);
    pthread_join(normErrorThread, nullptr);

    // results
    std::cout << "Норма невязки: " << std::scientific << residualNorm << std::endl;
    std::cout << "Норма погрешности: " << std::scientific << normError << std::endl;

    // time
    std::cout << "Время решения: " << elapsed.count() << " секунд" << std::endl;

    return 0;
}
