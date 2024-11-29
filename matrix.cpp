#include "matrix.h"

int readMatrixFromFile(const std::string& filename, std::vector<std::vector<double>>& A, int n) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        return 0;
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            file >> A[i][j];
        }
    }
    file.close();
    return 1;
}

void printMatrix(const std::vector<std::vector<double>>& A, int m) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            std::cout << std::setw(10) << std::setprecision(3) << std::scientific << A[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void printVector(const std::vector<double>& vec, int m) {
    for (int i = 0; i < m; ++i) {
        std::cout << std::setw(10) << std::setprecision(3) << std::scientific << vec[i] << " ";
    }
    std::cout << std::endl;
}

double f(int k, int n, int i, int j) {
    switch (k) {
    case 1:
        return n - std::max(i, j) + 1;
    case 2:
        return std::max(i, j);
    case 3:
        return std::abs(i - j);
    case 4:
        return 1.0 / (i + j - 1);
    default:
        throw std::invalid_argument("Error: Wrong formula number");
    }
}

// data for one thread
struct ThreadData {
    std::vector<std::vector<double>>* A;
    int k, n;
    int start_row, end_row; 
};

void* initializeMatrixThread(void* arg) {
    ThreadData* data = static_cast<ThreadData*>(arg);

    for (int i = data->start_row; i < data->end_row; ++i) {
        for (int j = 0; j < data->n; ++j) {
            (*data->A)[i][j] = f(data->k, data->n, i + 1, j + 1);
        }
    }
    return nullptr;
}

void initializeMatrix(std::vector<std::vector<double>>& A, int k, int n, int num_threads) {
    std::vector<pthread_t> threads(num_threads);
    std::vector<ThreadData> thread_data(num_threads);

    int rows_per_thread = n / num_threads;
    int remaining_rows = n % num_threads;
    int start_row = 0;

    // Creating threads with individual range
    for (int tid = 0; tid < num_threads; ++tid) {
        int end_row = start_row + rows_per_thread + (tid < remaining_rows ? 1 : 0);
        thread_data[tid] = {&A, k, n, start_row, end_row};  // Используем указатель на A
        pthread_create(&threads[tid], nullptr, initializeMatrixThread, &thread_data[tid]);
        start_row = end_row;
    }

    for (int tid = 0; tid < num_threads; ++tid) {
        pthread_join(threads[tid], nullptr);
    }
}