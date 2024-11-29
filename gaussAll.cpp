#include "gaussAll.h"

std::mutex mtx;

struct Barrier {
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    int count;
    int trip_count;

    void init(int num_threads) {
        pthread_mutex_init(&mutex, nullptr);
        pthread_cond_init(&cond, nullptr);
        count = 0;
        trip_count = num_threads;
    }

    void wait() {
        pthread_mutex_lock(&mutex);
        count++;
        if (count >= trip_count) {
            count = 0; // Reset for reuse
            pthread_cond_broadcast(&cond);
        } else {
            pthread_cond_wait(&cond, &mutex);
        }
        pthread_mutex_unlock(&mutex);
    }

    void destroy() {
        pthread_mutex_destroy(&mutex);
        pthread_cond_destroy(&cond);
    }
};

Barrier barrier;

struct ThreadData {
    std::vector<std::vector<double>>& A;
    std::vector<double>& b;
    std::vector<double>& x;
    int n;
    int tid; // Thread ID
    int num_threads;
    int* globalMaxRow; // Pointer to the global max row index
    int* globalMaxCol; // Pointer to the global max column index
    std::vector<int>& columnOrder; // Column order tracker

    ThreadData(std::vector<std::vector<double>>& A_, std::vector<double>& b_, std::vector<double>& x_, int n_, int tid_, int num_threads_, int* globalMaxRow_, int* globalMaxCol_, std::vector<int>& columnOrder_)
        : A(A_), b(b_), x(x_), n(n_), tid(tid_), num_threads(num_threads_), globalMaxRow(globalMaxRow_), globalMaxCol(globalMaxCol_), columnOrder(columnOrder_) {}
};

void* gaussianStep(void* arg) {
    ThreadData* data = static_cast<ThreadData*>(arg);
    int n = data->n;
    int tid = data->tid;
    int num_threads = data->num_threads;

    for (int step = 0; step < n; ++step) {
        *data->globalMaxRow = step;
        *data->globalMaxCol = step;

        int localMaxRow = step, localMaxCol = step;
        for (int i = step + tid; i < n; i += num_threads) {
            for (int j = step; j < n; ++j) {
                if (std::fabs(data->A[i][j]) > std::fabs(data->A[localMaxRow][localMaxCol])) {
                    localMaxRow = i;
                    localMaxCol = j;
                }
            }
        }

        // Update global maximum row and column
        mtx.lock();
        if (std::fabs(data->A[localMaxRow][localMaxCol]) > std::fabs(data->A[*data->globalMaxRow][*data->globalMaxCol])) {
            *data->globalMaxRow = localMaxRow;
            *data->globalMaxCol = localMaxCol;
        }
        mtx.unlock();

        barrier.wait();

        if (tid == 0) {
            if (std::fabs(data->A[*data->globalMaxRow][*data->globalMaxCol]) < 1e-19) {
                pthread_exit((void*)0); // Singular matrix detected
            }

            std::swap(data->A[step], data->A[*data->globalMaxRow]);
            std::swap(data->b[step], data->b[*data->globalMaxRow]);

            // Swap columns in A and track in columnOrder
            for (int i = 0; i < n; ++i) {
                std::swap(data->A[i][step], data->A[i][*data->globalMaxCol]);
            }
            std::swap(data->columnOrder[step], data->columnOrder[*data->globalMaxCol]);
        }

        barrier.wait();

        for (int i = step + 1 + tid; i < n; i += num_threads) {
            double factor = data->A[i][step] / data->A[step][step];
            data->b[i] -= factor * data->b[step];
            for (int j = step; j < n; ++j) {
                data->A[i][j] -= factor * data->A[step][j];
            }
        }

        barrier.wait();
    }

    for (int i = n - 1; i >= 0; --i) {
        if (tid == 0) {
            data->x[data->columnOrder[i]] = data->b[i] / data->A[i][i];
        }

        barrier.wait();

        for (int j = i - 1 - tid; j >= 0; j -= num_threads) {
            data->b[j] -= data->A[j][i] * data->x[data->columnOrder[i]];
        }

        barrier.wait();
    }

    return nullptr;
}

int gaussianElimination(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<double>& x, int n, int num_threads) {
    x.resize(n);
    std::vector<std::vector<double>> A_copy = A;
    std::vector<double> b_copy = b;

    // Initialize the custom barrier
    barrier.init(num_threads);

    int globalMaxRow = 0;
    int globalMaxCol = 0;

    std::vector<int> columnOrder(n);
    for (int i = 0; i < n; ++i) columnOrder[i] = i;

    std::vector<pthread_t> threads(num_threads);
    std::vector<ThreadData> thread_data;
    thread_data.reserve(num_threads);

    for (int tid = 0; tid < num_threads; ++tid) {
        thread_data.emplace_back(A_copy, b_copy, x, n, tid, num_threads, &globalMaxRow, &globalMaxCol, columnOrder);
        pthread_create(&threads[tid], nullptr, gaussianStep, &thread_data[tid]);
    }

    for (int tid = 0; tid < num_threads; ++tid) {
        pthread_join(threads[tid], nullptr);
    }

    barrier.destroy();

    return 1;
}
