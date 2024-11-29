#ifndef GAUSSAll_H
#define GAUSSAll_H

#include <vector>
#include <pthread.h>
#include <cmath>
#include <algorithm>
#include <mutex>

void* gaussianStep(void* arg);
int gaussianElimination(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<double>& x, int n, int num_threads);

#endif
