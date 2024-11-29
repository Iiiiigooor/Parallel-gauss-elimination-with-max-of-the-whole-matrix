#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <cmath>
#include <pthread.h>
#include <vector>

void* calculateAx(void* arg);
void* calculateNorms(void* arg);
double calculateResidualNorm(const std::vector<std::vector<double>>& A, const std::vector<double>& x, const std::vector<double>& b, const int n, const int num_threads);
void* calculatePartialNormError(void* arg);
double calculateNormError(const std::vector<double>& x, const int n, const int num_threads);

#endif
