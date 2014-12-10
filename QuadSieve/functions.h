//
//  functions.h
//  QuadSieve
//
//  Created by Christina Black on 12/10/14.
//  Copyright (c) 2014 Christina Black. All rights reserved.
//

#ifndef QuadSieve_functions_h
#define QuadSieve_functions_h
#include <vector>

struct tonelli_pair {
    long long r;
    long long p_r;
    long long prime;
};
long mod(unsigned long long a, unsigned long long b);
long long modpow(unsigned long long base, unsigned long long exponent, unsigned long long modulus);
long long power(long long i, long long j);
int jacobi(long long a, long long m);
tonelli_pair tonelli(long long a, long long p);
void loadPrimes(long *primes);
long ** createMatrix(long rows, long cols);
void initialize(long ** mat, long rows, long cols);
void destroyMatrix(long **matrix, long size);
void printToFile(long ** matrix, long rows, long cols);
void print(long **matrix, long rows, long cols);
void print(std::vector<long long> matrix, long rows);


#endif
