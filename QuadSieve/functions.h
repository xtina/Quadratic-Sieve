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
#include "BigInteger.hh"

struct tonelli_pair {
    BigInteger r;
    BigInteger p_r;
    BigInteger prime;
};
BigInteger intSqrt (BigInteger remainder);
BigInteger mod(BigInteger a, BigInteger b);
BigInteger modpow(BigInteger base, BigInteger exponent, BigInteger modulus);
BigInteger power(BigInteger i, BigInteger j);
BigInteger jacobi(BigInteger a, BigInteger m);
tonelli_pair tonelli(BigInteger a, BigInteger p);
void loadPrimes(BigInteger *primes);
BigInteger ** createMatrix(unsigned long rows, unsigned long cols);
void initialize(BigInteger ** mat, unsigned long rows, unsigned long cols);
void destroyMatrix(BigInteger **matrix, BigInteger size);
void printToFile(BigInteger ** matrix, BigInteger rows, BigInteger cols);
void print(BigInteger **matrix, BigInteger rows, BigInteger cols);
void print(std::vector<tonelli_pair> pair, BigInteger size);
void print(std::vector<BigInteger> matrix, BigInteger rows);
void print(BigInteger *matrix, int num);
int ** intMatrix(unsigned long rows, unsigned long cols);
void initialize(int **mat, unsigned long rows, unsigned long cols);
void initialize(int *mat, unsigned long rows);
void print(int **mat, unsigned long rows, unsigned long cols);
void destroyMatrix(int **matrix, unsigned long size);
void swap_rows(int **A, unsigned long i, unsigned long k, unsigned long size);
void to_reduced_row_echelon_form(int **A, unsigned long rows, unsigned long cols);
BigInteger findZeroRow(int ** mat, unsigned long rows, unsigned long cols, unsigned long startIndex=0);
void addRows(int* a, int* b, unsigned long i);
void copyMatrix(int** a, int**b, unsigned long rows, unsigned long cols);
BigInteger gcd(BigInteger a, BigInteger b);


#endif
