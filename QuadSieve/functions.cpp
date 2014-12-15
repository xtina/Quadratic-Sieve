//
//  functions.cpp
//  QuadSieve
//
//  Created by Christina Black on 12/10/14.
//  Copyright (c) 2014 Christina Black. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "BigInteger.hh"
#include "functions.h"

using namespace std;

BigInteger mod(BigInteger a, BigInteger b) {
    return (a % b + b) % b;
}

//modular exponentiation
BigInteger modpow(BigInteger base, BigInteger exponent, BigInteger modulus) {
    BigInteger b = base;
    BigInteger n = 1;
    while (exponent != 0) {
        if (mod(exponent, 2) == 1) {
            n = mod(n * base, modulus);
        }
        exponent = (unsigned long)floor((exponent / BigInteger(2)).toUnsignedLong());
        base = (base * base) % modulus;
    }
    
    return n;
    
}

//raise i to the power j
BigInteger power(BigInteger i, BigInteger j) {
    BigInteger temp = 1;
    for (BigInteger k = 0; k < j; k++) {
        temp *= i;
    }
    return temp;
}

//jacobi symbol calculator
BigInteger jacobi(BigInteger a, BigInteger m) {
    a = mod(a, m);
    BigInteger t = 1;
    BigInteger c;
    BigInteger temp;
    
    while (a != 0) {
        c = 0;
        while (a % 2 == 0) {
            a /= 2;
            c = BigInteger(1) - c; //if c is 1, exponent is odd
        }
        if (c == 1)
            if (m % 8 == 3 || m % 8 == 5)
                t *= -1;
        if (a % 4 == 3 && m % 4 == 3) {
            t *= -1;
        }
        temp = m;
        m = a;
        a = temp;
        a = mod(a, m);
    }
    
    if (m == 1)
        return t;
    
    return 0;
}

//tonelli
tonelli_pair tonelli(BigInteger a, BigInteger p) {
    BigInteger b = 0;
    BigInteger t = p - 1;
    BigInteger s = 0;
    BigInteger temp;
    tonelli_pair pair;
    
    if (jacobi(a, p)  == -1) {
        cout << "uh oh" << endl;
    }
    
    //find quad non res
    temp = 0;
    while (temp == 1 || temp == 0) {
        b++;
        temp = jacobi(b, p);
    }
    
    //pull out factors of 2
    while (mod(t, 2) == 0) {
        t /= 2;
        s++;
    }
    BigInteger i = 2;
    BigInteger c = mod(a * b * b, p);
    for (BigInteger k = 1; k <= s - 1; k++) {
        if (modpow(c, power(2, s - k - 1)*t, p) == p - 1) {
            i += power(2, k);
            c = mod(c * power(b, power(2, k)), p);
        }
    }
    
    pair.r = mod(modpow(b, (i * t) / 2, p) * modpow(a, (t + 1) / 2, p), p);
    pair.p_r = p - pair.r;
    
    while (pair.p_r < 0) {
        pair.p_r += p;
    }
    
    pair.prime = p;
    return pair;
}

void loadPrimes(BigInteger *primes) {
    ifstream infile;
    unsigned long temp;
    int i = 0;
    infile.open("/Users/Christina/Documents/QuadSieve/QuadSieve/primes.txt");
    while (infile >> temp) {
        
        primes[i++] = BigInteger(temp);
    }
    infile.close();
}

BigInteger **createMatrix(unsigned long rows, unsigned long cols) {
    BigInteger **mat;
    mat = new BigInteger*[rows];
    for (BigInteger i = 0; i < rows; ++i) {
        mat[i.toInt()] = new BigInteger[cols];
    }
    return mat;
}

long long ** intMatrix(unsigned long rows, unsigned long cols) {
    long long **mat;
    mat = new long long *[rows];
    for (int i = 0; i < rows; ++i) {
        mat[i] = new long long[cols];
    }
    return mat;
}

void initialize(long **mat, unsigned long rows, unsigned long cols) {
    for (BigInteger m = 0; m < rows; m++)
        for (BigInteger n = 0; n < cols; n++) {
            mat[m.toInt()][n.toInt()] = 0;
        }
}

void initialize(BigInteger ** a, unsigned long b, unsigned long c) {
    for(int i=0; i < b; i++) {
        for(int j=0; j < c; j++) {
            a[i][j] = 0;
        }
    }
}

void initialize(long long *mat, unsigned long rows) {
    for (unsigned long m = 0; m < rows; m++)
            mat[m] = 0;
}

void initialize(long long **mat, unsigned long a, unsigned long b) {
    for(int i=0; i<a; i++) {
        for(int j=0; j<b; j++) {
            mat[i][j] = 0;
        }
    }
}

void initialize(int *mat, unsigned long rows) {
    for (unsigned long m = 0; m < rows; m++)
            mat[m] = 0;
}

void destroyMatrix(BigInteger **matrix, BigInteger size) {
    for (BigInteger i = 0; i < size; i++)
        delete [] matrix[i.toInt()];
    delete [] matrix;
}

void destroyMatrix(long long **matrix, unsigned long size) {
    for (BigInteger i = 0; i < size; i++)
        delete [] matrix[i.toInt()];
    delete [] matrix;
}

void printToFile(BigInteger **matrix, BigInteger rows, BigInteger cols) {
    ofstream out;
    out.open("trial-division-matrix.txt");
    for (BigInteger i = 0; i < rows; i++) {
        for (BigInteger j = 0; j < cols; j++) {
            out << matrix[i.toInt()][j.toInt()].toUnsignedLong() << " ";
        }
        out << endl;
    }
    out.close();
}

void print(BigInteger **matrix, BigInteger rows, BigInteger cols) {
    for (BigInteger i = 0; i < rows; i++) {
        for (BigInteger j = 0; j < cols; j++) {
            cout << matrix[i.toInt()][j.toInt()].toLong() << " ";
        }
        cout << endl;
    }
}

void print(vector<tonelli_pair> pair, BigInteger size) {
    for (BigInteger i = 0; i < size; i++) {
        cout << pair[i.toInt()].r.toLong() << " " << pair[i.toInt()].p_r.toUnsignedLong() << " " << pair[i.toInt()].prime.toUnsignedLong() << endl;
    }
}

void print(BigInteger *matrix, int num) {
    for (int i = 0; i < num; i++) {
        cout << matrix[i].toLong() << " ";
    }
}

void print(long long **mat, unsigned long rows, unsigned long cols) {
    for(int i=0; i < rows; i++){
        for(int j=0; j < cols; j++) {
            cout << mat[i][j];
        }
        cout << endl;
    }
}

void print(vector<BigInteger> matrix, BigInteger rows) {
    for (BigInteger i = 0; i < rows; i++) {
        cout << matrix[i.toInt()].toLong() << " ";
    }
}

// Swap rows i and k of a matrix A
// Note that due to the reference, both dimensions are preserved for
// built-in arrays
void swap_rows(long long **A, unsigned long i, unsigned long k, unsigned long size) {
    
    // check indices
    for (unsigned long col = 0; col < size; ++col) {
        std::swap(A[i][col], A[k][col]);
    }
}

// convert A to reduced row echelon form
void to_reduced_row_echelon_form(long long **A, unsigned long rows, unsigned long cols) {
    

    unsigned long lead = 0;
    
    for (unsigned long j = 0;j < rows; ++j) {
        unsigned long i = j;
        while (A[i][lead] == 0 || A[i][lead]%2 == 0) {
            ++i;
            if (i > rows-1) {
                i = j;
                ++lead;
                if (lead > cols)
                    return;
            }
        }
        
        if(i != j) //if they're not the same row, swap them
            swap_rows(A, i, j, cols);
        
        for (unsigned long k = j+1; k < rows; k++) {
            if(A[k][lead]%2 == 1 && A[k][lead] != 0) { //add rows w/ leading 1 together
                for(int a=0; a < cols; a++) {
//                    A[k][a] = (A[k][a] + A[j][a])%2;
                    A[k][a] = A[k][a] + A[j][a];
                }
            }
        }
    }
}


//add row b to row a
void addRows(long long* a, long long * b, unsigned long i) {
    for(unsigned long j=0; j< i; j++) {
        a[j] += b[j];
    }
}

//copy b into a
void copyMatrix(long long ** a, long long **b, unsigned long rows, unsigned long cols) {
    for(int i=0; i < rows; i++){
        for(int j=0; j < cols; j++) {
            a[i][j] = b[i][j];
        }
    }
}

BigInteger findZeroRow(long long ** mat, unsigned long rows, unsigned long cols, unsigned long startIndex) {
    unsigned long j=0;
    for(unsigned long i=startIndex; i < rows; i++){
        while(j < cols) { 
            if (mat[i][j]%2 == 1)
                break;
            j++;
        }
        if(j == cols)
            return i;
        j=0;
    }
    cout << "no zero row found" << endl;
    return -1;
}

BigInteger gcd1(BigInteger a, BigInteger b){
    if(a==0 && b==0){
        std::cout << "a and b cannot both be 0" << std::endl;
        return 0;
    }
    BigInteger temp;
    while (b != 0) {
        temp = b;
        b = a%b;
        a=temp;
    }
    if(a < 0)
        a *=-1;
    return a;
}
