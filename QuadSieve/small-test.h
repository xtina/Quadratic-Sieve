//
//  small-test.h
//  QuadSieve
//
//  Created by Christina Black on 12/10/14.
//  Copyright (c) 2014 Christina Black. All rights reserved.
//

#ifndef __QuadSieve__small_test__
#define __QuadSieve__small_test__

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

struct tonelli_pair {
    long long r;
    long long p_r;
    long long prime;
};

long mod(unsigned long long a, unsigned long long b){
    return (a%b+b)%b;
}

//modular exponentiation
long long modpow(unsigned long long base, unsigned long long exponent, unsigned long long modulus) {
    
    unsigned long long result = 1;
    
    while (exponent > 0) {
        if ((exponent & 1) == 1) {
            // multiply in this bit's contribution while using modulus to keep result small
            result = (result * base) % modulus;
        }
        // move to the next bit of the exponent, square (and mod) the base accordingly
        exponent >>= 1;
        base = (base * base) % modulus;
    }
    
    return result;
}

//raise i to the power j
long long power(long long i, long long j){
    long long int temp=1;
    for(int k=0; k<j; k++){
        temp*=i;
    }
    return temp;
}

//jacobi symbol calculator
int jacobi(long long a, long long m){
    a = mod(a,m);
    int t = 1;
    long c;
    long temp;
    
    while(a != 0){
        c = 0;
        while(a%2 == 0){
            a /= 2;
            c=1-c; //if c is 1, exponent is odd
        }
        if(c == 1)
            if(m%8 == 3 || m%8 == 5)
                t *= -1;
        if(a%4 == 3 && m%4 == 3){
            t *= -1;
        }
        temp = m;
        m = a;
        a = temp;
        a = mod(a, m);
    }
    
    if(m==1)
        return t;
    
    return 0;
}

//tonelli
tonelli_pair tonelli(long long a, long long p){
    long b=1;
    long long t = p-1;
    long long s=0;
    int temp;
    tonelli_pair pair;
    
    if(jacobi(a, p)  == -1) {
        cout << "uh oh" << endl;
    }
    
    //find quad non res
    temp = 0;
    while(temp == 1 || temp == 0){
        temp = jacobi(b, p);
        b++;
    }
    
    //pull out factors of 2
    while(t%2 == 0){
        t /= 2;
        s++;
    }
    
    long long i = 2;
    long long c = mod(a*b*b, p);
    for(int k=1; k < s-1; k++){
        if(modpow(c, power(2, s-k-1)*t, p) == -1){
            i += power(2, k);
            c = mod(c*power(b, power(2, k)), p);
        }
    }
    pair.r = modpow(b, (i*t)/2, p)*modpow(a, (t+1)/2, p);
    pair.p_r = p-pair.r;
    
    while (pair.p_r < 0) {
        pair.p_r += p;
    }
    
    pair.prime = p;
    
    return pair;
}

void loadPrimes(long *primes){
    ifstream infile;
    int prime;
    int i=0;
    infile.open("/Users/Christina/Documents/QuadSieve/QuadSieve/primes.txt");
    while (infile >> prime) {
        primes[i++]=prime;
    }
    infile.close();
}

long ** createMatrix(long rows, long cols) {
    long **mat;
    mat = new long*[rows];
    for(int i = 0; i < rows; ++i) {
        mat[i] = new long[cols];
    }
    return mat;
}

void initialize(long ** mat, long rows, long cols) {
    for(int m=0; m < rows; m++)
        for(int n=0; n < cols; n++) {
            mat[m][n] = 0;
        }
}

void destroyMatrix(long **matrix, long size) {
    for(int i=0; i < size; i++)
        delete [] matrix[i];
    delete [] matrix;
}

void print(long ** matrix, long rows, long cols){
    ofstream out;
    out.open("trial-division-matrix.txt");
    for(int i=0; i < rows; i++){
        for(int j=0; j < cols; j++){
            out << matrix[i][j] << " ";
        }
        out << endl;
    }
    out.close();
}

void print(vector<long long> matrix, long rows) {
    for(int i=0; i < rows; i++) {
        cout << matrix[i] << " ";
    }
}


#endif /* defined(__QuadSieve__small_test__) */
