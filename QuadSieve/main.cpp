//
//  main.cpp
//  QuadSieve
//
//  Created by Christina on 12/2/14.
//  Copyright (c) 2014 Christina. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <math.h>
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

void print(vector<long long> matrix, long rows){
    for(int i=0; i < rows; i++) {
            cout << matrix[i] << " ";
    }
}


int main(int argc, const char * argv[]) {
    
    long long n = 13592675504123;
    long primes[18383];
    list<long> quadPrimes; //primes that are quad res of my prime
    vector<tonelli_pair> tonelliNums;
   
    long M=10000;
    tonelli_pair temp;
    long long sr;
    
    loadPrimes(primes);
    
    //find quad res n mod p
    for(int i=1; i < 18383; i++){
        if(jacobi(n,primes[i])==1){
            quadPrimes.push_back(primes[i]);
        }
    }
    
    //find tonelli pairs w.r.t. n mod p
    for(list<long>::const_iterator
        iterator = quadPrimes.begin(),
        end = quadPrimes.end();
        iterator != end; ++iterator)
    { //j starts at 1 to skip prime 2
        temp = tonelli(n, *iterator);
        tonelliNums.push_back(temp);
        //cout << "(" << temp.r << ", " << temp.p_r << ", " << temp.prime << ") ";
    }

    
    long long x;
    long **gaussian;
    long tonelliNumsSize = tonelliNums.size();
    gaussian = createMatrix(2*M, tonelliNumsSize);
    initialize(gaussian, 2*M, tonelliNumsSize);
    
    //gaussian.resize(quadPrimes.size(), vector<long>(2*M, 3)); //i is row, j is col
    sr = floor(sqrt(13592675504123));
    M=quadPrimes.size();
    
    for(int k=0; k < 2*M; k++) {
        for(int l=0; l < tonelliNums.size(); l++){
            x=power(sr-M+k+1, 2)-n;
            //sieve_temp = mod(x,tonelliNums[l].prime);
            if((x-tonelliNums[l].r)/tonelliNums[l].prime == 0) {
                gaussian[k][l] += floor(0.5+log(tonelliNums[l].prime));
            }
            if((x-tonelliNums[l].p_r)/tonelliNums[l].prime == 0) {
                gaussian[l][k] += floor(0.5+log(tonelliNums[l].prime));
            }
        }
    }
    
    //find rows with value > .5*log(n)+log(M) - TlogB
    long long temp_gaussian=0;
    vector<long long> trial_division;
    for(int m=0; m < 2*M; m++){
        for(int n=0; n < tonelliNumsSize; n++) {
            temp_gaussian += gaussian[m][n];
        }
        if(temp_gaussian >= (.5*log(n)+log(M)-1.5*log(tonelliNumsSize))){
            trial_division.push_back(power(sr-M+m+1, 2)-n);
        }
        temp_gaussian=0;
    }
    //print(trial_division, trial_division.size());
    
    //trial division
    long ** exponent = createMatrix(trial_division.size(), tonelliNumsSize);
    initialize(exponent, trial_division.size(), tonelliNumsSize);
    for(int v=0; v < trial_division.size(); v++) {
        for(int b=0; b < tonelliNumsSize && trial_division[v] != 0; b++){
            while(mod(trial_division[v], tonelliNums[b].prime) == 0){
                exponent[v][b]++;
                trial_division[v]=exponent[v][b]/tonelliNums[b].prime;
            }
        }
    }
    
    //look for 1's in columns 1...n
    for(int i=0; i < trial_division.size(); i++){
        for(int j=0; j < tonelliNumsSize; j++){
            if(exponent[i]j] > 0){
                cout << exponent[i]j] << " ";
                break;
            }
        }
    }
    
    
    return 0;
}













