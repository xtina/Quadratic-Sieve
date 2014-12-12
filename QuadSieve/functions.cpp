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
    return (a%b+b)%b;
}

//modular exponentiation
BigInteger modpow(BigInteger base, BigInteger exponent, BigInteger modulus) {
    /*
    BigInteger result = 1;
    
    while (exponent > 0) {
        if ((exponent & 1) == 1) {
            // multiply in this bit's contribution while using modulus to keep result small
            result = (result * base) % modulus;
        }
        // move to the next bit of the exponent, square (and mod) the base accordingly
        exponent >>= 1;
        base = (base * base) % modulus;
    }
    
    return result;*/
    
    BigInteger b = base;
    BigInteger n = 1;
    while (exponent != 0){
    if(mod(exponent, 2) == 1){
        n = mod(n*base, modulus);
    }
        exponent = (unsigned long)floor((exponent/BigInteger(2)).toUnsignedLong());
        base = (base*base)%modulus;
    }
    
    return n;
    
}

//raise i to the power j
BigInteger power(BigInteger i, BigInteger j){
    BigInteger temp=1;
    for(BigInteger k=0; k<j; k++){
        temp*=i;
    }
    return temp;
}

//jacobi symbol calculator
BigInteger jacobi(BigInteger a, BigInteger m){
    a = mod(a,m);
    BigInteger t = 1;
    BigInteger c;
    BigInteger temp;
    
    while(a != 0){
        c = 0;
        while(a%2 == 0){
            a /= 2;
            c=BigInteger(1)-c; //if c is 1, exponent is odd
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
tonelli_pair tonelli(BigInteger a, BigInteger p){
    BigInteger b=0;
    BigInteger t = p-1;
    BigInteger s=0;
    BigInteger temp;
    tonelli_pair pair;
    
    if(jacobi(a, p)  == -1) {
        cout << "uh oh" << endl;
    }
    
    //find quad non res
    temp = 0;
    while(temp == 1 || temp == 0){
        b++;
        temp = jacobi(b, p);
    }
    
    //pull out factors of 2
    while(mod(t, 2) == 0){
        t /= 2;
        s++;
    }
    /*
    if(s==1){
        pair.r = modpow(a, (p+1)/4, p);
        pair.p_r = p-modpow(a, (p+1)/4, p);
        while (pair.p_r < 0) {
            pair.p_r += p;
        }
        pair.prime = p;
        return pair;
    }
    */
    BigInteger i = 2;
    BigInteger c = mod(a*b*b, p);
    for(BigInteger k=1; k <= s-1; k++){
        if(modpow(c, power(2, s-k-1)*t, p) == p-1){
            i += power(2, k);
            c = mod(c*power(b, power(2, k)), p);
        }
    }

    cout << modpow(b, (i*t)/2, p).toInt() << " " << modpow(a, (t+1)/2, p).toInt() << " " << p.toUnsignedLong() << " ";
    pair.r = mod(modpow(b, (i*t)/2, p)*modpow(a, (t+1)/2, p), p);
    pair.p_r = p-pair.r;
    cout << pair.r.toUnsignedLong() << " ";
    while (pair.p_r < 0) {
        pair.p_r += p;
    }
    
    pair.prime = p;
    
    return pair;
}
/*

tonelli_pair tonelli(BigInteger a, BigInteger p){
    BigInteger s=0;
    BigInteger i=2;
    BigInteger temp_p;
    BigInteger quad=1;
    temp_p = p-1;
    tonelli_pair pair;
    //pull out factors of 2
    while(temp_p%2 == 0){
        temp_p=temp_p/2;
        s++;
    }
    //Q is the remainder
    BigInteger Q = temp_p;
    
    if(s == 1){
        pair.r = modpow(a, (p+1)/4, p);
        pair.p_r = p - pair.r;
        pair.prime = p;
        return pair;
    }
    
    //find b
    BigInteger z=0;
    BigInteger k=1; //counter
    while(quad==1 || quad==0){
        quad=jacobi(k, p);
        z=k;
        k++;
    }
    
    BigInteger c = power(z, Q);
    BigInteger R = power(a, (Q+1)/2);
    BigInteger t = power(a, Q);
    BigInteger M = s;
    BigInteger b;
    
    while(true){
        
        if(t%p == 1){
            pair.r = R;
            pair.p_r = p-R;
            pair.prime = p;
            return pair;
        }
        
        while(modpow(t, power(2,i) ,p) != 1 && i < M){
            i++;
        }
        b = modpow(c, power(2,M-i-1), p);
        R = (R*b)%p;
        t = (t * b* b)%p;
        c = (b*b)%p;
        M = i;
    }
    
}
*/
void loadPrimes(BigInteger *primes){
    ifstream infile;
    unsigned long temp;
    int i=0;
    infile.open("/Users/Christina/Documents/QuadSieve/QuadSieve/primes.txt");
    while (infile >> temp) {
        
        primes[i++]= BigInteger(temp);
    }
    infile.close();
}

BigInteger ** createMatrix(unsigned long rows, unsigned long cols) {
    BigInteger **mat;
    mat = new BigInteger*[rows];
    for(BigInteger i = 0; i < rows; ++i) {
        mat[i.toInt()] = new BigInteger[cols];
    }
    return mat;
}

void initialize(BigInteger ** mat, unsigned long rows, unsigned long cols) {
    for(BigInteger m=0; m < rows; m++)
        for(BigInteger n=0; n < cols; n++) {
            mat[m.toInt()][n.toInt()] = BigInteger(1);
        }
}

void destroyMatrix(BigInteger **matrix, BigInteger size) {
    for(BigInteger i=0; i < size; i++)
        delete [] matrix[i.toInt()];
    delete [] matrix;
}

void printToFile(BigInteger ** matrix, BigInteger rows, BigInteger cols){
    ofstream out;
    out.open("trial-division-matrix.txt");
    for(BigInteger i=0; i < rows; i++){
        for(BigInteger j=0; j < cols; j++){
            out << matrix[i.toInt()][j.toInt()].toUnsignedLong() << " ";
        }
        out << endl;
    }
    out.close();
}

void print(BigInteger ** matrix, BigInteger rows, BigInteger cols){
    for(BigInteger i=0; i < rows; i++){
        for(BigInteger j=0; j < cols; j++){
            cout << matrix[i.toInt()][j.toInt()].toUnsignedLong() << " ";
        }
        cout << endl;
    }
}

void print(vector<tonelli_pair> pair, BigInteger size){
    for(BigInteger i=0; i < size; i++){
        cout << pair[i.toInt()].r.toLong() << " " << pair[i.toInt()].p_r.toUnsignedLong() << " " << pair[i.toInt()].prime.toUnsignedLong() << endl;
    }
}

void print(vector<BigInteger> matrix, BigInteger rows) {
    for(BigInteger i=0; i < rows; i++) {
        cout << matrix[i.toInt()].toUnsignedLong() << " ";
    }
}

