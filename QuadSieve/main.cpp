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
#include "BigInteger.hh"
#include "functions.h"
using namespace std;
const BigInteger n = 87463;
const long long long_n = 87463;

// Swap rows i and k of a matrix A
// Note that due to the reference, both dimensions are preserved for
// built-in arrays
void swap_rows(BigInteger **A,
               BigInteger i,
               BigInteger k, BigInteger size) {
    
    // check indices
    BigInteger temp;
    for (BigInteger col = 0; col <= size; ++col) {
        std::swap(A[i.toUnsignedLong()][col.toUnsignedLong()], A[k.toUnsignedLong()][col.toUnsignedLong()]);
    }
}

// divide row i of matrix A by v
void divide_row(BigInteger ** A,
                BigInteger i,
                BigInteger v, BigInteger size) {
    
    for (BigInteger col = 0; col <= size; ++col)
        A[i.toUnsignedLong()][col.toUnsignedLong()] /= v;
}

// in matrix A, add v times row k to row i
void add_multiple_row(BigInteger **A,
                      BigInteger i,
                      BigInteger k,
                      BigInteger v, BigInteger size) {
    
    for (BigInteger col = 0; col <= size; ++col)
        A[i.toUnsignedLong()][col.toUnsignedLong()] += v * A[k.toUnsignedLong()][col.toUnsignedLong()];
}

// convert A to reduced row echelon form
void to_reduced_row_echelon_form(BigInteger **A, BigInteger rows, BigInteger cols) {
    
    BigInteger lead = 0;
    
    for (BigInteger j = 0;j <= rows; ++j) {
        if (lead > cols)
            return;
        BigInteger i = j;
        while (A[i.toLong()][lead.toLong()] == 0) {
            ++i;
            if (i > rows) {
                i = j;
                ++lead;
                if (lead > cols)
                    return;
            }
        }
        swap_rows(A, i, j, cols);
        divide_row(A, j, A[j.toLong()][lead.toLong()], cols); cout << endl;
        for (BigInteger i = 0; i < rows; ++i) {
            cout << "lead: " << lead.toInt() << " " << A[i.toLong()][lead.toLong()].toLong()*(-1) << endl;
            if (i != j) {
                add_multiple_row(A, i, j, A[i.toLong()][lead.toLong()].negative, cols);
                cout << "after: "; print(A[i.toInt()], cols.toInt()); cout << endl;;
            }
        }
        cout << endl; print(A[j.toInt()], cols.toInt()); cout << endl;
    }
}


int main(int argc, const char * argv[]) {
    
//    BigInteger n = 13592675504123;

    BigInteger primes[18383];
    vector<BigInteger> quadPrimes; //primes that are quad res of my prime
    vector<tonelli_pair> tonelliNums; //tonelli numbers of factor base
   
    BigInteger M=30;
    tonelli_pair temp;
    BigInteger sr;
    BigInteger length = BigInteger(2)*M;
    
    loadPrimes(primes);
    
    //find quad res n mod p
    BigInteger numOfPrimes=0;
    quadPrimes.push_back(2);
    for(BigInteger i=1; i < 18383 && numOfPrimes < 5; i++){ //skip 2
        if(jacobi(n,primes[i.toInt()])==1){
            quadPrimes.push_back(primes[i.toInt()]);
            numOfPrimes++;
        }
    }
    //print(quadPrimes, quadPrimes.size());
    temp.r = 1; temp.p_r = 1; temp.prime = 2;
    //find tonelli pairs w.r.t. n mod p
    tonelliNums.push_back(temp);
    for(BigInteger i=1; i < quadPrimes.size(); i++) {
        temp = tonelli(n, quadPrimes[i.toInt()]);
        tonelliNums.push_back(temp);
    }

    //print(tonelliNums, tonelliNums.size());
    BigInteger x;

    BigInteger **gaussian;
    BigInteger tonelliNumsSize = tonelliNums.size();
    gaussian = createMatrix(length.toInt(), tonelliNumsSize.toInt());
    initialize(gaussian, length.toInt(), tonelliNumsSize.toInt());
    
    sr = BigInteger((unsigned long)floor(sqrt(long_n)));
    //cout << sr.toLong() << endl << endl;
    
    for(int k=0; BigInteger(k) < length; k++) {
        for(int l=0; BigInteger(l) < tonelliNumsSize; l++){
            
            x=mod(sr-M+k, tonelliNums[l].prime);
            
            if(x==mod(tonelliNums[l].r, tonelliNums[l].prime)) {
                gaussian[k][l] += (long)floor(0.5+log(tonelliNums[l].prime.toUnsignedLong()));
            }
            else if(x==mod(tonelliNums[l].p_r, tonelliNums[l].prime)) {
                gaussian[k][l] += (long)floor(0.5+log(tonelliNums[l].prime.toUnsignedLong()));
            }
        }
    }

    //print(gaussian, length, tonelliNumsSize);
    //find rows with value > .5*log(n)+log(M) - TlogB
    BigInteger temp_gaussian=0;
    BigInteger negative = BigInteger(-1);
    unsigned long limit = .5*log(long_n)+log(M.toLong())-1.5*log(tonelliNums[tonelliNumsSize.toInt()-1].prime.toInt());
    vector<BigInteger> trial_division;
    for(BigInteger m=0; m < length; m++){
        for(BigInteger n=0; n < tonelliNumsSize; n++) {
            temp_gaussian += gaussian[m.toInt()][n.toInt()];
        }
        if (temp_gaussian >= limit) {
            
            trial_division.push_back(sr-M+m);
            //trial_division.push_back(M*negative+m);
        }
        temp_gaussian=0;
    }
    //print(trial_division, trial_division.size()); cout << endl << endl;
    
    //trial division
    BigInteger ** exponent = createMatrix(trial_division.size(), tonelliNumsSize.toInt());
    initialize(exponent, trial_division.size(), tonelliNumsSize.toInt());
    BigInteger **factorMe = createMatrix(trial_division.size(), tonelliNumsSize.toInt());
    int counter=0;
    BigInteger r, lala;
    for(BigInteger v=0; v < trial_division.size(); v++) {
        lala = trial_division[v.toInt()];
        r = trial_division[v.toInt()]*trial_division[v.toInt()]-n;
        for(BigInteger b=0; b < tonelliNumsSize && r != 1 && r!= -1; b++){
            while(r%tonelliNums[b.toInt()].prime == 0){
                if(exponent[v.toInt()][b.toInt()] == 1)
                    exponent[v.toInt()][b.toInt()]=0;
                else
                    exponent[v.toInt()][b.toInt()] = 1;
                r /= tonelliNums[b.toInt()].prime;
            }
        }
        if(r == 1 || r == -1){
            //cout << lala.toUnsignedLong() << " " << (lala-sr+M).toLong() << endl;
            factorMe[counter] = exponent[v.toInt()];
            counter++;//cout << (lala-M).toLong() << " ";
        }
    }
//    cout << counter << endl;
    for(unsigned long i=trial_division.size()-counter; i < trial_division.size(); i ++){
        delete [] factorMe[i];
    }
    print(factorMe, counter, tonelliNumsSize);
    cout << endl << endl;
    to_reduced_row_echelon_form(factorMe, counter, tonelliNumsSize.toInt());
    print(factorMe, counter, tonelliNumsSize);
    /*
    //look for 1's in columns 1...n
    for(BigInteger i=0; i < trial_division.size(); i++){
        for(BigInteger j=0; j < tonelliNumsSize; j++){
            if(exponent[i][j] > 0){
                cout << exponent[i][j] << " ";
                break;
            }
        }
    }
*/
    
    return 0;
}













