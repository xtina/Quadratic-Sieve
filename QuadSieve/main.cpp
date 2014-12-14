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
#include "BigIntegerAlgorithms.hh"
#include "functions.h"
using namespace std;
const BigInteger n = 4999486012441;//13592675504123;
const long long long_n = 4999486012441;//13592675504123;
const BigInteger M=5000;
const BigInteger numOfPrimes=28;


int main(int argc, const char * argv[]) {

    BigInteger primes[18383];
    vector<BigInteger> quadPrimes; //primes that are quad res of my prime
    vector<tonelli_pair> tonelliNums; //tonelli numbers of factor base
   
    
    tonelli_pair temp;
    BigInteger sr;
    BigInteger length = BigInteger(2)*M;
    
    loadPrimes(primes);
    
    //find quad res n mod p
    BigInteger cnt=0;
    quadPrimes.push_back(2);
    for(BigInteger i=1; i < 18383 && cnt < numOfPrimes; i++){ //skip 2
        if(jacobi(n,primes[i.toInt()])==1){
            quadPrimes.push_back(primes[i.toInt()]);
            cnt++;
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
    BigInteger gassian_temp, prime_temp, neg;
    BigInteger **matrix;
    BigInteger tonelliNumsSize = tonelliNums.size();
    matrix = createMatrix(length.toInt(), tonelliNumsSize.toInt());
    initialize(matrix, length.toInt(), tonelliNumsSize.toInt());
    
    //begin sieving
    sr = BigInteger((unsigned long)floor(sqrt(long_n)));
    
    for(int k=0; BigInteger(k) < length; k++) {
        for(int l=0; BigInteger(l) < tonelliNumsSize; l++){
            
            gassian_temp=mod(sr-M+k, tonelliNums[l].prime);

            if(gassian_temp==tonelliNums[l].r || gassian_temp==tonelliNums[l].p_r) {
                matrix[k][l] += (long)floor(0.5+log(tonelliNums[l].prime.toUnsignedLong()));
            }
        }
    }

    //print(gaussian, length, tonelliNumsSize);
    //find rows with value > .5*log(n)+log(M) - TlogB
    BigInteger temp_gaussian=0;
    BigInteger negative = BigInteger(-1);
    unsigned long limit = .5*log(long_n)+log(M.toLong())-1.5*log(tonelliNums[tonelliNumsSize.toUnsignedLong()-1].prime.toUnsignedLong());
    vector<BigInteger> trial_division;
    for(BigInteger m=0; m < length; m++){
        for(BigInteger n=0; n < tonelliNumsSize; n++) {
            temp_gaussian += matrix[m.toInt()][n.toInt()];
        }
        if (temp_gaussian >= limit) {
            trial_division.push_back(sr-M+m);
            //trial_division.push_back(M*negative+m);
        }
        temp_gaussian=0;
    }
    //print(trial_division, trial_division.size()); cout << endl << endl;
    
    //trial division
    unsigned long rows = trial_division.size();
    unsigned long cols = tonelliNumsSize.toUnsignedLong();
    
    int ** exponent = intMatrix(rows, 2*cols); //twice as many columns for identity matrix
    initialize(exponent, rows, 2*cols); //initialize array to 0
    int **factored = intMatrix(rows, 2*cols);
    initialize(factored, rows, 2*cols);
    unsigned long counter=0;
    vector<BigInteger> fullyFactored;
    BigInteger r; unsigned long curI;
    for(int v=0; v < rows; v++) {
        curI = (trial_division[v]-sr+M).toUnsignedLong();
        r = trial_division[v]*trial_division[v]-n;
        for(int b=0; b < cols && r != 1 && r!= -1; b++){ //while r isn't fully factored...
            while(r%tonelliNums[b].prime == 0){
                exponent[v][b]++;
                r /= tonelliNums[b].prime;
            }
        }
        if(r == 1 || r == -1){
            factored[counter] = exponent[v];
            factored[counter][cols+counter] = 1; //adding identity matrix
            fullyFactored.push_back(trial_division[v]);
            cout << (trial_division[v]-sr+M).toUnsignedLong() << " ";
            counter++;
        }
    }

//    destroyMatrix(exponent, rows);
    for(unsigned long i=rows-counter; i < rows; i++){ //delete unnecessary rows
        delete [] factored[i];
    }
    
    rows=counter;
    //print(factored, counter, 2*cols);
    int ** row_reduced = intMatrix(rows, 2*cols); //copy matrix into another to be row reduced
    copyMatrix(row_reduced, factored, counter, 2*cols);
    cout << endl;
    to_reduced_row_echelon_form(row_reduced, counter, 2*cols);
    print(row_reduced, counter, 2*cols);
    cout << endl;
    
    //find zero row
    BigInteger x=1, y=1, result, zeroRow;
    int *ex = new int[cols];
    unsigned long startRow =0;
    initialize(ex, cols);
    while(true) {
        zeroRow = findZeroRow(row_reduced, rows, cols, startRow);
        if(zeroRow == -1) {
            return 0;
        }
        //calculate x
        for(unsigned long i=cols; i < 2*cols; i++) {
            if(row_reduced[zeroRow.toUnsignedLong()][i] == 1) {
                x *= fullyFactored[i-cols] % n;
            }
        }
        
        //calculate y
        for(unsigned long i=cols; i < 2*cols; i++) {
            if(row_reduced[zeroRow.toUnsignedLong()][i] == 1) {
                addRows(ex, factored[i-cols], cols);
            }
            else
                ex[i] = 0;
        }

        for(int i=0; i < cols; i++)
            y *= power(tonelliNums[i].prime, ex[i]/2);
        result=gcd(x-y, n);
        if(result != 1 && result != 0)
            break;
        startRow=zeroRow.toUnsignedLong()+1;
        if(zeroRow > rows) {
            cout << "factorization not found" << endl;
            break;
        }
    }
    cout << result.toUnsignedLong();
    
    return 0;
}













