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


int main(int argc, const char * argv[]) {
    
//    BigInteger n = 13592675504123;
    BigInteger n = 4999486012441;
    long long long_n = 4999486012441;
    BigInteger primes[18383];
    vector<BigInteger> quadPrimes; //primes that are quad res of my prime
    vector<tonelli_pair> tonelliNums; //tonelli numbers of factor base
   
    BigInteger M=10000;
    tonelli_pair temp;
    BigInteger sr;
    BigInteger length = BigInteger(2)*M;
    
    loadPrimes(primes);
    
    //find quad res n mod p
    BigInteger numOfPrimes=0;
    for(BigInteger i=1; i < 18383 && numOfPrimes < 28; i++){
        if(jacobi(n,primes[i.toInt()])==1){
            quadPrimes.push_back(primes[i.toInt()]);
            numOfPrimes++;
        }
    }
    
    //find tonelli pairs w.r.t. n mod p
    for(BigInteger i=0; i < quadPrimes.size(); i++) {
        temp = tonelli(n, quadPrimes[i.toInt()]);
        tonelliNums.push_back(temp);
        cout << endl;
    }

    print(tonelliNums, tonelliNums.size());
    BigInteger x;

    BigInteger **gaussian;
    BigInteger tonelliNumsSize = tonelliNums.size();
    gaussian = createMatrix(length.toInt(), tonelliNumsSize.toInt());
    initialize(gaussian, length.toInt(), tonelliNumsSize.toInt());
    
    sr = BigInteger((long)floor(sqrt(13592675504123)));
    M=quadPrimes.size();
    
    for(int k=0; BigInteger(k) < length; k++) {
        for(int l=0; BigInteger(l) < tonelliNumsSize; l++){
            
            x=mod(sr-M+k, tonelliNums[l].prime);
            
            if(x==mod(tonelliNums[l].r, tonelliNums[l].prime)) {
                gaussian[k][l] += (long)floor(0.5+log(tonelliNums[l].prime.toUnsignedLong()));
            }
            if(x==mod(tonelliNums[l].p_r, tonelliNums[l].prime)) {
                gaussian[k][l] += (long)floor(0.5+log(tonelliNums[l].prime.toUnsignedLong()));
            }
        }
    }

    print(gaussian, length, tonelliNumsSize);
    cout << endl << endl;
    //find rows with value > .5*log(n)+log(M) - TlogB
    BigInteger temp_gaussian=0;
    vector<BigInteger> trial_division;
    for(BigInteger m=0; m < length; m++){
        for(BigInteger n=0; n < tonelliNumsSize; n++) {
            temp_gaussian += gaussian[m.toInt()][n.toInt()];
        }
        if(temp_gaussian >= (long)(.5*log(long_n)+log(M.toUnsignedLong())-1.5*log(tonelliNumsSize.toInt()))){
            trial_division.push_back(sr-M+m);
        }
        temp_gaussian=0;
    }
    print(trial_division, trial_division.size());
    cout << endl << endl;
    //trial division
    BigInteger ** exponent = createMatrix(trial_division.size(), tonelliNumsSize.toInt());
    initialize(exponent, trial_division.size(), tonelliNumsSize.toInt());
    for(BigInteger v=0; v < trial_division.size(); v++) {
        for(BigInteger b=0; b < tonelliNumsSize && trial_division[v.toInt()] != 0; b++){
            while(mod(trial_division[v.toInt()], tonelliNums[b.toInt()].prime) == 0){
                exponent[v.toInt()][b.toInt()]++;
                trial_division[v.toInt()]=exponent[v.toInt()][b.toInt()]/tonelliNums[b.toInt()].prime.toUnsignedLong();
            }
        }
    }
    //print(trial_division, trial_division.size());
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













