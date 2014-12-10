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
#include "functions.h"
using namespace std;


int main(int argc, const char * argv[]) {
    
    long long n = 13592675504123;
    long primes[18383];
    list<long> quadPrimes; //primes that are quad res of my prime
    vector<tonelli_pair> tonelliNums; //tonelli numbers of factor base
   
    long M=10000;
    tonelli_pair temp;
    long long sr;
    
    loadPrimes(primes);
    
    //find quad res n mod p
    int numOfPrimes=0;
    for(int i=1; i < 18383 && numOfPrimes < 50; i++){
        if(jacobi(n,primes[i])==1){
            quadPrimes.push_back(primes[i]);
            numOfPrimes++;
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
    
    sr = floor(sqrt(13592675504123));
    M=quadPrimes.size();
    
    for(int k=0; k < 2*M; k++) {
        for(int l=0; l < tonelliNumsSize; l++){
            x=mod(sr-M+k, tonelliNums[l].prime);
            //sieve_temp = mod(x,tonelliNums[l].prime);
            if(x==mod(tonelliNums[l].r, tonelliNums[l].prime)) {
                gaussian[k][l] += floor(0.5+log(tonelliNums[l].prime));
            }
            if(x==mod(tonelliNums[l].p_r, tonelliNums[l].prime)) {
                gaussian[l][k] += floor(0.5+log(tonelliNums[l].prime));
            }
        }
    }
    
    print(gaussian, 2*M, tonelliNumsSize);
    
    //find rows with value > .5*log(n)+log(M) - TlogB
    long long temp_gaussian=0;
    vector<long long> trial_division;
    for(int m=0; m < 2*M; m++){
        for(int n=0; n < tonelliNumsSize; n++) {
            temp_gaussian += gaussian[m][n];
        }
        if(temp_gaussian >= (.5*log(n)+log(M)-1.5*log(tonelliNumsSize))){
            trial_division.push_back(sr-M+m);
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
    print(trial_division, trial_division.size());
    /*
    //look for 1's in columns 1...n
    for(int i=0; i < trial_division.size(); i++){
        for(int j=0; j < tonelliNumsSize; j++){
            if(exponent[i][j] > 0){
                cout << exponent[i][j] << " ";
                break;
            }
        }
    }
    */
    
    return 0;
}













