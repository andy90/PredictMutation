#include <iostream>
#include <cstdlib>
#include <cmath>
#include <random>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List get_distribution(NumericMatrix J, NumericVector h){
    int nsite = J.ncol();

    NumericVector s(nsite), sav(nsite);
    NumericMatrix scoupleav(nsite, nsite);

    for (int i = 0; i < nsite; i++) //initialization
    {
        s[i] = 0;
        sav[i] = 0;
        for (int j = 0; j < nsite; j++)
        {
            scoupleav(i, j) = 0;
        }
    }

    int nsteps = 200000;
    int nsteps_dump = 10000;
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0, 1);
    for (int iter = 0; iter < nsteps; iter++)
    {
        int isite = std::rand() % nsite;
        int s_new = 1 - s[isite]; //flip the spin
        double deltaE = 0;
        deltaE = deltaE + h[isite] * (s_new - s[isite]);
        for (int i = 0; i < isite; i++)
        {
            if (i != isite)
            {
                deltaE = deltaE + J(i, isite) * s[i] * (s_new - s[isite]);
            }
        }

        if (deltaE < 0)
        {
            s[isite] = s_new;
        }
        else
        {
            double prob = exp(-deltaE);
            double p = dis(gen);
            if (p < prob)
            {
                s[isite] = s_new;
            }
        }
        
        if (iter >= nsteps_dump)
        {
            for (int i = 0; i < nsite; i++)
            {
                sav[i] += s[i];
            }
    
            for (int i = 0; i < nsite; i++)
            {
                for (int j = 0; j < nsite; j++)
                {
                    scoupleav(i, j) += s[i] * s[j];
                }
            }
        }
    }

    for (int i = 0; i < nsite; i++)
    {
        sav[i] = sav[i] / (nsteps - nsteps_dump);
    }

    for (int i = 0; i < nsite; i++)
    {
        for (int j = 0; j < nsite; j++)
        {
            scoupleav(i, j) = scoupleav(i, j) / (nsteps - nsteps_dump);
        }
    }

    List ret;
    ret["sav"] = sav;
    ret["scoupleav"] = scoupleav;

    return ret;
}
