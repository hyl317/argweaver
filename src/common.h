#ifndef DLCOAL_COMMON_H
#define DLCOAL_COMMON_H

/*=============================================================================

    SPIDIR - SPecies Informed DIstanced-based Reconstruction
    
    Matt Rasmussen
    Wed Jun 13 22:09:24 EDT 2007


=============================================================================*/

// headers c++ 
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <algorithm>

using namespace std;

namespace spidir {

// constants
#ifndef INFINITY
#   define INFINITY 1e1000
#endif



//=============================================================================
// Math


inline float frand(float max=1.0)
{ return rand() / float(RAND_MAX) * max; }

inline float frand(float min, float max)
{ return min + (rand() / float(RAND_MAX) * (max-min)); }

inline int irand(int max)
{
    const int i = int(rand() / float(RAND_MAX) * max); 
    return (i == max) ? max - 1 : i;
}

inline int irand(int min, int max)
{
    const int i = min + int(rand() / float(RAND_MAX) * (max - min)); 
    return (i == max) ? max - 1 : i;
}

inline float expovariate(float lambda)
{ return -log(frand()) / lambda; }



// computes log(a + b) given log(a) and log(b)
inline double logadd(double lna, double lnb)
{
    double diff = lna - lnb;
    if (lna == 1.0)
        return lnb;
    if (lnb == 1.0)
        return lna;
    if (diff < 500.0)
        return log(exp(diff) + 1.0) + lnb;
    else
        return lna;
}


inline double logsum(double *vals, int nvals)
{
    const double SUM_LOG_THRESHOLD = -15;
    double maxval = vals[0];
    int maxi = 0;

    // find maxval
    for (int i=1; i<nvals; i++) {
        if (vals[i] > maxval) {
            maxval = vals[i];
            maxi = i;
        }
    }

    double expsum = 1.0;
    for (int i=0; i<nvals; i++)
        if (i != maxi && vals[i] - maxval > SUM_LOG_THRESHOLD)
            expsum += exp(vals[i] - maxval);
  
    return maxval + log(expsum);        
}


template <class T>
T ipow(T val, int expo)
{
    T result = 1.0;
    unsigned int e = expo;

    if ((int)e < 0) {
	e = -e;
	val = 1.0 / val;
    }

    while (true) {
	if (e & 1)
	    result *= val;
	if ((e >>= 1) == 0)
	    break;
	val *= val;
    }

    return result;
}


template <class T>
T **new_matrix(int nrows, int ncols)
{
    T** mat = new T* [nrows];
    for (int i=0; i<nrows; i++)
        mat[i] = new T [ncols];
    return mat;
}


template <class T>
void delete_matrix(T **mat, int nrows)
{
    for (int i=0; i<nrows; i++)
        delete [] mat[i];
    delete [] mat;
}


int choose(int n, int k);

double fchoose(int n, int k);



} // namespace spidir

#endif // DLCOAL_COMMON_H