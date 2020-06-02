#pragma once
#include "fftw3.h"
#include <math.h>
inline void iexp(fftw_complex &res, fftw_complex inp)
{
	res[0] = exp(inp[0]) * cos(inp[1]);
	res[1] = exp(inp[0]) * sin(inp[1]);
}

inline void idiv(fftw_complex &n, double d)
{
	n[0] /= d;
	n[1] /= d;
}

inline void imul(fftw_complex &n, double d)
{
	n[0] *= d;
	n[1] *= d;
}

inline void imul(fftw_complex &res, fftw_complex l, fftw_complex r)
{
	fftw_complex tmp;
	tmp[0] = l[0] * r[0] - l[1] * r[1];
	tmp[1] = l[0] * r[1] + l[1] * r[0];
	res[0] = tmp[0]; 
	res[1] = tmp[1];
}

inline void imul_p(fftw_complex& res, fftw_complex l, fftw_complex r)
{
	res[0] = l[0] * r[0];
	res[1] = l[1] * r[1];
}

inline double iangle(fftw_complex inp)
{
	double res = atan2(inp[1], inp[0]);
	return res;
}


