#pragma once

#include "../../includes/header_common.h"
#include "../../Multi_array/Multi_array.h"
// #include "../../Omp_extra/includes/omp_extra.h"
#include "../../Interpolation/includes/Interpolation.h"
#include "../../FFTW_extra/includes/fftw_extra.h"

typedef unsigned int uint ;
typedef std::complex<double> complex_d;

typedef Multi_array<double,1,uint> double_1D ;
typedef Multi_array<complex_d,1,uint> complex_d_1D ;

/*
	Functions
*/

double_1D Aline_raw( const double_1D& spectrum );
double_1D Aline_resampling( const double_1D& spectrum , const double_1D& new_abs);
double_1D Aline_corrected( const double_1D& spectrum , const double_1D& new_abs , const complex_d_1D& disp_corr );

/*
	Classes 
*/