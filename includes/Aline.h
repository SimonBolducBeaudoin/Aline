#pragma once

#include "../../Interpolation/includes/Interpolation.h"
#include "../../Multi_array/Multi_array.h"
#include "../../FFTW_extra/includes/fftw_extra.h"

typedef unsigned int uint ;
typedef std::complex<double> complex_d;

typedef Multi_array<double,1,uint> double_1D ;
typedef Multi_array<complex_d,1,uint> complex_d_1D ;

typedef Multi_array<double,2,uint> double_2D ;
typedef Multi_array<complex_d,2,uint> complex_d_2D ;

typedef Multi_array<double,3,uint64_t> double_3D ;
typedef Multi_array<complex_d,3,uint64_t> complex_d_3D ;

/*
	Functions
*/

double_1D Aline_raw( const double_1D& spectrum ) ;
double_1D Aline_resampling( const double_1D& spectrum , const double_1D& new_abs );
double_1D Aline_corrected( const double_1D& spectrum , const double_1D& new_abs , const complex_d_1D& disp_corr );

double_2D Bscan_raw( const double_2D& spectrum );
double_2D Bscan_resampling( const double_2D& spectrum , const double_1D& new_abs);
// double_1D Bscan_corrected( const double_1D& spectrum , const double_1D& new_abs , const complex_d_1D& disp_corr );

/*
	Classes 
*/