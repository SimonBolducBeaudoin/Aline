#pragma once

#include "Aline.h"

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
typedef py::array_t<double,py::array::c_style> np_double ; 
typedef py::array_t<complex_d,py::array::c_style> np_complex_d ; 


/*
No checks are made for now
*/

// Numpy compatible functions
np_double Aline_raw_py( np_double& np_spectrum );
np_double Aline_resampling_py( np_double& np_spectrum , np_double& np_new_abs );
np_double Aline_corrected_py( np_double& np_spectrum , np_double& np_new_abs , np_complex_d& np_disp_corr );

np_double Bscan_raw( np_double& spectrum );

void init_Aline (py::module &m);
