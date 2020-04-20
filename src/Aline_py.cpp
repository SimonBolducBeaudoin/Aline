#include "../includes/Aline_py.h"

// Numpy compatible functions
np_double Aline_raw_py( const np_double& np_spectrum )
{
	/* checks ?? */
	double_1D spectrum = double_1D::numpy(np_spectrum);
	
	double_1D aline = Aline_raw( spectrum ) ;
	
	return aline.move_py(spectrum.get_n_i()/2+1); // Returns only the valid part to numpy
};

np_double Aline_resampling_py( const np_double& np_spectrum , const np_double& np_new_abs )
{
	/* more checks ?? */
	
	double_1D spectrum = double_1D::numpy(np_spectrum) ;
	double_1D new_abs = double_1D::numpy(np_new_abs) ;
	
	if ( spectrum.get_n_i() != new_abs.get_n_i() )
	{
		printf("line number %d in file %s\n", __LINE__, __FILE__);
		throw std::runtime_error(" spectrum.get_n_i() != new_abs.get_n_i() ");
	}
	
	double_1D aline = Aline_resampling( spectrum , new_abs ) ;
	
	return aline.move_py(spectrum.get_n_i()/2+1); // Returns only the valid part to numpy
};

np_double Aline_corrected_py( const np_double& np_spectrum , const np_double& np_new_abs , const np_complex_d& np_disp_corr )
{
	double_1D spectrum = double_1D::numpy(np_spectrum) ;
	double_1D new_abs = double_1D::numpy(np_new_abs) ;
	complex_d_1D disp_corr = complex_d_1D::numpy(np_disp_corr) ;
	
	if ( spectrum.get_n_i() != new_abs.get_n_i() )
	{
		printf("line number %d in file %s\n", __LINE__, __FILE__);
		throw std::runtime_error(" spectrum.get_n_i() != new_abs.get_n_i() ");
	}
	if ( spectrum.get_n_i() != disp_corr.get_n_i() )
	{
		printf("line number %d in file %s\n", __LINE__, __FILE__);
		throw std::runtime_error(" spectrum.get_n_i() != disp_corr.get_n_i() ");
	}
	
	double_1D aline = Aline_corrected( spectrum , new_abs , disp_corr );
	
	return aline.move_py(spectrum.get_n_i()/2+1); // Returns only the valid part to numpy
};



void init_Aline(py::module &m)
{
	m.def("Aline_raw", &Aline_raw_py , "spectrum"_a.noconvert() ) ;
	m.def("Aline_resampling", &Aline_resampling_py , "spectrum"_a.noconvert() , "new_abs"_a.noconvert() ) ;
	m.def("Aline_corrected", &Aline_corrected_py , "spectrum"_a.noconvert() , "new_abs"_a.noconvert() , "disp_corr"_a.noconvert() ) ;
}


