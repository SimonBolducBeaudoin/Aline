#include "../includes/Aline_py.h"

/* Are built-in Multi_array */

template<uint8_t Dim>
void check_dim( np_double& X )
{
	py::buffer_info info = X.request();

	if (info.ndim != Dim)
	{
		throw std::runtime_error("Wrong dimensions !");
	}		
};

template<class np_TYPE_X, class np_TYPE_Y>
void check_same_shape_0(  np_TYPE_X& X ,  np_TYPE_Y& Y )
{
	py::buffer_info Xinfo = X.request();
	py::buffer_info Yinfo = Y.request();
	
	if ( Xinfo.shape[0] != Yinfo.shape[0] )
	{
		throw std::runtime_error(" Shape[0] discordance ");
	}
};

template<class np_TYPE_X, class np_TYPE_Y>
void check_same_shape_1( np_TYPE_X& X ,  np_TYPE_Y& Y )
{
	py::buffer_info Xinfo = X.request();
	py::buffer_info Yinfo = Y.request();
	
	if ( Xinfo.shape[1] != Yinfo.shape[1] )
	{
		throw std::runtime_error(" Shape[1] discordance ");
	}
};

template<class np_TYPE_X, class np_TYPE_Y>
void check_same_shape_10( np_TYPE_X& X ,  np_TYPE_Y& Y )
{
	py::buffer_info Xinfo = X.request();
	py::buffer_info Yinfo = Y.request();
	
	if ( Xinfo.shape[1] != Yinfo.shape[0] )
	{
		throw std::runtime_error(" Shape_10 discordance ");
	}
};


// Numpy compatible functions
np_double Aline_raw_py(  np_double& np_spectrum )
{
	double_1D spectrum = double_1D::numpy(np_spectrum);
	double_1D aline = Aline_raw( spectrum ) ;
	
	return aline.move_py(spectrum.get_n_i()/2+1); // Returns only the valid part to numpy
};

np_double Aline_resampling_py(  np_double& np_spectrum ,  np_double& np_new_abs )
{
	check_same_shape_0( np_spectrum , np_new_abs );
	
	double_1D spectrum = double_1D::numpy(np_spectrum) ;
	double_1D new_abs = double_1D::numpy(np_new_abs) ;
	
	double_1D aline = Aline_resampling( spectrum , new_abs ) ;
	
	return aline.move_py(spectrum.get_n_i()/2+1); // Returns only the valid part to numpy
};

np_double Aline_corrected_py(  np_double& np_spectrum ,  np_double& np_new_abs ,  np_complex_d& np_disp_corr )
{
	check_same_shape_0( np_spectrum , np_new_abs );
	check_same_shape_0( np_spectrum , np_disp_corr );
	
	double_1D spectrum = double_1D::numpy(np_spectrum) ;
	double_1D new_abs = double_1D::numpy(np_new_abs) ;
	complex_d_1D disp_corr = complex_d_1D::numpy(np_disp_corr) ;
		
	double_1D aline = Aline_corrected( spectrum , new_abs , disp_corr );
	
	return aline.move_py(spectrum.get_n_i()/2+1); // Returns only the valid part to numpy
};

np_double Bscan_raw_py( np_double& np_spectrum )
{	
	double_2D spectrum = double_2D::numpy(np_spectrum);
	double_2D bscan = Bscan_raw( spectrum ) ;
	
	return bscan.move_py( spectrum.get_n_j() , spectrum.get_n_i()/2+1 ); // Returns only the valid part to numpy
};

np_double Bscan_resampling_py( np_double& np_spectrum , np_double& np_new_abs )
{	
	check_same_shape_10( np_spectrum , np_new_abs );
	double_2D spectrum = double_2D::numpy(np_spectrum);
	double_1D new_abs = double_1D::numpy(np_new_abs) ;
	
	double_2D bscan = Bscan_resampling( spectrum , new_abs ) ;
	
	return bscan.move_py( spectrum.get_n_j() , spectrum.get_n_i()/2+1 ); // Returns only the valid part to numpy
};

void init_Aline(py::module &m)
{
	m.def("Aline_raw", &Aline_raw_py , "spectrum"_a.noconvert() ) ;
	m.def("Aline_resampling", &Aline_resampling_py , "spectrum"_a.noconvert() , "new_abs"_a.noconvert() ) ;
	m.def("Aline_corrected", &Aline_corrected_py , "spectrum"_a.noconvert() , "new_abs"_a.noconvert() , "disp_corr"_a.noconvert() ) ;
	m.def("Bscan_raw", &Bscan_raw_py , "spectrum"_a.noconvert() ) ;
	m.def("Bscan_resampling", &Bscan_resampling_py , "spectrum"_a.noconvert() , "new_abs"_a.noconvert() ) ;
}


