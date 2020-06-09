#include "../includes/Aline.h"

/* BEGIN Vectorization land */

inline void absolute_SIMD_loop (uint n_i, double* out, complex_d* in )
{
	#pragma GCC ivdep
	for( int i=0; i < n_i ; i++ )
	{
		out[i] = std::abs(in[i]);
	}
}

/* END  Vectorization land */

double_1D Aline_raw( const double_1D& spectrum )
{
	int l_fft = (int)spectrum.get_n_i();
	double_1D aline ( l_fft + 2 , fftw_malloc , fftw_free ) ; 
	complex_d_1D aline_complex( (complex_d*)aline.get_ptr() , l_fft/2 + 1 ); 
	
	fftw_import_wisdom_from_filename("FFTW_Wisdom.dat");
	fftw_plan plan = fftw_plan_dft_r2c_1d( l_fft , &spectrum[0] , reinterpret_cast<fftw_complex*>(&aline_complex[0]) , FFTW_EXHAUSTIVE );
	fftw_export_wisdom_to_filename("FFTW_Wisdom.dat");	
	
	fftw_execute( plan );
	
	for( int i=0; i < l_fft/2+1 ; i++ )
	{
		aline[i] = std::abs(aline_complex[i]);
	}
	return aline; // Seul les l_fft/2 + 1 premier éléments sont pertinents
};

double_1D Aline_resampling( const double_1D& spectrum , const double_1D& new_abs )
{
	double_1D resamp_spectrum =  lin_interpol_1D( spectrum , new_abs ); // Calls move constructor ;
	return Aline_raw( resamp_spectrum ); // Calls move constructor ;	
};

double_1D Aline_corrected( const double_1D& spectrum , const double_1D& new_abs , const complex_d_1D& disp_corr )
{
	uint l_spectrum = spectrum.get_n_i() ; 
	complex_d_1D anal = analytic( spectrum ) ; // Could be change to avoid extra memory allocation
	
	double_1D spectrum_corr( l_spectrum ); 
	
	for( int i=0; i < (int)l_spectrum ; i++ )
	{
		spectrum_corr[i] = std::real( anal[i] * disp_corr[i] ) ; // This could be changed for in-place calculations
	}
	
	return Aline_resampling( spectrum_corr , new_abs ) ;
};

double_2D Bscan_raw( const double_2D& spectrum )
{
	int l_fft = (int)spectrum.get_n_i();
	uint n_bs =  spectrum.get_n_j();
	
	double_2D bscan ( l_fft + 2 , n_bs , fftw_malloc , fftw_free ) ; 
	complex_d_2D bscan_complex( (complex_d*)bscan.get_ptr() , l_fft/2 + 1 , n_bs ); 
	
	fftw_import_wisdom_from_filename("FFTW_Wisdom.dat");
	int n[] = {l_fft} ;
	int howmany =  n_bs ;
	fftw_plan plan = 
		fftw_plan_many_dft_r2c // In place ifft with strides 
		( 
			1 , // rank
			n , //  list of dimensions 
			howmany , // howmany (to do many ffts on the same core)
			spectrum[0] , // input
			NULL , // inembed
			1 , // istride
			0 , // idist
			reinterpret_cast<fftw_complex*>(bscan_complex[0]) , // output pointer
			NULL , //  onembed
			1 , // ostride
			0 , // odist
			FFTW_EXHAUSTIVE
		);
	fftw_export_wisdom_to_filename("FFTW_Wisdom.dat");	
	
	fftw_execute( plan );
	
	for ( uint j=0 ; j < n_bs ; j++ )
	{
		for( int i=0; i < l_fft/2+1 ; i++ )
		{
			bscan(j,i) = std::abs(bscan_complex(j,i));
		}
	}
	return bscan; // Seul les l_fft/2 + 1 premier éléments sont pertinents
};

double_2D Bscan_resampling( const double_2D& spectrum , const double_1D& new_abs )
{
	double_2D resamp_spectrum =  lin_interpol_1D( spectrum , new_abs ); // Calls move constructor ;

	int l_fft = (int)spectrum.get_n_i();
	uint n_bs =  spectrum.get_n_j();
	
	double_2D bscan ( l_fft + 2 , n_bs , fftw_malloc , fftw_free ) ; 
	complex_d_2D bscan_complex( (complex_d*)bscan.get_ptr() , l_fft/2 + 1 , n_bs ); 
	
	fftw_import_wisdom_from_filename("FFTW_Wisdom.dat");
	int n[] = {l_fft} ;
	int howmany =  n_bs ;
	fftw_plan plan = 
		fftw_plan_many_dft_r2c // In place ifft with strides 
		( 
			1 , // rank
			n , //  list of dimensions 
			howmany , // howmany (to do many ffts on the same core)
			spectrum[0] , // input
			NULL , // inembed
			1 , // istride
			0 , // idist
			reinterpret_cast<fftw_complex*>(bscan_complex[0]) , // output pointer
			NULL , //  onembed
			1 , // ostride
			0 , // odist
			FFTW_EXHAUSTIVE
		);
	fftw_export_wisdom_to_filename("FFTW_Wisdom.dat");	
	
	fftw_execute( plan );
	
	for ( uint j=0 ; j < n_bs ; j++ )
	{
		for( int i=0; i < l_fft/2+1 ; i++ )
		{
			bscan(j,i) = std::abs(bscan_complex(j,i));
		}
	}
	return bscan; // Seul les l_fft/2 + 1 premier éléments sont pertinents
	/**/
};

double_2D Bscan_corrected( const double_2D& spectrum , const double_1D& new_abs )
{
	double_2D resamp_spectrum =  lin_interpol_1D( spectrum , new_abs ); // Calls move constructor ;

	/**/
	int l_fft = (int)spectrum.get_n_i();
	uint n_bs =  spectrum.get_n_j();
	
	double_2D bscan ( l_fft + 2 , n_bs , fftw_malloc , fftw_free ) ; 
	complex_d_2D bscan_complex( (complex_d*)bscan.get_ptr() , l_fft/2 + 1 , n_bs ); 
	
	fftw_import_wisdom_from_filename("FFTW_Wisdom.dat");
	int n[] = {l_fft} ;
	int howmany =  n_bs ;
	fftw_plan plan = 
		fftw_plan_many_dft_r2c // In place ifft with strides 
		( 
			1 , // rank
			n , //  list of dimensions 
			howmany , // howmany (to do many ffts on the same core)
			spectrum[0] , // input
			NULL , // inembed
			1 , // istride
			0 , // idist
			reinterpret_cast<fftw_complex*>(bscan_complex[0]) , // output pointer
			NULL , //  onembed
			1 , // ostride
			0 , // odist
			FFTW_EXHAUSTIVE
		);
	fftw_export_wisdom_to_filename("FFTW_Wisdom.dat");	
	
	fftw_execute( plan );
	
	for ( uint j=0 ; j < n_bs ; j++ )
	{
		for( int i=0; i < l_fft/2+1 ; i++ )
		{
			bscan(j,i) = std::abs(bscan_complex(j,i));
		}
	}
	return bscan; // Seul les l_fft/2 + 1 premier éléments sont pertinents
	/**/
};





