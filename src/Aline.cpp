#include "../includes/Aline.h"

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
	double_1D resamp_spectrum =  interpolation_linear_index( spectrum , new_abs ); // Calls move constructor ;
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








