#include <assert.h>
#include <cmath>

#include "fftfilter.h"
#include "kiss_fft.h"


FFTFilter::FFTFilter(int filter_len) :
		filter_len_(filter_len),
		filter_state_(filter_len, 0.0f),
		fft_time_domain_buffer_(filter_len * 2, 0.0f),
		forward_fft_(nullptr),
		inverse_fft_(nullptr){

	forward_fft_=kiss_fft_alloc(filter_len*2,0,NULL,NULL);
	inverse_fft_=kiss_fft_alloc(filter_len*2,0,NULL,NULL);

}

FFTFilter::~FFTFilter()  {

}

