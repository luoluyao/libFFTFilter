#ifndef FFT_FILTER_H_
#define FFT_FILTER_H_

#include <vector>

// Real->Complex version of kiss_fft
#include "kiss_fftr.h"

using std::vector;

class FFTFilter {
public:
	FFTFilter(int filter_len);
	virtual ~FFTFilter();

	void SetKernel(const vector<float>& kernel);

	void AddSignalBlock(const vector<float>& signal_block);

	void GetResult(vector<float>* signal_block);


private:
	void Init();

	int filter_len_;
	int fft_len_;

	vector<float> filter_state_;
	vector<float> window_;

	bool kernel_defined_;
	vector<kiss_fft_scalar> kernel_time_domain_buffer_;
	vector<kiss_fft_cpx> kernel_freq_domain_buffer_;

	int buffer_selector_;
	vector< vector<kiss_fft_scalar> > signal_time_domain_buffer_;
	vector< vector<kiss_fft_cpx> > signal_freq_domain_buffer_;

	vector<kiss_fft_cpx> filtered_freq_domain_buffer_;

	kiss_fftr_cfg forward_fft_;
	kiss_fftr_cfg inverse_fft_;
};

#endif
