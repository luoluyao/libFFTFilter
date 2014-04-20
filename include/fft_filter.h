#include <vector>

#include "kiss_fft.h"

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
	vector<float> filter_state_;
	vector<float> window_;

	vector<float> fft_time_domain_buffer_;

	bool kernel_defined_;
	vector<kiss_fft_cpx> kernel_time_domain_buffer_;
	vector<kiss_fft_cpx> kernel_freq_domain_buffer_;

	int buffer_selector_;
	vector< vector<kiss_fft_cpx> > signal_time_domain_buffer_;
	vector< vector<kiss_fft_cpx> > signal_freq_domain_buffer_;

	vector<kiss_fft_cpx> filtered_freq_domain_buffer_;

	kiss_fft_cfg forward_fft_;
	kiss_fft_cfg inverse_fft_;
};

