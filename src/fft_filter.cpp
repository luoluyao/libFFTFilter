#include <assert.h>
#include <cmath>

#include "fft_filter.h"

using namespace std;

FFTFilter::FFTFilter(int filter_len) :
		filter_len_(filter_len),
		fft_len_(filter_len * 2),
		filter_state_(filter_len, 0.0f),
		kernel_defined_(false),
		kernel_time_domain_buffer_(fft_len_),
		kernel_freq_domain_buffer_(fft_len_ / 2 + 1),
		buffer_selector_(0),
		signal_time_domain_buffer_(2, vector<kiss_fft_scalar>(fft_len_)),
		signal_freq_domain_buffer_(2, vector<kiss_fft_cpx>(fft_len_ / 2 + 1)),
		filtered_freq_domain_buffer_(fft_len_ / 2 + 1) {
	bool is_power_of_two = ((fft_len_ != 0) && !(fft_len_ & (fft_len_ - 1)));
	assert(is_power_of_two && "Filter length must be a power of 2");

	forward_fft_ = kiss_fftr_alloc(fft_len_, 0, 0, 0);
	inverse_fft_ = kiss_fftr_alloc(fft_len_, 1, 0, 0);

	Init();
}

FFTFilter::~FFTFilter() {
	kiss_fft_free(forward_fft_);
	kiss_fft_free(inverse_fft_);
}

void FFTFilter::Init() {
	// Initialize all buffers with zeros.
	memset(&kernel_time_domain_buffer_[0], 0, sizeof(kiss_fft_scalar) * fft_len_);
	memset(&kernel_freq_domain_buffer_[0], 0, sizeof(kiss_fft_cpx) * (fft_len_ / 2 + 1));
	for (int i = 0; i < 2; ++i) {
		memset(&signal_time_domain_buffer_[i][0], 0,
				sizeof(kiss_fft_scalar) * fft_len_);
		memset(&signal_freq_domain_buffer_[i][0], 0,
				sizeof(kiss_fft_cpx) * (fft_len_ / 2 + 1));
	}
}

void FFTFilter::SetKernel(const std::vector<float>& kernel) {
	assert(kernel.size()==filter_len_ && "Kernel size must match filter length");

	for (int i = 0; i < filter_len_; ++i) {
		kernel_time_domain_buffer_[i] = kernel[i];
	}
	// Zero padding
	for (int i = filter_len_; i < fft_len_; ++i) {
		kernel_time_domain_buffer_[i] = 0.0f;
	}

	// Perform forward FFT transform
	kiss_fftr(forward_fft_, &kernel_time_domain_buffer_[0],
			&kernel_freq_domain_buffer_[0]);

	kernel_defined_ = true;
}

void FFTFilter::AddSignalBlock(const vector<float>& signal_block) {
	assert(signal_block.size()==filter_len_ && "Signal block size must match filter length");
	assert(kernel_defined_ && "No suitable kernel defined");

	// Switch buffer selector
	buffer_selector_ = !buffer_selector_;

	vector<kiss_fft_scalar>& time_domain_buffer = signal_time_domain_buffer_[buffer_selector_];
	vector<kiss_fft_cpx>& freq_domain_buffer = signal_freq_domain_buffer_[buffer_selector_];


	for (int i = 0; i < filter_len_; ++i) {
		time_domain_buffer[i] = signal_block[i];
	}
	// Zero padding
	for (int i = filter_len_; i < fft_len_; ++i) {
		time_domain_buffer[i] = 0.0f;
	}

	// Perform forward FFT transform
	kiss_fftr(forward_fft_, &time_domain_buffer[0], &freq_domain_buffer[0]);

	// Complex multiplication in frequency domain with transformed kernel.
	for (int i = 0; i < filtered_freq_domain_buffer_.size(); ++i) {
		filtered_freq_domain_buffer_[i].r =
				freq_domain_buffer[i].r * kernel_freq_domain_buffer_[i].r -
				freq_domain_buffer[i].i * kernel_freq_domain_buffer_[i].i;
		filtered_freq_domain_buffer_[i].i =
				freq_domain_buffer[i].r * kernel_freq_domain_buffer_[i].i +
				freq_domain_buffer[i].i * kernel_freq_domain_buffer_[i].r;
	}

	// Perform inverse FFT transform of filtered_freq_domain_buffer_ and store result back in signal_time_domain_buffer_
	kiss_fftri(inverse_fft_, &filtered_freq_domain_buffer_[0], &time_domain_buffer[0]);

	// Invert FFT scaling
	for (int i = 0; i < fft_len_; ++i) {
		time_domain_buffer[i] /= fft_len_;
	}
}

void FFTFilter::GetResult(vector<float>* signal_block) {
	assert(signal_block);
	signal_block->resize(filter_len_);

	int curr_buf = buffer_selector_;
	int prev_buf = !buffer_selector_;
	for (int i = 0; i < filter_len_; ++i) {
		(*signal_block)[i] = signal_time_domain_buffer_[curr_buf][i]
				+ signal_time_domain_buffer_[prev_buf][i + filter_len_]; // Add overlap from previous FFT transform.
	}
}

