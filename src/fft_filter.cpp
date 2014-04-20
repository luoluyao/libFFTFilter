#include <assert.h>
#include <cmath>

#include "fft_filter.h"

using namespace std;

FFTFilter::FFTFilter(int filter_len) :
		filter_len_(filter_len),
		filter_state_(filter_len, 0.0f),
		kernel_defined_(false),
		kernel_time_domain_buffer_(2 * filter_len),
		kernel_freq_domain_buffer_(2 * filter_len),
		buffer_selector_(0),
		signal_time_domain_buffer_(2, vector<kiss_fft_cpx>(2 * filter_len)),
		signal_freq_domain_buffer_(2, vector<kiss_fft_cpx>(2 * filter_len)),
		filtered_freq_domain_buffer_(2 * filter_len) {

	bool is_power_of_two = ((filter_len_ != 0)
			&& !(filter_len_ & (filter_len_ - 1)));
	assert(is_power_of_two && "Filter length must be a power of 2");

	forward_fft_ = kiss_fft_alloc(filter_len_ * 2, 0, 0, 0);
	inverse_fft_ = kiss_fft_alloc(filter_len_ * 2, 0, 0, 0);

	Init();
}

FFTFilter::~FFTFilter() {
	kiss_fft_free(forward_fft_);
	kiss_fft_free(inverse_fft_);
}

void FFTFilter::Init() {
	// Initialize all buffers with zeros.
	memset(&kernel_time_domain_buffer_[0], 0,
			sizeof(kiss_fft_cpx) * 2 * filter_len_);
	memset(&kernel_freq_domain_buffer_[0], 0,
			sizeof(kiss_fft_cpx) * 2 * filter_len_);
	for (int i = 0; i < 2; ++i) {
		memset(&signal_time_domain_buffer_[i][0], 0,
				sizeof(kiss_fft_cpx) * 2 * filter_len_);
		memset(&signal_freq_domain_buffer_[i][0], 0,
				sizeof(kiss_fft_cpx) * 2 * filter_len_);
	}
}

void FFTFilter::SetKernel(const std::vector<float>& kernel) {
	assert(kernel.size()==filter_len_ && "Kernel size must match filter length");

	for (int i = 0; i < filter_len_; ++i) {
		kernel_time_domain_buffer_[i].r = kernel[i];
		kernel_time_domain_buffer_[i].i = 0.0f;
	}
	// Zero padding
	for (int i = filter_len_; i < filter_len_ * 2; ++i) {
		kernel_time_domain_buffer_[i].r = 0.0f;
		kernel_time_domain_buffer_[i].i = 0.0f;
	}

	// Perform forward FFT transform
	kiss_fft(forward_fft_, &kernel_time_domain_buffer_[0],
			&kernel_freq_domain_buffer_[0]);

	kernel_defined_ = true;
}

void FFTFilter::AddSignalBlock(const vector<float>& signal_block) {
	assert(signal_block.size()==filter_len_ && "Signal block size must match filter length");
	assert(!kernel_defined_ && "No suitable kernel defined");

	// Switch buffer selector
	buffer_selector_ = !buffer_selector_;

	vector<kiss_fft_cpx>& time_domain_buffer = signal_time_domain_buffer_[buffer_selector_];
	vector<kiss_fft_cpx>& freq_domain_buffer = signal_time_domain_buffer_[buffer_selector_];


	for (int i = 0; i < filter_len_; ++i) {
		time_domain_buffer[i].r = signal_block[i];
		time_domain_buffer[i].i = 0.0f;
	}
	// Zero padding
	for (int i = filter_len_; i < filter_len_ * 2; ++i) {
		time_domain_buffer[i].r = 0.0f;
		time_domain_buffer[i].i = 0.0f;
	}

	// Perform forward FFT transform
	kiss_fft(forward_fft_, &time_domain_buffer[0], &freq_domain_buffer[0]);

	// Complex multiplication in frequency domain with transformed kernel.
	for (int i = 0; i < filter_len_ * 2; ++i) {
		filtered_freq_domain_buffer_[i].r =
				freq_domain_buffer[i].r * kernel_freq_domain_buffer_[i].r -
				freq_domain_buffer[i].i * kernel_freq_domain_buffer_[i].i;
		filtered_freq_domain_buffer_[i].i =
				freq_domain_buffer[i].r * kernel_freq_domain_buffer_[i].i +
				freq_domain_buffer[i].i * kernel_freq_domain_buffer_[i].r;
	}

	// Perform inverse FFT transform of filtered_freq_domain_buffer_ and store result back in signal_time_domain_buffer_
	kiss_fft(inverse_fft_, &filtered_freq_domain_buffer_[0], &time_domain_buffer[0]);
}

void FFTFilter::GetResult(vector<float>* signal_block) {
	assert(signal_block);
	signal_block->resize(filter_len_);

	int curr_buf = buffer_selector_;
	int prev_buf = !buffer_selector_;
	for (int i = 0; i < filter_len_; ++i) {
		(*signal_block)[i] = signal_time_domain_buffer_[curr_buf][i].r
				+ signal_time_domain_buffer_[prev_buf][i + filter_len_].r; // Add overlap from previous FFT transform.
	}
}

