#include <vector>

struct kiss_fft_state;

class FFTFilter {
public:
	FFTFilter(int filter_len);

private:
	int filter_len_;
	std::vector<float> filter_state_;
	std::vector<float> window_;

	std::vector<float> fft_time_domain_buffer_;

	kiss_fft_state* forward_fft_;
	kiss_fft_state* inverse_fft_;
};

