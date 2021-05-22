#include "CRM_utils.h"
using namespace CRM_utils;

// some methods implementations
array& array::operator=(const array& rhs) {
	delete[] arr;
	size = rhs.size;
	arr = new complex[size];

	for (int i = 0; i < size; i++)
		arr[i] = rhs.arr[i];

	return *this;
}

array& array::operator=(array&& rhs) noexcept {
	arr = rhs.arr;
	rhs.arr = nullptr;

	return *this;
}

array::array(const array& a) {
	size = a.size;
	arr = new complex[size];

	for (int i = 0; i < size; i++)
		arr[i] = a.arr[i];

}

array::array(const std::initializer_list<complex>& l) {
	size = l.size();
	arr = new complex[size];

	for (unsigned i = 0; i < size; i++)
		arr[i] = *(l.begin() + i);
}

array::array(const std::vector<double>& l) {
	size = l.size();
	arr = new complex[size];

	for (unsigned i = 0; i < size; i++)
		arr[i] = l[i];
}

array array::operator*(const array& x) {
	array res(size);

	for (unsigned i = 0; i < size; i++)
		res[i] = arr[i] * x.arr[i];

	return res;
}

// algortihm implementations

int nearest_2_pow(int x) noexcept {
	int ans = 1;

	for (; x > 0; x >>= 1) ans *= 2;

	return ans;
}

array array::bit_reversal() {
	array y(size);

	auto rev = [](int i, int n) { // this function returns the reversed bit sequence of i, n is the size of bit sequence
		int len = 0; int res = 0; int power = 1; for (; n - 1 > 0; n >>= 1) len++;

		for (int k = 1; k <= len; k++) {
			res += ((i >> (len - k)) % 2) * power;
			power *= 2;
		}
		return res;
	};

	for (int i = 0; i < size; i++)
		y[rev(i, size)] = *(arr + i);

	return y;
}

array array::fft(bool inverse) {
	array y(size);

	y = bit_reversal();

	int n = size;
	int m = 1;

	for (int s = 1; s <= log2(n); s++) {
		m *= 2;
		complex wn = complex{ 1 }.root(m, inverse ? m - 1 : 1);

		for (int k = 0; k < n; k += m) {
			complex w = 1;

			for (int j = 0; j < m / 2; j++) {
				complex t = w * y[k + j + m / 2];
				complex u = y[k + j];

				//if (inverse)
				//	std::cout << y[k + j].re << " + i*" << y[k + j].im << ' ' << t.re << " + i*" << t.im << ' ' << ' ' << u.re << " + i*" << u.im << ' ' << j << ' ' << m << ' ' << k << '\n';

				y[k + j] = u + t;
				y[k + j + m / 2] = u - t;
				w = w * wn;
				//if(inverse)
				//std::cout << y[k+j].re << " + i*" << y[k+j].im << ' ' << j << ' ' << m << ' ' << k << '\n';
			}
		}
	}

	if (inverse)
		for (int i = 0; i < size; i++) {
			//std::cout << y[i].re << ' ';

			y[i] = y[i] / size;
		}

	return y;
}