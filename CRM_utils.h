#ifndef _CRM_utils_
#define _CRM_utils_

#include <initializer_list>
#include <stdexcept>
#include <algorithm>
#include <vector>

namespace CRM_utils {
	struct array;

	struct complex {
		double re;
		double im;
		complex(double a = 0., double b = 0.) : re{ a }, im{ b } {}
		complex(const complex& z) : re{ z.re }, im{ z.im } {}
		complex operator* (const complex& b) const noexcept;
		complex operator+ (const complex& b) const noexcept;
		complex operator- (const complex& b) const noexcept;
		complex operator/ (const complex& b) const noexcept;
		complex& operator= (const complex& b) { re = b.re; im = b.im; return *this; }
		array roots(int n);
		complex root(int n, int k);
		inline double& operator[](int i) {
			if (i == 0)
				return re;
			else if (i == 1)
				return im;
			else
				throw std::out_of_range{ "a" };
		}
	};

	struct array {
		complex* arr;
		unsigned size;
		array(const std::initializer_list<complex>& l);
		array(unsigned n) : size{ n }, arr{ new complex[n] } {}
		array(array&& a) noexcept { size = a.size; arr = a.arr; a.arr = nullptr; }
		array(const array& a);
		array(const std::vector<double>& a);
		~array() { delete[] arr; };
		inline complex& operator[](const unsigned long long n) const { if ((0 <= n) && (n < size)) return *(arr + n); else throw std::out_of_range{ "a" }; }
		complex* begin() const { return arr; }
		complex* end() const { return arr + size; }
		array operator*(const array& x);
		array& operator=(const array& rhs);
		array& operator=(array&& rhs) noexcept;
		array fft(bool inverse = false);
		array bit_reversal();
		std::vector<double> re() const { std::vector<double> ans(size); for (unsigned long long i = 0; i < size; i++) ans[i] = arr[i].re; return ans; }
	};

	namespace complex_const {
		const double PI = 3.14159265358979323846;
		const complex j{ 0, 1 };
	}
}

int nearest_2_pow(int x) noexcept;

namespace cc = CRM_utils::complex_const;

#endif _CRM_utils_
