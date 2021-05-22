#include "CRM.h"
#include <cmath>
#include "CRM_utils.h"

const long long _n = 10000;
const long long c = 10, s = 2;
const long long ln = c * pow(static_cast<long long>(log(_n)), s);

std::vector<double> CRM_c::discretize_pdf(long long n) const {
	std::vector<double> g(n+1);
	double t = 0;
	double lim = eta->limit;

	long long ln_n = c * pow(static_cast<long long>(log(n)), s);
	double ratio;

	if (lim < INFINITY) {
		ratio = lim / static_cast<double>(n);

		for (int i = 0; i < n; i++) {
			g[i] = eta->pdf((i + 1 / 2) * ratio) * ratio;
			t += g[i];
		}
	}
	else {
		ratio = static_cast<double>(ln_n) / static_cast<double>(n);

		for (int i = 0; i < n; i++) {
			g[i] = eta->pdf((i + 1 / 2) * ratio) * ratio;
			t += g[i];
		}
	}
	g[n] = 1 - t;
	
	return g;
}

std::vector<double> CRM_c::cdf(double x) const {
	long long _x;

	if (eta->limit != INFINITY)
		_x = static_cast<long long>(x * static_cast<double>(_n) / eta->limit);
	else
		_x = static_cast<long long>(x * static_cast<double>(_n) / static_cast<double>(ln));

	double ans_r = 0, ans_i = 0;

	std::vector<double> f_r = recursive(_x);

	for (auto x : f_r)
		ans_r += x;

	long long r = 1;

	while ((r < _n) or (r < _x))
		r *= 2;

	std::vector<double> f_i = inverse(8 * r);

	for (int i = 0; i <= _x; i++) 
		ans_i += f_i[i];

	return std::vector<double>{ans_r, ans_i};
}

std::vector<double> CRM_c::inverse(long long x) const {
	std::vector<double> t = discretize_pdf(_n);

	if (x > t.size()) {
		std::vector<double> t1(x - t.size(), 0);
		t.reserve(x);
		t.insert(t.end(), t1.begin(), t1.end());
	}

	CRM_utils::array g{ t };
	
	g = g.fft();

	for (int i = 0; i < g.size; i++)
		g[i] = N->gf(g[i]);

	return g.fft(true).re();
}

std::vector<double> CRM_c::recursive(long long x) const {
	std::vector<double> g = discretize_pdf(_n);

	double p0 = N->pmf(0), p1 = N->pmf(1), a = N->a, b = N->b;

	if (x + 1 > g.size()) {
		std::vector<double> t1(x + 1 - g.size(), 0);
		g.reserve(x+1);
		g.insert(g.end(), t1.begin(), t1.end());
	}

	std::vector<double> f(x+1);

	//add support for posion and log
	f[0] = (p1 + (a + b) * (1 - p0)) / (a + b) * pow((1 - a) / (1 - a * g[0]), 1 + b / a) - (p1 - (a + b) * p0) / (a + b);

	for (int k = 1; k <= x; k++) {
		double t = 0;
		for (int j = 1; j <= k; j++) {
			double q = static_cast<double>(j) / static_cast<double>(k);
			t += (a + q * b) * g[j] * f[k - j];
		}

		f[k] = (g[k] * (p1 - p0 * (a + b)) + t) / (1 - a * g[0]);
	}

	return f;
}

