#include "CRM_utils.h"
#include <algorithm>

using namespace CRM_utils;
using namespace cc;

complex complex::operator*(const complex& b) const noexcept {
	double r = re * b.re - im * b.im;
	double i = re * b.im + im * b.re;

	return complex{ r, i };
}

complex complex::operator/(const complex& b) const noexcept {
	double r = (re*b.re - im*b.im)/(b.re*b.re+b.im*b.im);
	double i = (im * b.re - re * b.im) / (b.re * b.re + b.im * b.im);

	return complex{ r, i };
}

complex complex::operator+(const complex& b) const noexcept {
	double r = re + b.re;
	double i = im + b.im;

	return complex{ r, i };
}

complex complex::operator-(const complex& b) const noexcept {
	double r = re - b.re;
	double i = im - b.im;

	return complex{ r, i };
}

array complex::roots(int n) {
	double r = pow(re * re + im * im, 1 / (2 * n));
	double phi;

	if ((re > 0) && (im >= 0))
		phi = atan(im / re);
	else if ((re > 0) && (im < 0))
		phi = 2 * PI + atan(im / re);
	else if (re == 0)
		phi = PI / 2 + (im > 0 ? 0 : PI);
	else
		phi = PI + atan(im / re);

	double sphi = phi / n; // arguments of roots

	array rs(n);

	std::for_each(rs.begin(), rs.end(), [r, phi, &sphi, n](complex& w) {
		w.re = r * cos(sphi); w.im = r * sin(sphi); sphi += 2 * PI / n;
	});

	return rs;
}

complex complex::root(int n, int k) {
	double r = pow(re * re + im * im, 1 / static_cast<double>(2 * n));
	double phi;

	if ((re >= 0) && (im >= 0))
		phi = atan(im / re);
	else if ((re >= 0) && (im >= 0))
		phi = 2 * PI + atan(im / re);
	else
		phi = PI + atan(im / re);

	double sphi = (phi + 2 * PI * k) / n;

	return complex{ r * cos(sphi), r * sin(sphi) };
}
