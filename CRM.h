#ifndef _CRM_
#define _CRM_
#include "CRM_utils.h"
#include "distr.h"
#include <vector>


class CRM_c {
	distr::cont_distr* eta;
	distr::discrete_distr* N;
	std::vector<double> discretize_pdf(long long n) const;
	std::vector<double> recursive(long long x) const;
	std::vector<double> inverse(long long x) const;
public:
	CRM_c(distr::cont_distr* _eta, distr::discrete_distr* _N) : eta{ _eta }, N{ _N } { }
	std::vector<double> cdf(double x) const;
};


#endif _CRM_