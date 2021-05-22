#ifndef _DISTR_
#define _DISTR_

#include "CRM_utils.h"
#include <cmath>

namespace {
	double binom(double n, int k) {
		double res = 1;

		for (long long i = 0; i < k; i++)
			res *= (n - i) / (i + 1L);

		return res;
	}
}

namespace distr {
	class cont_distr {
	public:
		double limit;
		cont_distr(double _limit) : limit { _limit } {}
		virtual double pdf(double x) const = 0;
		virtual double cdf(double x) const = 0;
		virtual double mean() const = 0;
	};

	class discrete_distr {
	public:
		double a;
		double b;
		discrete_distr(double _a, double _b) : a{ _a }, b{ _b } {}
		virtual double pmf(int n) const = 0;
		virtual CRM_utils::complex gf(CRM_utils::complex x) const = 0;
	};

	class expon : public cont_distr {
		double l;
	public:
		expon(double lam, double lim) : l {lam}, cont_distr(lim) {}
		expon(double lam) : l{ lam }, cont_distr(INFINITY) { }
		double pdf(double x) const { if ((x > 0) and (x <= limit)) return 1/l * exp(-x / l); else return 0; }
		double cdf(double x) const { if ((x > 0) and (x <= limit)) return 1 - exp(-x / l); else if (x > limit) return 1; else return 0; }
		double mean() const { if (limit == INFINITY) return l; else return (l - exp(-limit / l) * (limit + l) + limit * (1 - cdf(limit))); }
	};

	class nb : public discrete_distr {
		double r;
		double p;
	public:
		nb(double _r, double _p) : r{ _r }, p{ _p }, discrete_distr(_p, _p * (_r - 1)) {  }
		double pmf(int n) const { return binom(r + n - 1, n-1) * pow(p, n) * pow(1 - p, r); }
		CRM_utils::complex gf(CRM_utils::complex x) const {
			if (x.re * x.re + x.im * x.im > 1) 
				return CRM_utils::complex{ 0, 0 };
			else {
				CRM_utils::complex z = p - 1;
				
				return z / (x * p - 1); //extend this to whole nb distr, not just geometric
			}
		}
	};

	class gamma : public cont_distr {
		double alpha;
		double beta;
	public:
		gamma(double _a, double _b) : alpha{ _a }, beta{ _b }, cont_distr(INFINITY) { }
		gamma(double _a, double _b, double limit) : alpha{ _a }, beta{ _b }, cont_distr(limit) { }
		double pdf(double x) const { if (x <= limit) return pow(beta, alpha) * pow(x, alpha - 1) * exp(-beta * x) / tgamma(alpha); else return 0; }
		double cdf(double x) const { return 1; }
		double mean() const { return alpha / beta; }
	};

	class pareto : public cont_distr {
		double xm;
		double alpha;
	public:
		pareto(double _x, double _a) : xm{ _x }, alpha{ _a }, cont_distr(INFINITY) {}
		pareto(double _x, double _a, double limit) : xm{ _x }, alpha{ _a }, cont_distr(limit) {}
		double pdf(double x) const { if ((xm <= x) and (x <= limit)) return alpha * pow(xm, alpha) / pow(x, alpha + 1); else return 0; }
		double cdf(double x) const { if ((xm <= x) and (x <= limit)) return 1 - pow(xm / x, alpha); else if (x < xm) return 0; else return 1; }
		double mean() const { return alpha * xm / (alpha - 1); }
	};

	class lognormal : public cont_distr {
		double mu;
		double sigma;
	public:
		lognormal(double _m, double _s) : mu{ _m }, sigma{ _s }, cont_distr(INFINITY) { }
		lognormal(double _m, double _s, double limit) : mu{ _m }, sigma{ _s }, cont_distr(limit) { }
		double pdf(double x) const { if (x <= limit) return 1/(sigma*x*sqrt(2*cc::PI))*exp(-pow(log(x)-mu, 2)/(2*pow(sigma,2))); else return 0; }
		double cdf(double x) const { return 1/2*(1+erf((log(x)-mu)/(sigma*sqrt(2)))); }
		double mean() const { return exp(mu+pow(sigma, 2)/2); }
	};

	class ruin_special : public cont_distr {
		cont_distr* eta;
	public:
		ruin_special(cont_distr* _eta) : eta{ _eta }, cont_distr( _eta->limit ) { }
		double pdf(double x) const { return (1 / eta->mean()) * (1 - eta->cdf(x)); }
		double cdf(double x) const { return -1; }
		double mean() const noexcept { return -1; }
	};
}

#endif _DISTR_