#include <iostream>
#include "CRM_utils.h"
#include "distr.h"
#include "CRM.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace distr;

namespace py = pybind11;

class Py_cont_distr : public cont_distr {
public:
	using cont_distr::cont_distr;

	virtual double pdf(double x) const override {
		PYBIND11_OVERRIDE_PURE(
			double, cont_distr, pdf, x
		);
	}
	virtual double cdf(double x) const override {
		PYBIND11_OVERRIDE_PURE(
			double, cont_distr, cdf, x
		);
	}
	virtual double mean() const override {
		PYBIND11_OVERRIDE_PURE(
			double, cont_distr, mean
		);
	}
};

class Py_discrete_distr : public discrete_distr {
public:
	using discrete_distr::discrete_distr;

	virtual double pmf(int n) const override {
		PYBIND11_OVERRIDE_PURE(
			double, discrete_distr, pmf, n
		);
	}
	virtual CRM_utils::complex gf(CRM_utils::complex x) const override {
		PYBIND11_OVERRIDE_PURE(
			CRM_utils::complex, discrete_distr, gf, x
		);
	}
};

PYBIND11_MODULE(CRM, m) {
	py::class_<CRM_c>(m, "CRM")
		.def(py::init<cont_distr*, discrete_distr*>())
		.def("cdf", &CRM_c::cdf);

	py::class_<discrete_distr, Py_discrete_distr>(m, "discrete_distr")
		.def(py::init<double, double>())
		.def("pmf", &discrete_distr::pmf)
		.def("gf", &discrete_distr::gf);

	py::class_<nb, discrete_distr>(m, "nb")
		.def(py::init<double, double>());

	py::class_<cont_distr, Py_cont_distr>(m, "cont_distr")
		.def(py::init<double>())
		.def("pdf", &cont_distr::pdf)
		.def("cdf", &cont_distr::cdf)
		.def("mean", &cont_distr::mean);

	py::class_<ruin_special, cont_distr>(m, "ruin_special")
		.def(py::init<cont_distr*>());

	py::class_<expon, cont_distr>(m, "expon")
		.def(py::init<double>())
		.def(py::init<double, double>());

	py::class_<gamma, cont_distr>(m, "gamma")
		.def(py::init<double, double>())
		.def(py::init<double, double, double>());

	py::class_<lognormal, cont_distr>(m, "lognormal")
		.def(py::init<double, double>())
		.def(py::init<double, double, double>());

	py::class_<pareto, cont_distr>(m, "pareto")
		.def(py::init<double, double>())
		.def(py::init<double, double, double>());
}
