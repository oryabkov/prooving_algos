#ifndef __T_MODEL_H__
#define __T_MODEL_H__

#include <boost/numeric/interval.hpp>
#include "polymul.h"

//supposed to be defined on [-1;1]^nvars cube
template<int nvars, int ord, class T = double>
class TModel
{
public:
	typedef T                               base_type;
	typedef boost::numeric::interval<T>	interval_type;
	typedef polynomial<T,nvars,ord>		polynomial_type;
	static const int			vars_num = nvars;
	static const int			taylor_order = ord;
private:
	typedef typename boost::numeric::interval_lib::unprotect<interval_type>::type	interval_type_f;
	typedef typename interval_type::traits_type::rounding				interval_rounding;

	typedef polynomial<interval_type,nvars,ord>					t_interv_polynomial;
	typedef polynomial<interval_type_f,nvars,ord>					t_interv_polynomial_f;

	static int	*_exp_table;
	static int	&exp_table(int i_term, int i_var)
	{
		return _exp_table[i_term*nvars + i_var];
	}
	
	static bool	ignore_remainder;

	polynomial_type		p;
	interval_type		R;

	//j could be greater than out model order
	interval_type	mon_bnds(int j)const
	{
		return (j == 0?interval_type(1.,1.):interval_type(-1.,1.));
	}
	interval_type_f	mon_bnds_f(int j)const
	{
		return (j == 0?interval_type_f(1.,1.):interval_type_f(-1.,1.));
	}
public:
	static void	init_exponents_table()
	{
		_exp_table = new int[polynomial_type::size*nvars];
		int	exp[nvars];
		polynomial_type::exponents(0, exp);
		for (int j = 0;j < polynomial_type::size;++j) {
			for (int k = 0;k < nvars;++k) exp_table(j, k) = exp[k];
			//for (int k = 0;k < nvars;++k) exp_table(j, k) = 0;
			polynomial_type::next_exponents(nvars, exp);
		}
	}
	static void	free_exponents_table()
	{
		delete []_exp_table;
	}
	static void	set_ignore_remainder_flag(bool _ignore_remainder)
	{
		ignore_remainder = _ignore_remainder;
	}

	explicit TModel(const polynomial_type &_p = polynomial_type(0.f), const interval_type &_R = interval_type(0.f)) : p(_p), R(_R) {}

	polynomial_type		&polynomial() { return p; }
	interval_type		&remainder() { return R; }
	const polynomial_type	&polynomial()const { return p; }
	const interval_type	&remainder()const { return R; }

	interval_type	bnds()const
	{
		interval_type	res(0.);
		res += mon_bnds(0)*interval_type(p[0]);
		do {
			interval_rounding	rnd;
			interval_type_f		tmp(0.);
			for (int j = 1;j < polynomial_type::size;++j) {
				tmp += interval_type_f(p[j])*mon_bnds_f(j);
			}
			res += interval_type(tmp);
		} while(0);
		return res;
	}

	TModel operator+(const interval_type &x)const
	{
		TModel			res;
		res.p = p;
		res.p[0] += median(x);
		//experimental feature
		if (ignore_remainder) return res;
		res.R = R;
		res.R += ((interval_type(p[0]) + x) - interval_type(res.p[0]))*mon_bnds(0);
		return res;
	}
	TModel operator+(const TModel &m)const
	{
		TModel			res;
		for (int j = 0;j < polynomial_type::size;++j) {
			res.p[j] = p[j] + m.p[j];
		}
		//experimental feature
		if (ignore_remainder) return res;
		res.R = R + m.R;
		res.R += ((interval_type(p[0]) + interval_type(m.p[0])) - interval_type(res.p[0]))*mon_bnds(0);
		do {
			interval_rounding	rnd;
			interval_type_f		tmp(0.);
			for (int j = 1;j < polynomial_type::size;++j) {
				tmp += ((interval_type_f(p[j]) + interval_type_f(m.p[j])) - interval_type_f(res.p[j]))*mon_bnds_f(j);
			}
			res.R += interval_type(tmp);
		} while (0);
		return res;
	}
	// TODO ignore_remainder and optimization (interval_rounding)
	TModel operator-(const TModel &m)const
	{
		TModel			res;
		res.R = R - m.R;
		for (int j = 0;j < polynomial_type::size;++j) {
			res.p[j] = p[j] - m.p[j];
			res.R += ((interval_type(p[j]) - interval_type(m.p[j])) - interval_type(res.p[j]))*mon_bnds(j);
		}
		return res;
	}
	// TODO ignore_remainder and optimization (interval_rounding)
	TModel operator*(const T &x)const
	{
		TModel			res;
		res.R = R*interval_type(x);
		for (int j = 0;j < polynomial_type::size;++j) {
			res.p[j] = p[j]*x;
			res.R += ((interval_type(p[j])*interval_type(x)) - interval_type(res.p[j]))*mon_bnds(j);
		}
		return res;
	}
	TModel operator*(const interval_type &x)const
	{
		TModel		res;
		T		med = median(x);
		for (int j = 0;j < polynomial_type::size;++j) {
			res.p[j] = p[j]*med;
		}
		//experimental feature
		if (ignore_remainder) return res;
		res.R = R*x;
		res.R += ((interval_type(p[0])*interval_type(x)) - interval_type(res.p[0]))*mon_bnds(0);
		do {
			interval_rounding	rnd;

			interval_type_f	tmp(0.), x_f(x);
			for (int j = 1;j < polynomial_type::size;++j) {
				tmp += ((interval_type_f(p[j])*x_f) - interval_type_f(res.p[j]))*mon_bnds_f(j);
			}
			res.R += interval_type(tmp);
		} while(0);
		return res;
	}
	TModel operator*(const TModel &m)const
	{
		TModel		res;
		taylormul(res.p, p, m.p);
		//experimental feature
		if (ignore_remainder) return res;

		res.R = R*m.R;
		res.R += bnds()*m.R + m.bnds()*R;

		t_interv_polynomial	ipoly1,ipoly2;
		for (int j = 0;j < polynomial_type::size;++j) {
			ipoly1[j] = p[j];
			ipoly2[j] = m.p[j];
		}

		//low order error estimation
		t_interv_polynomial	polymul_exact;
		do {
			interval_rounding	rnd;
			t_interv_polynomial_f	ipoly1_f, ipoly2_f,
						polymul_exact_f;
			for (int j = 0;j < polynomial_type::size;++j) {
				ipoly1_f[j] = ipoly1[j];
				ipoly2_f[j] = ipoly2[j];
			}
			taylormul(polymul_exact_f, ipoly1_f, ipoly2_f);
			for (int j = 0;j < polynomial_type::size;++j) {
				polymul_exact[j] = polymul_exact_f[j];
			}
		} while(0);

		//TODO check whether j correspondence in polymul holds true
		res.R += (polymul_exact[0] - interval_type(res.p[0]))*mon_bnds(0);
		do {
			interval_rounding	rnd;
			interval_type_f		tmp(0.);
			for (int j = 1;j < t_interv_polynomial::size;++j) {
				tmp += (interval_type_f(polymul_exact[j]) - interval_type_f(res.p[j]))*mon_bnds_f(j);
			}
			res.R += interval_type(tmp);
		} while(0);

		int		polylens[ord+1];
		for (int curr_ord = 0;curr_ord <= ord;++curr_ord) {
			polylens[curr_ord] = polylen(nvars, curr_ord);
		}
		int		curr_ord;

		//high order error estimation
		interval_type	poly2_bnds_by_ord[ord+1];
		poly2_bnds_by_ord[0] = mon_bnds(0)*ipoly2[0];
		do {
			interval_rounding	rnd;
			interval_type_f		poly2_bnds_by_ord_f[ord+1];
			for (int curr_ord = 1;curr_ord <= ord;++curr_ord) poly2_bnds_by_ord_f[curr_ord] = interval_type_f(0.);
			curr_ord = 1;
			for (int j = 1;j < polynomial_type::size;++j) {
				if (j == polylens[curr_ord]) ++curr_ord;

				poly2_bnds_by_ord_f[curr_ord] += interval_type_f(ipoly2[j])*mon_bnds_f(j);
			}
			for (int curr_ord = 1;curr_ord <= ord;++curr_ord) poly2_bnds_by_ord[curr_ord] = poly2_bnds_by_ord_f[curr_ord];
		} while (0);
		curr_ord = 0;
		for (int poly2_ord = ord - curr_ord + 1;poly2_ord <= ord;++poly2_ord) {
			res.R += mon_bnds(0)*ipoly1[0]*poly2_bnds_by_ord[poly2_ord];
		}
		do {
			interval_rounding	rnd;
			interval_type_f		tmp(0.);
			curr_ord = 1;
			for (int j = 1;j < polynomial_type::size;++j) {
				if (j == polylens[curr_ord]) ++curr_ord;

				for (int poly2_ord = ord - curr_ord + 1;poly2_ord <= ord;++poly2_ord) {
					tmp += interval_type_f(ipoly1[j])*mon_bnds_f(j)*interval_type_f(poly2_bnds_by_ord[poly2_ord]);
				}
			}
			res.R += interval_type(tmp);
		} while (0);

		return res;
	}

	TModel	part_eval(int ivar, const T &x)
	{
		TModel			res;
		t_interv_polynomial	res_exact;
		res.R = R;
		for (int j = 0;j < polynomial_type::size;++j) {
			res.p[j] = T(0.f);
			res_exact[j] = interval_type(0.f);
		}
		int	exp[nvars];
		for (int j = 0;j < polynomial_type::size;++j) {
			for (int kk = 0;kk < nvars;++kk) exp[kk] = exp_table(j, kk);

			T		x_pow = (exp[ivar] == 0 ? T(1.f) : pow(x,exp[ivar]));
			interval_type	x_pow_exact = (exp[ivar] == 0 ? interval_type(1.f) : pow(interval_type(x),exp[ivar]));
			exp[ivar] = 0;
			int	j_res = polynomial_type::term_index(exp);
			interval_type	res_coeff;
			res.p[j_res] += p[j]*x_pow;
			res_exact[j_res] += interval_type(p[j])*x_pow_exact;
		}
		for (int j = 0;j < polynomial_type::size;++j) {
			res.R += (res_exact[j] - interval_type(res.p[j]))*mon_bnds(j);
		}
		return res;
	}
	TModel	part_integr(int ivar, const T &lower_pnt)
	{
		//TODO how to check whether we are right here??
		TModel	res;
		//first, work with polynomial part only
		res.R = interval_type(0.f);
		for (int j = 0;j < polynomial_type::size;++j) res.p[j] = T(0.f);
		int	exp[nvars];
		for (int j = 0;j < polynomial_type::size;++j) {
			for (int kk = 0;kk < nvars;++kk) {
				exp[kk] = exp_table(j, kk);
			}

			exp[ivar]++;
			int		j_new = polynomial_type::term_index(exp);
			interval_type	res_coeff;
			if (j_new < polynomial_type::size) {
				res.p[j_new] = p[j]/T(exp[ivar]);
				res_coeff = interval_type(res.p[j_new]);
			} else
				res_coeff = interval_type(0.f);
			res.R += (interval_type(p[j])/interval_type(T(exp[ivar])) - res_coeff)*mon_bnds(j_new);
		}
		res = res - res.part_eval(ivar, lower_pnt);
		//add remainder integral
		//TODO we can make it more fittable (get R.center() as constant and integrate it)
		res.R += R*(interval_type(-1.f,1.f)-interval_type(lower_pnt));
		return res;
	}
};

template<int nvars, int ord, class T>
int *TModel<nvars, ord, T>::_exp_table;

template<int nvars, int ord, class T>
bool TModel<nvars, ord, T>::ignore_remainder = false;

#endif
