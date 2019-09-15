#ifndef __T_MODEL_PRINT_H__
#define __T_MODEL_PRINT_H__

#include <cstdio>
#include "TModel.h"

template<int nvars, int ord, class T>
void TModel_fprint(FILE *f, const TModel<nvars,ord,T> &m)
{
	typedef TModel<nvars,ord,T>		t_TM;
	typedef typename t_TM::polynomial_type	t_polynomial;

        int	exp[nvars];
	t_polynomial::exponents(0, exp);

	for (int i = 0;i < t_polynomial::size;++i) {
		//fprintf(f, "a%d = %0.15e", i, m.polynomial()[i]);
		fprintf(f, "a%d = %18.15g", i, m.polynomial()[i]);
		for (int j = 0;j < nvars;++j) {
			fprintf(f, " %d", exp[j]);
		}
		fprintf(f, "\n");
		t_polynomial::next_exponents(nvars, exp);
	}
        //fprintf(f, "remainder = [%0.15e,%0.15e]\n", m.remainder().lower(), m.remainder().upper());
        fprintf(f, "remainder = [%0.15g,%0.15g]\n", m.remainder().lower(), m.remainder().upper());
}

template<int nvars, int ord, class T>
void TModel_print(const TModel<nvars,ord,T> &m)
{
	TModel_fprint<nvars,ord,T>(stdout, m);
}

template<int nvars, int ord, class T>
void TModel_fprint_sc(FILE *f, const TModel<nvars,ord,T> &m, T sc[nvars])
{
	typedef TModel<nvars,ord,T>		t_TM;
	typedef typename t_TM::polynomial_type	t_polynomial;

        int	exp[nvars];
	t_polynomial::exponents(0, exp);

	for (int i = 0;i < t_polynomial::size;++i) {
		T	scale(1.f);
                for (int j = 0;j < nvars;++j) scale *= pow(sc[j],exp[j]);
		//fprintf(f, "a%d = %0.15e", i, m.polynomial()[i]);
		fprintf(f, "a%d = %18.15g", i, m.polynomial()[i]/scale);
		for (int j = 0;j < nvars;++j) {
			fprintf(f, " %d", exp[j]);
		}
		fprintf(f, "\n");
		t_polynomial::next_exponents(nvars, exp);
	}
        //fprintf(f, "remainder = [%0.15e,%0.15e]\n", m.remainder().lower(), m.remainder().upper());
        fprintf(f, "remainder = [%0.15g,%0.15g]\n", m.remainder().lower(), m.remainder().upper());
}

template<int nvars, int ord, class T>
void TModel_print_sc(const TModel<nvars,ord,T> &m, T sc[nvars])
{
	TModel_fprint_sc<nvars,ord,T>(stdout, m, sc);
}

#endif
