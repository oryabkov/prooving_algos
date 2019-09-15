
#include <cstdio>
#include <boost/numeric/interval.hpp>
#include "t_vec_tml.h"
#include "TModel.h"
#include "TModelPrint.h"
#include "TModelPicardIteration.h"

template<class T>
struct t_quadratic_ode
{
	typedef t_vec_tml<2,T>		vec_type;
	static const int		dim = 2;

	void operator()(const vec_type& x, const T &t, vec_type& f)const
	{
		f.x[0] = x.x[1];
		f.x[1] = x.x[0]*x.x[0];
	}
};

typedef TModel<3,3,double>						t_TM;
typedef t_vec_tml<2,t_TM>						t_TM_vec;
typedef t_quadratic_ode<t_TM>						t_TM_ode;
typedef boost::numeric::interval<double>				t_interval;
typedef TModelPicardIteration<t_TM_ode,t_TM_vec,t_TM,t_interval,double>	t_picard_iteration;

int main()
{
	t_TM_vec		x;
        t_TM_ode		ode;
        t_picard_iteration      picard_iteration;

	x[0].polynomial()[0] = 1.;
	int	exp1[3] = {0,1,0};
	x[0].polynomial()[t_TM::polynomial_type::term_index(exp1)] = 0.05;
	x[1].polynomial()[0] = -1.;
	int	exp2[3] = {0,0,1};
	x[1].polynomial()[t_TM::polynomial_type::term_index(exp2)] = 0.05;

        t_interval		t(0.);
        double			dt = 0.1;

        picard_iteration.doStep(ode, x, t, dt);

	return 0;
}