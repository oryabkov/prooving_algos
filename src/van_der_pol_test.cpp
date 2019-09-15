
#include "rk4.h"

struct t_vec2
{
	double x,y;
};

struct t_van_der_pol_ode
{
	double	mu;

	void verify_size(t_vec2& x)const {}
	void operator()(const t_vec2& x,double t,t_vec2& f)const
	{
		f.x = x.y;
		f.y = mu*(1.f-x.x*x.x)*x.y - x.x;
	}
	void addMul(t_vec2& B, double mulB, const t_vec2& A, double mulA)const
	{
		B.x = B.x*mulB + A.x*mulA;
		B.y = B.y*mulB + A.y*mulA;
	}
	void addMul(t_vec2& C, double mulC, const t_vec2& A, double mulA, const t_vec2& B, double mulB)const
	{
		C.x = C.x*mulC + A.x*mulA + B.x*mulB;
		C.y = C.y*mulC + A.y*mulA + B.y*mulB;
	}
};

typedef RK4Stepper<t_van_der_pol_ode,t_vec2,double> t_stepper;

int main(int argc,char **args)
{
        t_vec2			st;
        t_van_der_pol_ode       ode;
	t_stepper		stepper;
	
	st.x = 0.01f; st.y = 0.01f;
	ode.mu = 1.f;

	double	T = 100.f,
		t = 0.f,
		//dt = 1e-2f;
		dt = 1e-3f;
	while (t < T) {
		stepper.doStep(ode, st, t, dt);
		printf("%f %f %f\n", t, st.x, st.y);
	}

	return 0;
}