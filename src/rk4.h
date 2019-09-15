#ifndef __RKCD_H__
#define __RKCD_H__

#include <cmath>
#include <cstdio>

#define RKCD_MAX_DEGREE	30

//VECTOR concept:
//
//ODE_T concept:
//void verify_size(VECTOR&)const
//void operator()(const VECTOR& x,BASE_T t,VECTOR& f)const
//addMul(VECTOR& B, BASE_T mulB, const VECTOR& A, BASE_T mulA) calc: B = B*mulB + A*mulA
//addMul(VECTOR& C, BASE_T mulC, const VECTOR& A, BASE_T mulA, const VECTOR& B, BASE_T mulB) calc: C = C*mulC + A*mulA + B*mulB
template<class ODE_T,class VECTOR,class BASE_T>
class RK4Stepper
{
	mutable VECTOR k1, k2, k3, k4, tmp;
public:
	void	doStep(const ODE_T &ode, VECTOR &y, BASE_T &t, BASE_T dt)const;
};

template<class BASE_T>
class	__RKCD_test_vec_t
{
public:
	BASE_T x;
};

template<class BASE_T>
class	__RKCD_test_ode_t
{
public:
	void verify_size(__RKCD_test_vec_t<BASE_T> &x)const { x.x = 0.; }
	void operator()(const __RKCD_test_vec_t<BASE_T>& x,BASE_T t,__RKCD_test_vec_t<BASE_T>& f)const { f.x = -x.x; }
	void addMul(__RKCD_test_vec_t<BASE_T>& B, BASE_T mulB, const __RKCD_test_vec_t<BASE_T>& A, BASE_T mulA)const { B.x = B.x*mulB + A.x*mulA; }
	void addMul(__RKCD_test_vec_t<BASE_T>& C, BASE_T mulC, const __RKCD_test_vec_t<BASE_T>& A, BASE_T mulA, const __RKCD_test_vec_t<BASE_T>& B, BASE_T mulB)const { C.x = C.x*mulC + A.x*mulA + B.x*mulB; }
};

template<class ODE_T,class VECTOR,class BASE_T>
void RK4Stepper<ODE_T,VECTOR,BASE_T>::doStep(const ODE_T &ode, VECTOR &y, BASE_T &t, BASE_T dt)const
{
	ode.verify_size(k1); ode.verify_size(k2); ode.verify_size(k3); ode.verify_size(k4);
	ode.verify_size(tmp);

    	ode.addMul(tmp, 0.f, y, 1.f);
    	ode(tmp, t        , k1);
    	ode.addMul(tmp, 0.f, y, 1.f, k1, 0.5f*dt);
    	ode(tmp, t+0.5f*dt, k2);
    	ode.addMul(tmp, 0.f, y, 1.f, k2, 0.5f*dt);
    	ode(tmp, t+0.5f*dt, k3);
    	ode.addMul(tmp, 0.f, y, 1.f, k3, dt);
    	ode(tmp, t+dt     , k4);

        ode.addMul(y, 1.f, k1, 1.*dt/6, k2, 2.*dt/6);
        ode.addMul(y, 1.f, k3, 2.*dt/6, k4, 1.*dt/6);

	t+=dt;
}

#endif
