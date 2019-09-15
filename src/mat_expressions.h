#ifndef __MAT_EXPRESSIONS_H__
#define __MAT_EXPRESSIONS_H__

#include "mat_expr_base.h"

#ifndef __CUDACC__
#define __DEVICE_TAG__
#else
#define __DEVICE_TAG__ __device__ __host__
#endif

//E is VECTOR_EXPRESSION
template<class E>
class mat_expr_negate : public mat_expr_base< mat_expr_negate<E> >
{
public:
	static const int		sz1 = E::sz1;
	static const int		sz2 = E::sz2;
	typedef typename E::value_type	value_type;

        __DEVICE_TAG__ mat_expr_negate(const E	&_e) : e(_e)
	{
		//printf("mat_expr_negate created\n");
        }
        /*__DEVICE_TAG__ ~mat_expr_negate()
        {
		printf("mat_expr_negate deleted\n");
	}*/
	__DEVICE_TAG__ value_type	operator()(int i1,int i2)const
	{
		return -e(i1,i2);
	}
private:
	const E	&e;
};

//E is VECTOR_EXPRESSION
template<class E>
class mat_expr_scalar_mul : public mat_expr_base< mat_expr_scalar_mul<E> >
{
public:
	static const int		sz1 = E::sz1;
	static const int		sz2 = E::sz2;
	typedef typename E::value_type	value_type;

	//ISSUE maybe allow different mul type (not value_type)
        __DEVICE_TAG__ mat_expr_scalar_mul(const E &_e, const value_type &_mul) : e(_e), mul(_mul)
	{
		//printf("mat_expr_negate created\n");
        }
        /*__DEVICE_TAG__ ~mat_expr_negate()
        {
		printf("mat_expr_negate deleted\n");
	}*/
	__DEVICE_TAG__ value_type	operator()(int i1,int i2)const
	{
		return e(i1,i2)*mul;
	}
private:
	value_type	mul;
	const E		&e;
};

//E1,E2 are VECTOR_EXPRESSION
template<class E1,class E2>
class mat_expr_sum : public mat_expr_base< mat_expr_sum<E1,E2> >
{
public:
	//TODO check size coincidence (statically)
	static const int		sz1 = E1::sz1;
	static const int		sz2 = E1::sz2;
	//TODO somehow deduce type of sum E1::value_type and E2::value_type
	typedef typename E1::value_type	value_type;

        __DEVICE_TAG__ mat_expr_sum(const E1 &_e1,const E2 &_e2) : e1(_e1), e2(_e2)
	{
        }
	__DEVICE_TAG__ value_type	operator()(int i1,int i2)const
	{
		return e1(i1,i2)+e2(i1,i2);
	}
private:
	const E1	&e1;
	const E2	&e2;
};

//E1,E2 are VECTOR_EXPRESSION
template<class E1,class E2>
class mat_expr_diff : public mat_expr_base< mat_expr_diff<E1,E2> >
{
public:
	//TODO check size coincidence (statically)
	static const int		sz1 = E1::sz1;
	static const int		sz2 = E1::sz2;
	//TODO somehow deduce type of sum E1::value_type and E2::value_type
	typedef typename E1::value_type	value_type;

        __DEVICE_TAG__ mat_expr_diff(const E1 &_e1,const E2 &_e2) : e1(_e1), e2(_e2)
	{
        }
	__DEVICE_TAG__ value_type	operator()(int i1,int i2)const
	{
		return e1(i1,i2)-e2(i1,i2);
	}
private:
	const E1	&e1;
	const E2	&e2;
};

//E1,E2 are VECTOR_EXPRESSION
template<class E1,class E2>
class mat_expr_prod : public mat_expr_base< mat_expr_prod<E1,E2> >
{
public:
	//TODO check size coincidence (statically)
	static const int		sz1 = E1::sz1;
	static const int		sz2 = E2::sz2;
	//TODO somehow deduce type of sum E1::value_type and E2::value_type
	typedef typename E1::value_type	value_type;

        __DEVICE_TAG__ mat_expr_prod(const E1 &_e1,const E2 &_e2) : e1(_e1), e2(_e2)
	{
        }
	__DEVICE_TAG__ value_type	operator()(int i1,int i2)const
	{
		value_type	res(0.f);
		for (int j = 0;j < E1::sz2;++j) {
			res += e1(i1,j)*e2(j,i2);
		}
		return res;
	}
private:
	const E1	&e1;
	const E2	&e2;
};

template<class E>
__DEVICE_TAG__ mat_expr_negate<E>	operator-(const mat_expr_base<E> &e)
{
	return mat_expr_negate<E>(e());
}

template<class E>
__DEVICE_TAG__ mat_expr_scalar_mul<E>	operator*(const mat_expr_base<E> &e, const typename E::value_type &mul)
{
	return mat_expr_scalar_mul<E>(e(), mul);
}

template<class E>
__DEVICE_TAG__ mat_expr_scalar_mul<E>	operator*(const typename E::value_type &mul, const mat_expr_base<E> &e)
{
	return mat_expr_scalar_mul<E>(e(), mul);
}

template<class E1,class E2>
__DEVICE_TAG__ mat_expr_sum<E1,E2>	operator+(const mat_expr_base<E1> &e1, const mat_expr_base<E2> &e2)
{
	return mat_expr_sum<E1,E2>(e1(),e2());
}

template<class E1,class E2>
__DEVICE_TAG__ mat_expr_diff<E1,E2>	operator-(const mat_expr_base<E1> &e1, const mat_expr_base<E2> &e2)
{
	return mat_expr_diff<E1,E2>(e1(),e2());
}

template<class E1,class E2>
__DEVICE_TAG__ mat_expr_prod<E1,E2>	operator*(const mat_expr_base<E1> &e1, const mat_expr_base<E2> &e2)
{
	return mat_expr_prod<E1,E2>(e1(),e2());
}

#endif
