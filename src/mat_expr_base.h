#ifndef __MAT_EXPR_BASE_H__
#define __MAT_EXPR_BASE_H__

#ifndef __CUDACC__
#define __DEVICE_TAG__
#else
#define __DEVICE_TAG__ __device__ __host__
#endif

//MATRIX_EXPRESSION concept E
//public base mat_expr_base<E>
//static const int sz1 - rows num
//static const int sz2 - cols num
//nested type value_type
//value_type operator()(int i1,int i2)
//OR
//const value_type &operator()(int i1,int i2)

template<class E>
class mat_expr_base
{
public:
	typedef	E	expr_type;

	__DEVICE_TAG__ E	&operator()()
	{
		return static_cast<E&>(*this);
	}
	__DEVICE_TAG__ const E	&operator()()const
	{
		return static_cast<const E&>(*this);
	}
};

#endif
