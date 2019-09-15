#ifndef __T_MAT_TML_H__
#define __T_MAT_TML_H__

#include "mat_expr_base.h"    
#include "mat_expressions.h"

#ifndef __CUDACC__
#define __DEVICE_TAG__
#else
#define __DEVICE_TAG__ __device__ __host__
#endif

template<class T,int dim1,int dim2>
class __no_alias_t_mat_tml;

template<class T,int dim1,int dim2>
struct t_mat_tml : public mat_expr_base<t_mat_tml<T,dim1,dim2> >
{
	//was used before
	//TODO delete it
	typedef T		base_type;

	//MATRIX_EXPRESSION concept
	typedef T		value_type;
	static const int	sz1 = dim1;
	static const int	sz2 = dim2;

	T			d[dim1][dim2];
	
	__DEVICE_TAG__  t_mat_tml()
	{
	}
	//TODO check size (statically)
	template<class E>
	__DEVICE_TAG__  t_mat_tml(const mat_expr_base<E> &e)
	{
		for (int i = 0;i < dim1;++i)
		for (int j = 0;j < dim2;++j)
			d[i][j] = T(e()(i,j));
		//return m;
	}

	__DEVICE_TAG__ T	&operator()(int i1,int i2) { return d[i1][i2]; }
	__DEVICE_TAG__ const T	&operator()(int i1,int i2)const { return d[i1][i2]; }
	__DEVICE_TAG__ int	size1()const { return dim1; }
	__DEVICE_TAG__ int	size2()const { return dim2; }

	/*__DEVICE_TAG__ t_mat_tml	&operator+=(const t_mat_tml &x)
        {
		for (int i = 0;i < dim1;++i)
		for (int j = 0;j < dim2;++j) {
			d[i][j] += x.d[i][j];
		}
		return *this;
	}
	__DEVICE_TAG__ t_mat_tml	&operator-=(const t_mat_tml &x)
        {
		for (int i = 0;i < dim1;++i)
		for (int j = 0;j < dim2;++j) {
			d[i][j] -= x.d[i][j];
		}
		return *this;
	}
	__DEVICE_TAG__ t_mat_tml	&operator*=(const T &x)
        {
		for (int i = 0;i < dim1;++i)
		for (int j = 0;j < dim2;++j) {
			d[i][j] *= x;
		}
		return *this;
	}
	__DEVICE_TAG__ t_mat_tml	&operator/=(const T &x)
        {
		for (int i = 0;i < dim1;++i)
		for (int j = 0;j < dim2;++j) {
			d[i][j] /= x;
		}
		return *this;
	}*/

	//__DEVICE_TAG__ int	size1()const { return dim1; }
	//__DEVICE_TAG__ int	size2()const { return dim2; }
	
	//TODO check size (statically)
	template<class E>
	__DEVICE_TAG__ t_mat_tml &operator=(const mat_expr_base<E> &e)
	{
		t_mat_tml<T,dim1,dim2>	res;
		for (int i = 0;i < dim1;++i)
		for (int j = 0;j < dim2;++j)
			res.d[i][j] = T(e()(i,j));
		for (int i = 0;i < dim1;++i)
		for (int j = 0;j < dim2;++j)
			d[i][j] = res.d[i][j];
		return *this;
	}
	template<class E>
	__DEVICE_TAG__ t_mat_tml &operator+=(const mat_expr_base<E> &e)
	{
		t_mat_tml<T,dim1,dim2>	res;
		for (int i = 0;i < dim1;++i)
		for (int j = 0;j < dim2;++j)
			res.d[i][j] = T(e()(i,j));
		for (int i = 0;i < dim1;++i)
		for (int j = 0;j < dim2;++j)
			d[i][j] += res.d[i][j];
		return *this;
	}
	template<class E>
	__DEVICE_TAG__ t_mat_tml &operator-=(const mat_expr_base<E> &e)
	{
		t_mat_tml<T,dim1,dim2>	res;
		for (int i = 0;i < dim1;++i)
		for (int j = 0;j < dim2;++j)
			res.d[i][j] = T(e()(i,j));
		for (int i = 0;i < dim1;++i)
		for (int j = 0;j < dim2;++j)
			d[i][j] -= res.d[i][j];
		return *this;
	}
	__DEVICE_TAG__ t_mat_tml &operator*=(const value_type &mul)
	{
		for (int i = 0;i < dim1;++i)
		for (int j = 0;j < dim2;++j)
			d[i][j] *= mul;
		return *this;
	}
	__DEVICE_TAG__ t_mat_tml &operator/=(const value_type &mul)
	{
		for (int i = 0;i < dim1;++i)
		for (int j = 0;j < dim2;++j)
			d[i][j] /= mul;
		return *this;
	}
	//TODO check size (statically)
	template<class E>
	__DEVICE_TAG__ t_mat_tml &operator*=(const mat_expr_base<E> &e)
	{
		t_mat_tml<T,E::sz1,E::sz2>	res;
		for (int i = 0;i < E::sz1;++i)
		for (int j = 0;j < E::sz2;++j)
			res.d[i][j] = T(e()(i,j));

		T	row_buf[dim2];
		for (int i = 0;i < dim1;++i) {
			for (int j = 0;j < dim2;++j) row_buf[j] = d[i][j];

			for (int j = 0;j < dim2;++j) {
				d[i][j] = T(0.f);
				for (int k = 0;k < dim2;++k)
					d[i][j] += row_buf[k]*res.d[k][j];
			}
		}
		return *this;
	}

        __DEVICE_TAG__ __no_alias_t_mat_tml<T,dim1,dim2> noalias();
};

template<class T,int dim1,int dim2>
class __no_alias_t_mat_tml
{
	t_mat_tml<T,dim1,dim2>	&m;
public:
	__DEVICE_TAG__ explicit __no_alias_t_mat_tml(t_mat_tml<T,dim1,dim2> &_m) : m(_m)
	{
	}
	//TODO check size (statically)
	template<class E>
	__DEVICE_TAG__ t_mat_tml<T,dim1,dim2> &operator=(const mat_expr_base<E> &e)
	{
		for (int i = 0;i < dim1;++i)
		for (int j = 0;j < dim2;++j)
			m.d[i][j] = T(e()(i,j));
		return m;
	}
	//TODO check size (statically)
	template<class E>
	__DEVICE_TAG__ t_mat_tml<T,dim1,dim2> &operator+=(const mat_expr_base<E> &e)
	{
		for (int i = 0;i < dim1;++i)
		for (int j = 0;j < dim2;++j)
			m.d[i][j] += T(e()(i,j));
		return m;
	}
	//TODO check size (statically)
	template<class E>
	__DEVICE_TAG__ t_mat_tml<T,dim1,dim2> &operator-=(const mat_expr_base<E> &e)
	{
		for (int i = 0;i < dim1;++i)
		for (int j = 0;j < dim2;++j)
			m.d[i][j] -= T(e()(i,j));
		return m;
	}
	__DEVICE_TAG__ t_mat_tml<T,dim1,dim2> &operator*=(const T &mul)
	{
		for (int i = 0;i < dim1;++i)
		for (int j = 0;j < dim2;++j)
			m.d[i][j] *= mul;
		return m;
	}
	__DEVICE_TAG__ t_mat_tml<T,dim1,dim2> &operator/=(const T &mul)
	{
		for (int i = 0;i < dim1;++i)
		for (int j = 0;j < dim2;++j)
			m.d[i][j] /= mul;
		return m;
	}
	//TODO check size (statically)
	template<class E>
	__DEVICE_TAG__ t_mat_tml<T,dim1,dim2> &operator*=(const mat_expr_base<E> &e)
	{
		T	row_buf[dim2];
		for (int i = 0;i < dim1;++i) {
			for (int j = 0;j < dim2;++j) row_buf[j] = m.d[i][j];

			for (int j = 0;j < dim2;++j) {
				m.d[i][j] = T(0.f);
				for (int k = 0;k < dim2;++k)
					m.d[i][j] += row_buf[k]*T(e()(k,j));
			}
		}
		return m;
	}
};

template<class T,int dim1,int dim2>
__no_alias_t_mat_tml<T,dim1,dim2> t_mat_tml<T,dim1,dim2>::noalias()
{
	return __no_alias_t_mat_tml<T,dim1,dim2>(*this);
}

template<class T,int dim1,int dim2>
struct t_mat_zero_tml : public mat_expr_base<t_mat_zero_tml<T,dim1,dim2> >
{
	//MATRIX_EXPRESSION concept
	typedef T		value_type;
	static const int	sz1 = dim1;
	static const int	sz2 = dim2;
	
	__DEVICE_TAG__ t_mat_zero_tml()
        {
	}
	//NOTE this strange help-constructor was added because 4 some strange reason
	//some times some compilers doesnot recognize initializations like this
	//t_mat_tml<float,5,5>	m(t_mat_zero_tml<float,5,5>());
	//and percipe this expression as some function signature (??)
	//instead we can use
	//t_mat_tml<float,5,5>	m(t_mat_zero_tml<float,5,5>(0));
	__DEVICE_TAG__ t_mat_zero_tml(const int &)
	{
	}

	__DEVICE_TAG__ T	operator()(int i1,int i2)const { return T(0.f); }
};

template<class T,int dim1,int dim2>
struct t_mat_ident_tml : public mat_expr_base<t_mat_ident_tml<T,dim1,dim2> >
{
	//MATRIX_EXPRESSION concept
	typedef T		value_type;
	static const int	sz1 = dim1;
	static const int	sz2 = dim2;
	
	__DEVICE_TAG__ t_mat_ident_tml()
        {
	}
	//NOTE see NOTE in t_mat_zero_tml
	__DEVICE_TAG__ t_mat_ident_tml(const int &)
	{
	}

	__DEVICE_TAG__ T	operator()(int i1,int i2)const { return (i1==i2?T(1.f):T(0.f)); }
};

template<class T,int dim1,int dim2>
struct t_mat_scalar_tml : public mat_expr_base<t_mat_scalar_tml<T,dim1,dim2> >
{
	//MATRIX_EXPRESSION concept
	typedef T		value_type;
	static const int	sz1 = dim1;
	static const int	sz2 = dim2;

	T	val;
        __DEVICE_TAG__ t_mat_scalar_tml()
        {
	}
	__DEVICE_TAG__ t_mat_scalar_tml(const T &_val) : val(_val)
	{
	}

	__DEVICE_TAG__ T	operator()(int i1,int i2)const { return (i1==i2?val:T(0.f)); }
};

template<class T,int dim1,int dim2,int dim3>
__DEVICE_TAG__ void mat_prod(const t_mat_tml<T,dim1,dim2> &m1, const t_mat_tml<T,dim2,dim3> &m2, t_mat_tml<T,dim1,dim3> &m_res)
{
	for (int i = 0;i < dim1;++i)
	for (int j = 0;j < dim3;++j) {
		m_res(i,j) = T(0.f);
		for (int k = 0;k < dim2;++k) {
			m_res(i,j) += m1(i,k)*m2(k,j);
		}
	}

}


#endif
