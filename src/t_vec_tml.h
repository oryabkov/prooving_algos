#ifndef __T_VEC_TML_H__
#define __T_VEC_TML_H__

template<int dim, class T>
struct t_vec_tml
{
	T x[dim];
	
	t_vec_tml(const T &x_ = T(0.f), const T &y_ = T(0.f), const T &z_ = T(0.f))
	{
		if (dim > 0) x[0] = x_;
		if (dim > 1) x[1] = y_;
		if (dim > 2) x[2] = z_;
	}

	T	&operator[](int j)
	{
		return x[j];
	}
	const T	&operator[](int j)const
	{
		return x[j];
	}
};

#endif
