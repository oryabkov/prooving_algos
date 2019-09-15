
#include <cstdio>
#include <boost/numeric/interval.hpp>
#include "t_vec_tml.h"
#include "TModel.h"
#include "TModelPrint.h"
#include "TModelPicardIteration.h"

#define TAYLOR_ORDER 4
//#define TAYLOR_ORDER 2

template<class T>
struct t_van_der_pol_ode
{
	typedef t_vec_tml<2,T>		vec_type;
	static const int		dim = 2;
	double				mu;

	//void verify_size(vec_type& x)const {}
	void operator()(const vec_type& x, const T &t, vec_type& f)const
	{
		f[0] = x[1];
		//f.y = mu*(t_int(1.f)-x.x*x.x)*x.y - x.x;
		//pow
		f[1] = (T(1.f)-x[0]*x[0])*x[1]*mu - x[0];
	}
	/*void addMul(vec_type& B, t_int mulB, const vec_type& A, t_int mulA)const
	{
		B.x = B.x*mulB + A.x*mulA;
		B.y = B.y*mulB + A.y*mulA;
	}
	void addMul(vec_type& C, t_int mulC, const vec_type& A, t_int mulA, const vec_type& B, t_int mulB)const
	{
		C.x = C.x*mulC + A.x*mulA + B.x*mulB;
		C.y = C.y*mulC + A.y*mulA + B.y*mulB;
	}*/
};

typedef TModel<3,TAYLOR_ORDER,double>					t_TM;
typedef t_vec_tml<2,t_TM>						t_TM_vec;
typedef t_van_der_pol_ode<t_TM>						t_TM_ode;
typedef boost::numeric::interval<double>				t_interval;
typedef TModelPicardIteration<t_TM_ode,t_TM_vec,t_TM,t_interval,double>	t_picard_iteration;

//old output style (parallelogramm + remaider rect)
/*
void	print_res_line(FILE *stream, double u1, double v1, double u2, double v2)
{
	fprintf( stream, "SL(%f, %f, %f, %f, %f, %f)",  u1, v1, 0.f,
							u2, v2, 0.f);
	fprintf( stream,"{");
	fprintf(stream, "%f,",1.f);
	fprintf(stream, "%f",1.f);
	fprintf(stream, "};\n");
}

void	print_res_line(FILE *stream, const t_TM_vec &x, double u1, double v1, double u2, double v2)
{
	double	p1[3] = {0.,u1,v1},
		p2[3] = {0.,u2,v2};
	fprintf( stream, "SL(%f, %f, %f, %f, %f, %f)",  x[0].polynomial().eval(p1), x[1].polynomial().eval(p1), 0.f,
							x[0].polynomial().eval(p2), x[1].polynomial().eval(p2), 0.f);
	fprintf( stream,"{");
	fprintf(stream, "%f,",1.f);
	fprintf(stream, "%f",1.f);
	fprintf(stream, "};\n");
}

void	print_res_pnt(FILE *stream, const t_TM_vec &x, double u1, double v1, double err_u1, double err_v1, double err_u2, double err_v2)
{
	double	p1[3] = {0.,u1,v1};
	fprintf( stream, "SP(%f, %f, %f)",  x[0].polynomial().eval(p1), x[1].polynomial().eval(p1), 0.f);
	fprintf( stream,"{");
	fprintf(stream, "%f",1.f);
	fprintf(stream, "};\n");

	print_res_line(stream, x[0].polynomial().eval(p1)+err_u1, x[1].polynomial().eval(p1)+err_v1, x[0].polynomial().eval(p1)+err_u2, x[1].polynomial().eval(p1)+err_v1);
	print_res_line(stream, x[0].polynomial().eval(p1)+err_u2, x[1].polynomial().eval(p1)+err_v1, x[0].polynomial().eval(p1)+err_u2, x[1].polynomial().eval(p1)+err_v2);
	print_res_line(stream, x[0].polynomial().eval(p1)+err_u2, x[1].polynomial().eval(p1)+err_v2, x[0].polynomial().eval(p1)+err_u1, x[1].polynomial().eval(p1)+err_v2);
	print_res_line(stream, x[0].polynomial().eval(p1)+err_u1, x[1].polynomial().eval(p1)+err_v2, x[0].polynomial().eval(p1)+err_u1, x[1].polynomial().eval(p1)+err_v1);
}

void	print_res(FILE *stream, const t_TM_vec &x, int view_num)
{
	fprintf( stream, "View");
	fprintf( stream, " '");
	fprintf( stream, "%d", view_num);
	fprintf( stream, "' {\n");
	fprintf( stream, "TIME{0};\n");

	print_res_line(stream, x, -1., -1.,  1., -1.);
	print_res_line(stream, x,  1., -1.,  1.,  1.);
	print_res_line(stream, x,  1.,  1., -1.,  1.);
	print_res_line(stream, x, -1.,  1., -1., -1.);
	
	print_res_pnt(stream, x, 0., 0., x[0].remainder().lower(), x[1].remainder().lower(), x[0].remainder().upper(), x[1].remainder().upper());

	fprintf( stream, "};\n");
}
*/

//new output style (something similar to convex hull)

void    print_res_line(FILE *stream, double u1, double v1, double u2, double v2)
{
        fprintf( stream, "SL(%0.15e, %0.15e, %0.15e, %0.15e, %0.15e, %0.15e)",  u1, v1, 0.f,
                                                                                u2, v2, 0.f);
        fprintf( stream,"{");
        fprintf(stream, "%f,",1.f);
        fprintf(stream, "%f",1.f);
        fprintf(stream, "};\n");
}

void    print_res_curve(FILE *stream, double *u, double *v, int steps_n = 1)
{
        for (int i = 0;i < steps_n;++i) {
                /*double  u1_ = u1 + i*(u2-u1)/steps_n, 
                        v1_ = v1 + i*(v2-v1)/steps_n, 
                        u2_ = u1 + (i+1)*(u2-u1)/steps_n, 
                        v2_ = v1 + (i+1)*(v2-v1)/steps_n;*/
                double  u1_ = u[i], 
                        v1_ = v[i], 
                        u2_ = u[i+1],
                        v2_ = v[i+1];
                fprintf( stream, "SL(%0.15e, %0.15e, %0.15e, %0.15e, %0.15e, %0.15e)",  u1_, v1_, 0.f,
                                                                                        u2_, v2_, 0.f);
                fprintf( stream,"{");
                fprintf(stream, "%f,",1.f);
                fprintf(stream, "%f",1.f);
                fprintf(stream, "};\n");
        }
}

void    print_res_rect(FILE *stream, double u1, double v1, double u2, double v2)
{
        print_res_line(stream, u1, v1, u2, v1);
        print_res_line(stream, u2, v1, u2, v2);
        print_res_line(stream, u2, v2, u1, v2);
        print_res_line(stream, u1, v2, u1, v1);
}

void    print_res_connect_rects(FILE *stream, double u11, double v11, double u12, double v12, double u21, double v21, double u22, double v22)
{
        print_res_line(stream, u11, v11, u21, v21);
        print_res_line(stream, u11, v12, u21, v22);
        print_res_line(stream, u12, v11, u22, v21);
        print_res_line(stream, u12, v12, u22, v22);
}

void    print_res_connect_rects(FILE *stream, double *u1, double *v1, double *u2, double *v2, int steps_n)
{
        print_res_curve(stream, u1, v1, steps_n);
        print_res_curve(stream, u1, v2, steps_n);
        print_res_curve(stream, u2, v1, steps_n);
        print_res_curve(stream, u2, v2, steps_n);
}

//void    print_res(const t_TM_ode &ode, FILE *stream, t_TM_space &x, int view_num)
void    print_res(FILE *stream, const t_TM_vec &x, int view_num, bool type)
{
        fprintf( stream, "View");
        fprintf( stream, " '");
        fprintf( stream, "%d", view_num);
        fprintf( stream, "' {\n");
        fprintf( stream, "TIME{0};\n");

        t_TM tmp;
        
        /*t_interval bnds1[2][2],bnds2[2][2];

        for (int s1 = -1;s1 <= 1;s1 += 2)
        for (int s2 = -1;s2 <= 1;s2 += 2) {
                tmp = x[0];
                tmp = tmp.part_eval(0, double(s1));
                tmp = tmp.part_eval(1, double(s2));
                t_interval bnd1;
                if (type)
                        bnd1 = tmp.bnds() + x[0].remainder();
                else
                        bnd1 = tmp.bnds() + tmp.remainder();
                bnds1[(s1+1)/2][(s2+1)/2] = bnd1;

                tmp = x[1];
                tmp = tmp.part_eval(0, double(s1));
                tmp = tmp.part_eval(1, double(s2));
                t_interval bnd2;
                if (type)
                        bnd2 = tmp.bnds() + x[1].remainder();
                else 
                        bnd2 = tmp.bnds() + tmp.remainder();
                bnds2[(s1+1)/2][(s2+1)/2] = bnd2;

                print_res_rect(stream, bnd1.lower(), bnd2.lower(), bnd1.upper(), bnd2.upper());
        }

        print_res_connect_rects(stream, bnds1[0][0].lower(), bnds2[0][0].lower(), bnds1[0][0].upper(), bnds2[0][0].upper(), bnds1[0][1].lower(), bnds2[0][1].lower(), bnds1[0][1].upper(), bnds2[0][1].upper());
        print_res_connect_rects(stream, bnds1[0][1].lower(), bnds2[0][1].lower(), bnds1[0][1].upper(), bnds2[0][1].upper(), bnds1[1][1].lower(), bnds2[1][1].lower(), bnds1[1][1].upper(), bnds2[1][1].upper());
        print_res_connect_rects(stream, bnds1[1][1].lower(), bnds2[1][1].lower(), bnds1[1][1].upper(), bnds2[1][1].upper(), bnds1[1][0].lower(), bnds2[1][0].lower(), bnds1[1][0].upper(), bnds2[1][0].upper());
        print_res_connect_rects(stream, bnds1[1][0].lower(), bnds2[1][0].lower(), bnds1[1][0].upper(), bnds2[1][0].upper(), bnds1[0][0].lower(), bnds2[0][0].lower(), bnds1[0][0].upper(), bnds2[0][0].upper());*/

        int steps_n = 20;

        double  *u1 = new double [steps_n+1], 
                *v1 = new double [steps_n+1], 
                *u2 = new double [steps_n+1], 
                *v2 = new double [steps_n+1];

        for (int l = 0;l < 4;++l) {
                double u_begin, v_begin, u_end, v_end;
                switch (l) {
                        case 0:
                                u_begin =-1.; v_begin =-1.; u_end = 1.; v_end =-1.;
                        break;
                        case 1:
                                u_begin = 1.; v_begin =-1.; u_end = 1.; v_end = 1.;
                        break;
                        case 2:
                                u_begin = 1.; v_begin = 1.; u_end =-1.; v_end = 1.;
                        break;
                        case 3:
                                u_begin =-1.; v_begin = 1.; u_end =-1.; v_end =-1.;
                        break;
                }
                for (int s = 0;s < steps_n+1;++s) {
                        double  s1 = u_begin + s*(u_end-u_begin)/steps_n, 
                                s2 = v_begin + s*(v_end-v_begin)/steps_n;
                        tmp = x[0];
                        tmp = tmp.part_eval(1, s1);
                        tmp = tmp.part_eval(2, s2);
                        t_interval bnd1;
                        if (type)
                                bnd1 = tmp.bnds() + x[0].remainder();
                        else
                                bnd1 = tmp.bnds() + tmp.remainder();
                        //printf("%f %f\n", bnd1.lower(), bnd1.upper());
                        //bnds1[(s1+1)/2][(s2+1)/2] = bnd1;

                        tmp = x[1];
                        tmp = tmp.part_eval(1, s1);
                        tmp = tmp.part_eval(2, s2);
                        t_interval bnd2;
                        if (type)
                                bnd2 = tmp.bnds() + x[1].remainder();
                        else 
                                bnd2 = tmp.bnds() + tmp.remainder();
                        //printf("%f %f\n", bnd2.lower(), bnd2.upper());
                        //bnds2[(s1+1)/2][(s2+1)/2] = bnd2;

                        if ((s == 0)||(s == steps_n))
                                print_res_rect(stream, bnd1.lower(), bnd2.lower(), bnd1.upper(), bnd2.upper());

                        u1[s] = bnd1.lower(); 
                        v1[s] = bnd2.lower();
                        u2[s] = bnd1.upper();
                        v2[s] = bnd2.upper();
                }
                print_res_connect_rects(stream, u1, v1, u2, v2, steps_n);
        }

        delete []u1;
        delete []v1;
        delete []u2;
        delete []v2;

        fprintf( stream, "};\n");
}

int main()
{
	t_TM::init_exponents_table();

	t_TM_vec		x, x_prev;
        t_TM_ode		ode;
        t_picard_iteration      picard_iteration;
        
        ode.mu = 1.;
        //x.x[0] = t_int(-2.0085,-2.0085); st.y = 0.;

	x[0].polynomial()[0] = -2.0085;
	//x[0].polynomial()[0] = -2.0100;
	int	exp1[3] = {0,1,0};
	x[0].polynomial()[t_TM::polynomial_type::term_index(exp1)] = 0.001;
	//x[0].polynomial()[t_TM::polynomial_type::term_index(exp1)] = 0.005;
	x[1].polynomial()[0] = 0.;
	int	exp2[3] = {0,0,1};
	x[1].polynomial()[t_TM::polynomial_type::term_index(exp2)] = 0.0001;

        t_interval		t(0.);
        double			dt = 0.0005;

        FILE    *f = fopen("van_der_pol.pos", "w");

	int	views_cnt = 0;

	for (int i = 0;i < ceil(7./dt);++i) {
		x_prev = x;

		//if (i%100 == 0) {
		if (i%10 == 0) {
			print_res(f, x, i, true);
                        //print_res(f, x, i, false);
			views_cnt++;
                        //views_cnt += 2;
		}

        	picard_iteration.doStep(ode, x, t, dt);
        	for (int j = 0;j < 2;++j)
	        	x[j] = x[j].part_eval(0,1.);

		/*if ((x_prev[1].polynomial()[0] < 0.)&&(x[1].polynomial()[0] > 0.)) {	//aapproximate section detection
			print_res(f, x_prev, i-1, false);
			print_res(f, x, i, false);
			views_cnt += 2;
		}*/

                /*double	sc[2+1];
		sc[0] = 0.1;
		sc[1] = 0.001;
		sc[2] = 0.001;
		printf("t = %f:\n", t.lower());
		for (int j = 0;j < 2;++j) {
			printf("TModel %d:\n", j);
			TModel_print_sc<2+1,3,double>(x[j],sc);
		}*/

		if (i%100 == 0) printf("t = %f\n", t.lower());
	}

        printf("final reminder:\n");
        for (int j = 0;j < 2;++j) {
                printf("j = %d: [%e,%e]\n", j, x[j].remainder().lower(), x[j].remainder().upper());
        }

	for (int i = 0;i < views_cnt;++i) {
		fprintf(f, "View[%d].ShowScale = 0;\n", i);
                fprintf(f, "View[%d].RangeType = 2;\n", i);
                fprintf(f, "View[%d].CustomMin = 1.;\n", i);
                fprintf(f, "View[%d].CustomMax = 2.;\n", i);
	}

 	fclose(f);
 	
 	t_TM::free_exponents_table();

	return 0;
}