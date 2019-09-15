
#include <cstdio>
#include <vector>
#include <boost/numeric/interval.hpp>
#include "t_vec_tml.h"
#include "t_mat_tml.h"
#include "TModel.h"
#include "TModelPrint.h"
#include "TModelPicardIteration.h"

#define TAYLOR_ORDER 2

template<class T>
struct t_van_der_pol_ode
{
        typedef t_vec_tml<2,T>          vec_type;
        static const int                dim = 2;
        double                          mu;

        //void verify_size(vec_type& x)const {}
        void operator()(const vec_type& x, const T &t, vec_type& f)const
        {
                f[0] = x[1];
                f[1] = (T(1.f)-x[0]*x[0])*x[1]*mu - x[0];
        }
};

typedef t_vec_tml<2,double>                                             t_vec;
typedef t_mat_tml<double,2,2>                                           t_mat;
typedef TModel<3,TAYLOR_ORDER,double>                                   t_TM;
typedef t_vec_tml<2,t_TM>                                               t_TM_vec;
typedef t_van_der_pol_ode<t_TM>                                         t_TM_ode;
typedef boost::numeric::interval<double>                                t_interval;
typedef TModelPicardIteration<t_TM_ode,t_TM_vec,t_TM,t_interval,double> t_picard_iteration;

void    print_res_line(FILE *stream, double u1, double v1, double u2, double v2, int c)
{
        fprintf( stream, "SL(%f, %f, %f, %f, %f, %f)",  u1, v1, 0.f,
                                                        u2, v2, 0.f);
        fprintf( stream,"{");
        //fprintf(stream, "%f,",1.f);
        if (c == 0) fprintf(stream, "%f,",1.f); else fprintf(stream, "%f,",2.f);
        if (c == 0) fprintf(stream, "%f",1.f); else fprintf(stream, "%f",2.f);
        fprintf(stream, "};\n");
}

void    print_res_line(FILE *stream, const t_TM_vec &x, double u1, double v1, double u2, double v2, int c)
{
        double  p1[3] = {0.,u1,v1},
                p2[3] = {0.,u2,v2};
        fprintf( stream, "SL(%f, %f, %f, %f, %f, %f)",  x[0].polynomial().eval(p1), x[1].polynomial().eval(p1), 0.f,
                                                        x[0].polynomial().eval(p2), x[1].polynomial().eval(p2), 0.f);
        fprintf( stream,"{");
        if (c == 0) fprintf(stream, "%f,",1.f); else fprintf(stream, "%f,",2.f);
        if (c == 0) fprintf(stream, "%f",1.f); else fprintf(stream, "%f",2.f);
        fprintf(stream, "};\n");
}

void    print_res_pnt(FILE *stream, const t_TM_vec &x, double u1, double v1, double err_u1, double err_v1, double err_u2, double err_v2, int c)
{
        double  p1[3] = {0.,u1,v1};
        fprintf( stream, "SP(%f, %f, %f)",  x[0].polynomial().eval(p1), x[1].polynomial().eval(p1), 0.f);
        fprintf( stream,"{");
        fprintf(stream, "%f",1.f);
        fprintf(stream, "};\n");

        print_res_line(stream, x[0].polynomial().eval(p1)+err_u1, x[1].polynomial().eval(p1)+err_v1, x[0].polynomial().eval(p1)+err_u2, x[1].polynomial().eval(p1)+err_v1, c);
        print_res_line(stream, x[0].polynomial().eval(p1)+err_u2, x[1].polynomial().eval(p1)+err_v1, x[0].polynomial().eval(p1)+err_u2, x[1].polynomial().eval(p1)+err_v2, c);
        print_res_line(stream, x[0].polynomial().eval(p1)+err_u2, x[1].polynomial().eval(p1)+err_v2, x[0].polynomial().eval(p1)+err_u1, x[1].polynomial().eval(p1)+err_v2, c);
        print_res_line(stream, x[0].polynomial().eval(p1)+err_u1, x[1].polynomial().eval(p1)+err_v2, x[0].polynomial().eval(p1)+err_u1, x[1].polynomial().eval(p1)+err_v1, c);
}

void    print_res(FILE *stream, const t_TM_vec &x, int view_num, int c)
{
        fprintf( stream, "View");
        fprintf( stream, " '");
        fprintf( stream, "%d", view_num);
        fprintf( stream, "' {\n");
        fprintf( stream, "TIME{0};\n");

        print_res_line(stream, x, -1., -1.,  1., -1., c);
        print_res_line(stream, x,  1., -1.,  1.,  1., c);
        print_res_line(stream, x,  1.,  1., -1.,  1., c);
        print_res_line(stream, x, -1.,  1., -1., -1., c);

        print_res_pnt(stream, x, 0., 0., x[0].remainder().lower(), x[1].remainder().lower(), x[0].remainder().upper(), x[1].remainder().upper(), c);

        fprintf( stream, "};\n");
}

std::vector<double>     cycle_times;
std::vector<t_vec>      cycle_pnts;
std::vector<t_mat>      transf_mats;

void    read_cycle_pnts()
{
        //read approximate cycle file
        FILE    *f_in = fopen("test.dat","r");
        double  t_,x_,y_;
        //double        t_prev,x_prev,y_prev;
        int     cnt = 0;
        while (fscanf(f_in,"%lf %lf %lf", &t_, &x_, &y_) == 3) {
                /*if ((cnt > 0)&&(cnt%dn == 0)) {
                        double  tang_x = (x_-x_prev)/sqrt(pow((x_-x_prev),2)+pow((y_-y_prev),2)),
                                tang_y = (y_-y_prev)/sqrt(pow((x_-x_prev),2)+pow((y_-y_prev),2)),
                                tran_x = tang_y,
                                tran_y = -tang_x;

                        t_TM_vec        x;
                        x[0].polynomial()[0] = x_;
                        x[1].polynomial()[0] = y_;
                        int     exp1[3] = {0,1,0};
                        int     exp2[3] = {0,0,1};
                        x[0].polynomial()[t_TM::polynomial_type::term_index(exp1)] = tang_x*r;
                        x[0].polynomial()[t_TM::polynomial_type::term_index(exp2)] = tran_x*r;
                        x[1].polynomial()[t_TM::polynomial_type::term_index(exp1)] = tang_y*r;
                        x[1].polynomial()[t_TM::polynomial_type::term_index(exp2)] = tran_y*r;

                        //integrate(f, x_, y_, r, 0.01*dn*30, views_cnt);
                        integrate(f, x, 0.01*dn*10, views_cnt);
                } */
                cycle_times.push_back(t_);
                cycle_pnts.push_back(t_vec(x_,y_));
                cnt++;
                //t_prev = t_; x_prev = x_; y_prev = y_;
                //if (cnt > 500)  break;
        }
        fclose(f_in);

        //leave one whole period
        /*double        min_dist2 = -1.;
        int     min_dist_i;
        for (int i = 100;i < cycle_pnts.size();++i) {
                double dist2 = pow(cycle_pnts[i][0]-cycle_pnts[0][0],2) + pow(cycle_pnts[i][1]-cycle_pnts[0][1],2);
                if ((min_dist2 < 0.)||(dist2 < min_dist2)) {
                        min_dist2 = dist2;
                        min_dist_i = i;
                }
        }*/
        double  eps = 1e-2;
        int     min_dist_i;
        for (int i = 5;i < cycle_pnts.size();++i) {
                double dist2 = pow(cycle_pnts[i][0]-cycle_pnts[0][0],2) + pow(cycle_pnts[i][1]-cycle_pnts[0][1],2);
                if (dist2 < eps*eps) {
                        min_dist_i = i;
                        break;
                }
        }
        cycle_pnts.erase(cycle_pnts.begin() + min_dist_i, cycle_pnts.end());
        cycle_times.erase(cycle_times.begin() + min_dist_i, cycle_times.end());
}

double  det22(const t_mat &m)
{
        return m(0,0)*m(1,1)-m(0,1)*m(1,0);
}

void    build_loc_coords_stupid()
{
        for (int i = 0;i < cycle_pnts.size();++i) {
                t_mat           transf_mat;
                transf_mat(0,0) = 1.; transf_mat(0,1) = 0.;
                transf_mat(1,0) = 0.; transf_mat(1,1) = 1.;
                transf_mats.push_back(transf_mat);
        }
}

void    build_loc_coords()
{
        //integrate linearized system with Euler method
        double                  mu = 1.;
        std::vector<t_mat>      X_mats;
        t_mat                   X = t_mat_ident_tml<double,2,2>();
        for (int i = 0;i < cycle_pnts.size();++i) {
                t_vec           p(cycle_pnts[i]);
                t_mat           A;
                double          dt = (i > 0?cycle_times[i]-cycle_times[i-1]:cycle_times[i+1] - cycle_times[i]);
                A(0,0) = 0.;                    A(0,1) = 1.;
                A(1,0) = -1.-2.*mu*p[0]*p[1];   A(1,1) = (1.-p[0]*p[0])*mu;
                X += A*X*dt;
                
                //printf("dt = %e\n", dt);
                
                X_mats.push_back(X);
        }

        double                  T = cycle_times.back() - cycle_times.front();
        t_mat                   X_T = X_mats.back();
        t_mat                   B;      //X_T = B*X_T_diag*B^(-1)
        //X_T right eigs
        double                  lambda[2];
        lambda[0] = (X_T(0,0)+X_T(1,1) - sqrt(pow(X_T(0,0)+X_T(1,1),2) + 4.*(X_T(0,1)*X_T(1,0)-X_T(0,0)*X_T(1,1))))/2.;
        lambda[1] = (X_T(0,0)+X_T(1,1) + sqrt(pow(X_T(0,0)+X_T(1,1),2) + 4.*(X_T(0,1)*X_T(1,0)-X_T(0,0)*X_T(1,1))))/2.;
        B(0,0) = 1.; B(1,0) = (lambda[0]-X_T(0,0))/X_T(0,1);
        B(0,1) = 1.; B(1,1) = (lambda[1]-X_T(0,0))/X_T(0,1);

        printf("lambda[0] = %e; lambda[1] = %e; \n", lambda[0], lambda[1]);

        /*t_mat test00 = X_T-t_mat_ident_tml<double,2,2>()*lambda[0],
                test01 = X_T-t_mat_ident_tml<double,2,2>()*lambda[1];

        printf("%e %e\n", det22(test00), det22(test01));*/

        /*t_mat test1;
        test1(0,0) = lambda[0]; test1(0,1) = 0.;
        test1(1,0) =        0.; test1(1,1) = lambda[1];
        t_mat   test2 = X_T*B - B*test1;
        printf("%e %e\n", test2(0,0), test2(0,1));
        printf("%e %e\n", test2(1,0), test2(1,1));*/

        double                  omega[2] = {log(lambda[0])/T, log(lambda[1])/T};

        for (int i = 0;i < cycle_pnts.size();++i) {
                double          t = cycle_times[i] - cycle_times[0];
                t_mat           exp_omt_inv;
                exp_omt_inv(0,0) = exp(-omega[0]*t); exp_omt_inv(0,1) = 0.;
                exp_omt_inv(1,0) =               0.; exp_omt_inv(1,1) = exp(-omega[1]*t);
                t_mat           transf_mat = t_mat(X_mats[i]*B)*exp_omt_inv;
                transf_mats.push_back(transf_mat);
        }
}

//void  integrate(FILE *f, double _x, double _y, double r, double T, int &views_cnt)
void    integrate(FILE *f, t_TM_vec &x, double T, int &views_cnt)
{
        //t_TM_vec              x, x_prev;
        t_TM_vec                x_prev;
        t_TM_ode                ode;
        t_picard_iteration      picard_iteration;

        ode.mu = 1.;

        print_res(f, x, views_cnt, 0);
        views_cnt++;

        t_interval              t(0.);
        double                  dt = 0.0005;
        //double                        dt = 0.0001;

        for (int i = 0;i < ceil(T/dt);++i) {
                x_prev = x;

                picard_iteration.doStep(ode, x, t, dt);
                for (int j = 0;j < 2;++j)
                        x[j] = x[j].part_eval(0,1.);

                //if (i%100 == 0) printf("t = %f\n", t.lower());
        }

        print_res(f, x, views_cnt, 1);
        views_cnt++;
}

int main()
{
        int     dn = 1;
        //int   dn = 10;
        double  r = 0.01;       //this is basic for good parallelograms picard test
        //double  r = 0.03;       //for dummy qudrangles
        //double        r = 0.2;

        int     views_cnt = 0;
        
        t_TM::init_exponents_table();

        read_cycle_pnts();
        //build_loc_coords_stupid();
        build_loc_coords();

        FILE    *f = fopen("van_der_pol_neib.pos", "w");

        for (int i = 0;i < cycle_pnts.size();++i) {
                if (i%dn == 0) {
                        t_TM_vec        x;
                        x[0].polynomial()[0] = cycle_pnts[i][0];
                        x[1].polynomial()[0] = cycle_pnts[i][1];
                        int     exp1[3] = {0,1,0};
                        int     exp2[3] = {0,0,1};

                        x[0].polynomial()[t_TM::polynomial_type::term_index(exp1)] = transf_mats[i](0,0)*r*0.1;
                        x[0].polynomial()[t_TM::polynomial_type::term_index(exp2)] = transf_mats[i](0,1)*r*0.1;
                        x[1].polynomial()[t_TM::polynomial_type::term_index(exp1)] = transf_mats[i](1,0)*r*0.1;
                        x[1].polynomial()[t_TM::polynomial_type::term_index(exp2)] = transf_mats[i](1,1)*r*0.1;

                        /*x[0].polynomial()[t_TM::polynomial_type::term_index(exp1)] = 0.001;
                        x[0].polynomial()[t_TM::polynomial_type::term_index(exp2)] = 0.;
                        x[1].polynomial()[t_TM::polynomial_type::term_index(exp1)] = 0.;
                        x[1].polynomial()[t_TM::polynomial_type::term_index(exp2)] = 0.001;*/

                        //print_res(f, x, 0, 0);
                        //views_cnt++;

                        //integrate(f, x_, y_, r, 0.01*dn*30, views_cnt);
                        //for dt = 0.01
                        //integrate(f, x, 0.01*dn*10, views_cnt);
                        integrate(f, x, 0.001*dn*100, views_cnt);

                        printf("point %d final reminder:\n", i);
                        for (int j = 0;j < 2;++j) {
                                printf("j = %d: [%e,%e]\n", j, x[j].remainder().lower(), x[j].remainder().upper());
                        }
                }
                //if (i > 10)  break;
        }

        for (int i = 0;i < views_cnt;++i) {
                fprintf(f, "View[%d].ShowScale = 0;\n", i);
        }

        fclose(f);
        
        t_TM::free_exponents_table();

        return 0;
}
