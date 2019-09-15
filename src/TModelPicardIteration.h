#ifndef __TMODEL_PICARD_ITERATION_H__
#define __TMODEL_PICARD_ITERATION_H__

template<class TM_ODE,class TM_VECTOR,class TM,class TINTERVAL,class T>
class TModelPicardIteration
{
public:
	//typedef
	void	doStep(const TM_ODE &ode, TM_VECTOR &y, TINTERVAL &t, const T &dt)const
	{
		const int tord = TM::taylor_order;
		const int ode_dim = TM_ODE::dim;
		//iteration to get polynomial part (remainder is meaningless here)
		TM_VECTOR	res(y), tmp;
		for (int i = 0;i < tord+2;++i) {
			//TODO t is not right
			ode(res, TM(0.), tmp);
			for (int j = 0;j < ode_dim;++j) {
                        	/*tmp[j] = tmp[j].part_integr(0,-1.);
                        	//TODO make mul more correct!
                        	tmp[j] = (tmp[j]*dt)*(0.5);
				res[j] = y[j].part_eval(0,-1.) + tmp[j];*/

				//TODO check!!
				tmp[j] = tmp[j].part_integr(0,0.);
                        	//TODO make mul more correct!
                        	tmp[j] = (tmp[j]*dt);
				res[j] = y[j].part_eval(0,0.) + tmp[j];
			}
			for (int j = 0;j < ode_dim;++j) {
				//printf("TModel %d:\n", j);
				//TModel_print(res[j]);
			}
		}
		//res[0].remainder() = TINTERVAL(-1.,1.);
		//res[1].remainder() = TINTERVAL(-1.,1.);
		res[0].remainder() = TINTERVAL(-0.5,0.5);
		res[1].remainder() = TINTERVAL(-0.5,0.5);


		//iterations to get remainder
		//printf("remai")
		//for (int i = 0;i < tord+2;++i) {
		//do {
                TINTERVAL	prev_rem[ode_dim];
		while (1) {
                        for (int j = 0;j < ode_dim;++j) {
				prev_rem[j] = res[j].remainder();
			}
			//TODO t is not right
			ode(res, TM(0.), tmp);
			for (int j = 0;j < ode_dim;++j) {
                        	/*tmp[j] = tmp[j].part_integr(0,-1.);
                        	//TODO make mul more correct!
                        	tmp[j] = (tmp[j]*dt)*(0.5);
				res[j] = y[j].part_eval(0,-1.) + tmp[j];*/

				//TODO check!!
				tmp[j] = tmp[j].part_integr(0,0.);
                        	//TODO make mul more correct!
                        	tmp[j] = (tmp[j]*dt);
				res[j] = y[j].part_eval(0,0.) + tmp[j];
			}
			for (int j = 0;j < ode_dim;++j) {
				//printf("TModel %d:\n", j);
				//TModel_print(res[j]);
				/*T	sc[ode_dim+1];
				sc[0] = 0.1;
				sc[1] = 0.05;
				sc[2] = 0.05;*/
				//TModel_print_sc<ode_dim+1,tord,T>(res[j],sc);
			}
			bool failed = false;
                        for (int j = 0;j < ode_dim;++j) {
				if ((prev_rem[j].lower() > res[j].remainder().lower())||(prev_rem[j].upper() < res[j].remainder().upper())) failed = true;
			}
			if (failed) {
				for (int j = 0;j < ode_dim;++j) {
					res[j].remainder() = prev_rem[j]*2.;
				}
			} else {
				bool	stop = true;
				T	eps = 1e-2;
                                for (int j = 0;j < ode_dim;++j) {
					if ((std::abs(prev_rem[j].lower() - res[j].remainder().lower()) > width(res[j].remainder())*eps)||
					   (std::abs(prev_rem[j].upper() - res[j].remainder().upper()) > width(res[j].remainder())*eps)) stop = false;
					//printf("%e %e %e\n", std::abs(prev_rem[j].lower() - res[j].remainder().lower()), std::abs(prev_rem[j].upper() - res[j].remainder().upper()), width(res[j].remainder()*5e-2));
				}
				if (stop) {
					//printf("change is too small; stop remainder iterations;\n");
					break;
				}
			}
		}
		y = res;
		t += dt;
	}
};

#endif
