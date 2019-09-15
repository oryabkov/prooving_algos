
#include "TModel.h"

void test1()
{
	printf("test1:\n");

        TModel<2,2>	m1(1.),m2(1.),m3(-1e-20), mm;

	mm = (m1 + m3) - m2;

	printf("mm a0 = %e\n", mm.polynomial()[0]);
	printf("mm R = [%e;%e]\n", mm.remainder().lower(), mm.remainder().upper());
}

void test2()
{
	printf("test2:\n");
	
	TModel<1,4>	X(0.);
	X.polynomial()[1] = 1.;
	//X is 'x' polynom
        TModel<1,4>	F1 = TModel<1,4>(1.)+X*2.+X*X*3.+X*X*X*4.,
			F2 = TModel<1,4>(5.)+X*6.+X*X*7.+X*X*X*8.,
			FF = F1*F2;

        for (int j = 0;j < 5;++j) {
		printf("F1 a%d = %e\n", j, F1.polynomial()[j]);
	}
	printf("F1 R = [%e;%e]\n", F1.remainder().lower(), F1.remainder().upper());
	
	for (int j = 0;j < 5;++j) {
		printf("F2 a%d = %e\n", j, F2.polynomial()[j]);
	}
	printf("F2 R = [%e;%e]\n", F2.remainder().lower(), F2.remainder().upper());
	
	for (int j = 0;j < 5;++j) {
		printf("FF a%d = %e\n", j, FF.polynomial()[j]);
	}
	printf("FF R = [%e;%e]\n", FF.remainder().lower(), FF.remainder().upper());

	/*mm = (m1 + m3) - m2;

	mmm = m1*m3;

	//printf("len = %d\n", TModel<2,2>::polynomial_type::size);
	//printf("a0 = %f\n", m1.polynomial()[0]);

	printf("mm a0 = %e\n", mm.polynomial()[0]);
	printf("mm R = [%e;%e]\n", mm.remainder().lower(), mm.remainder().upper());

	printf("mmm a0 = %e\n", mmm.polynomial()[0]);*/
}

int main()
{
	//polynomial<double,2,2>	p1(10.), p2(1.);
	//printf("a0 = %f\n", p1[0]);

        test1();
        test2();

	return 0;
}