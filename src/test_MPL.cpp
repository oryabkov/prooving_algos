
#include <boost/mpl/bool.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/at.hpp>

typedef boost::mpl::vector<boost::mpl::bool_<true>,boost::mpl::bool_<false>,boost::mpl::bool_<true> > t_CT_vec;

int main()
{
	printf("%d\n",boost::mpl::size<t_CT_vec>::type::value);
	//for (int i = 0;i < boost::mpl::size<t_CT_vec>::type::value;++i) {
	printf("%d ", boost::mpl::at<t_CT_vec,boost::mpl::int_<0> >::type::value);
	printf("%d ", boost::mpl::at<t_CT_vec,boost::mpl::int_<1> >::type::value);
	printf("%d ", boost::mpl::at<t_CT_vec,boost::mpl::int_<2> >::type::value);
	//}
	return 0;
}
