SRC=src
INCLUDE_DIRS = -I$(SRC) -I$(SRC)/common
van_der_pol_picard_test.bin:
	g++ $(INCLUDE_DIRS) -fpermissive -O3 $(SRC)/van_der_pol_picard_test.cpp -o van_der_pol_picard_test.bin
van_der_pol_picard_neib_test.bin:
	g++ $(INCLUDE_DIRS) -fpermissive -O3 $(SRC)/van_der_pol_picard_neib_test.cpp -o van_der_pol_picard_neib_test.bin