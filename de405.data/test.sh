#!/bin/sh
echo "Running DE405 ephemeris test"
echo "	*** Check ephemeris data\n"
if [ ! -f ascp1620.405 ]; then
	./get_de405.sh eph
fi

echo "	*** Check test data\n"
if [ ! \( -f asc2eph.f -a -f testeph.f -a -f testpo.405 \) ]; then
	./get_de405.sh test
fi

echo "	*** Check gfortran compiler"
if [ ! `which gfortran` ]; then
	echo "	!!! gfortran not found"
	echo "	You can get GNU Fortran compiler at http://gcc.gnu.org/wiki/GFortranBinaries"
	exit 1
fi

echo "	*** Prepare dir for tests\n"
rm -rf test
mkdir -p test
cp testpo.405 test/

echo "	*** Join ASCII ephemeris data to single file\n"
cat header.405 ascp*.405 > test/full.405

echo "	*** Compile testeph.f\n"
gfortran testeph.f -o test/testeph
if [ $? -ne 0 ]; then
	echo "\n	!!! Something is wrong"
	echo "	!!! You need to fix errors in testeph.f file"
	echo "	!!! You can use testeph.f.mod as example for UNIX (run 'diff testeph.f testeph.f.mod')"
	exit 1
fi

echo "	*** Compile asc2eph.f\n"
gfortran asc2eph.f -o test/asc2eph 
if [ $? -ne 0 ]; then
	echo "\n	!!! Something is wrong"
	echo "	!!! You need to fix errors in asc2eph.f file"
	echo "	!!! You can use asc2eph.f.mod as example for UNIX (run 'diff asc2eph.f asc2eph.f.mod')"
	exit 1
fi

echo "	*** Build JPLDE binary file from ASCII ephemeris\n"
(cd test; ./asc2eph < full.405 | tee asc2eph.out)
if [ ! -f ./test/JPLEPH ]; then
	echo "\n	!!! File JPLEPH was not built"
	echo "	!!! Try to fix asc2eph.f"
	echo "	!!! You can use asc2eph.f.mod as example for UNIX (run 'diff asc2eph.f asc2eph.f.mod')"
	exit 1
fi

echo "	*** Run tests\n"
(cd test; ./testeph < testpo.405 | tee testeph.out )

echo "	*** Tests finished"
echo "	You can find asc2eph and testeph output in files ./test/asc2eph.out and ./test/testeph.out"
