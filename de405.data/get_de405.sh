#!/bin/sh
case $1 in
	"eph")
		wget -c "ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de405/header.405"
		num=1600
		while [ $num -lt 2201 ]
		do
			wget -c "ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de405/ascp$num.405"
			num=$(($num+20))
		done
		;;
	"test")
		wget -c "ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de405/testpo.405"
		wget -c "ftp://ssd.jpl.nasa.gov/pub/eph/planets/fortran/asc2eph.f"
		wget -c "ftp://ssd.jpl.nasa.gov/pub/eph/planets/fortran/testeph.f"
		;;
	*)
		echo "Usage: $0 <data-type>"
		echo "	<data-type> can be:"
		echo "	eph	Get DE405 ephemeris data in ASCII format"
		echo "	test	Get test results file from NASA and Fortran programs to check downloaded ephemeris"
		;;
esac
		
