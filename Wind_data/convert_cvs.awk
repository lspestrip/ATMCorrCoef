#!/usr/bin/awk -f

##
##Parser inline for data reduction. Convert the raw style data in CSV format.  
##Usage: ./convert_cvs.awk inputfilename.dat > outputfilename.dat
##Date: 6-3-2016


BEGIN{
	TIME = 0	
	print "#x,y,z,t";
}

{
	if( $1 == "M:x" ){
		print TIME" "$6/100." "$3/100." "$9/100.","$12/100.
		TIME = TIME + 1
	}else{
	}
}
END{

}



