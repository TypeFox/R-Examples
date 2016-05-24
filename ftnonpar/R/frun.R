"frun"<-function(y,...,alpha=0.5,r=0,mr=0)
{
#	
# INPUTS:		
#	y:	the data
#    ... :	an optional argument which specifies
#		the approximate positions of the local extreme values. 
#		These should be consisten with the run length otherwise 
#		the result will be nonsensical. Should you wish to use
#		this option then you should first run the macro
#		without it. The item list$xb of the output list gives 
#		the acceptable limits of the local extreme values. You
#		can then specify the positions within these limits.
#   alpha:	Qunatile determining the acceptable run length.
#       r:	Acceptable run length:	 Overrides alpha if not 0
#      mr:	mr=0 minimizes the run length consistent with the
#		number of local extreme values found for the specified 
#		run length. mr=1 disables the option.
#
#		
	args<-list(y,...)
	k<-length(args)
	y<-args[[1]]
	nxl<-1
	il<-0
	if (k>1)
		{
		il<-args[[2]]
		nxl<-length(il)
		}
		
#	sample size n
	n<- length(y)
#
# fortran subroutine
#
    	tmp <- .Fortran(
		"frun", 
		as.double(y), 
		double(n), 
		double(n), 
		double(2*n),
		double(2*n),
		double(n),
		integer(n),
		integer(n),
		as.integer(n),
		as.integer(r),
		as.double(alpha),
		as.integer(mr),
		as.integer(nxl),
		integer(1),
                PACKAGE="ftnonpar"
		)
#
# OUTPUTS		
#
#	number of extremes
	nx<-tmp[[14]]
#
#	run length: may differ from specified if mr=1	 
	r<-tmp[[10]]
	
#    	lower bounds 
	l1<- tmp[[2]]
#
#       upper bound
	u1<- tmp[[3]]
#
#	bounds for location of extremes: the position of the ith
#	extreme value lies between xb[2*i-1] and xb[2*i]		 
	xb<-tmp[[7]]
	xb<-xb[1:(2*nx)]
#
#       lower bound with specified extremes: the default choices for
#	the positions of the local extreme values are the mid-points
#	of the intervals specified by xb above.			
	l2<- tmp[[4]]
	l2<-l2[1:n]
#
#	upper bound with specified extremes
 	u2<- tmp[[5]]
	u2<-u2[1:n]
#
#	function between l2 and u2 satisfying run condition
	f<- tmp[[6]]
#
#       IN GENERAL THE MEAN OF THE BOUNDS l2 AND u2 (l2+u2)/2
#	GIVES A BETTER REGRESSION FUNCTION THAN f. HOWEVER THIS
#	FUNCTION IS INFINITE AT THE TWO ENDPOINTS AND AT LOCAL EXTREME
#	VALUES. IN THESE INTERVALS IT CAN BE REPLACED BY ANY VALUES
#	WHICH DO NOT ALTER THE NUMBER OF LOCAL EXTREME VALUES. THE
#	MEDIAN OF THE y-VALUES IN THESE INTERVALS IS A REASONABLE 
#	CHOICE.
#			 				
	list(l1=l1,u1=u1,l2=l2,u2=u2,f=f,xb=xb,nx=nx,r=r)
#
}


