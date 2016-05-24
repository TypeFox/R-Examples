#' @export
WilksFormula <- function(alpha=0.95,beta=0.95,bilateral=FALSE,order=1){

#--------------------------------------------------------------
# Gives the minimal size of a i.i.d. sample in order to apply |
# the Wilks formula (quantile estimation with a given         |
# confidence level)                                           |
#                                                             |
#-------------------------------------------------------------|
#Input arguments                                              |
#              aplha: order of the quantile (between 0 and 1) |
# 	                  Default: 0.95                           |
#              beta:  value of the level of the confidence    |
#                       interval (between 0 and 1).           |
# 	                  Default: 0.95                           |
#              bilateral: TRUE for bilateral quantile         |
#                           FALSE for unilateral quantile     |
#                     Default: FALSE                          |
#              order: order of the Wilks formula              |
#                     Default: 1                              |
#-------------------------------------------------------------|
#Output: the value of the minimal N                           |
#-------------------------------------------------------------|
#      Authors: Paul Lemaitre, Bertrand Iooss                 |
#-------------------------------------------------------------|

# treatment of uni or bilateral cases: double the order if bilateral
if (bilateral=="TRUE")
 order=2*order

# we determine the Wilks formula
FWilks <- function(N,alpha,beta,order){
	coef=matrix(NA,nrow=order)
	coef[1]=1
	if (order>=2){
		for (i in 2:order){
		coef[i]=coef[i-1]*(N-i+2)/factorial(i-1)*factorial(i-2)}}
	else {}

FWilks=1;
	for (i in 1:order){
	FWilks=FWilks-coef[i]*alpha^(N-i+1)*(1-alpha)^(i-1)}
FWilks=FWilks-beta
return (FWilks)
}

nmax=1000

if (order!=1) {
	while (FWilks(nmax,alpha,beta,order)<0) {nmax=nmax+100}
	nmin=nmax-1000+1

	# searching the n such as FWilks=0 by dichotomy
	sortie=10

	while (sortie>1) {
		n=ceiling((nmax+nmin)/2)
    		FWilksn=FWilks(n,alpha,beta,order)

    		if (FWilksn<0) {nmin=n}
    		if (FWilksn>0) {nmax=n}
    		if (FWilksn==0) {sortie=0}
    			
    		FWilksn_1=FWilks(n-1,alpha,beta,order)
    		FWilksnp1=FWilks(n+1,alpha,beta,order)
   
    		if (FWilksn_1<0 && FWilksn>0) {
        		if (abs(FWilksn)>abs(FWilksn_1)) {n=n-1}
         	      	sortie=0 
						}#end if

    		if (FWilksn<0 && FWilksnp1>0)    {
        		if (abs(FWilksn)>abs(FWilksnp1)) {n=n+1}
           	    	sortie=0
						}#end if
    
	}# end while

	} #end if

if (order==1) {n=ceiling(log(1-beta)/log(alpha))}

return (N=n)
}
