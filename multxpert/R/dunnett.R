# PDunnett is a secondary function which computes the cumulative 
# distribution function of the Dunnett distribution in one-sided 
# hypothesis testing problems with a balanced one-way layout and 
# equally weighted null hypotheses
#library(mvtnorm)
pdunnett<-function(x,df,m)
# X, Argument
# DF, Number of degrees of freedom
# M, Number of comparisons
{
	# Correlation matrix
	corr<-matrix(0.5,m,m)
	for (i in 1:m) corr[i,i]<-1
	p<-pmvt(lower=rep(-Inf,m), upper=rep(x,m), delta=rep(0,m), df=df, corr=corr, algorithm=GenzBretz(maxpts=25000, abseps=0.00001, releps=0))[1]    
	return(p)
}
# End of pdunnett

# QDunnett is a secondary function which computes a quantile of the 
# Dunnett distribution in one-sided hypothesis testing problems 
# with a balanced one-way layout and equally weighted null hypotheses
#library(mvtnorm)
qdunnett<-function(x,df,m)
# X, Argument
# DF, Number of degrees of freedom
# M, Number of comparisons
{
	# Correlation matrix
	corr<-matrix(0.5,m,m)
	for (i in 1:m) corr[i,i]<-1
	temp<-qmvt(x,interval=c(0,4),tail="lower.tail",df=df, delta=rep(0,m),corr=corr, algorithm=GenzBretz(maxpts=25000, abseps=0.00001, releps=0))[1]
	return(temp$quantile)
}
# End of qdunnett

