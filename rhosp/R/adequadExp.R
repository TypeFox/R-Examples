"adequadExp" <-
function(T,toplot=FALSE)
{
# test the adequation of the variable T to an exponential distribution
# the result is a list of the result of a linear regression on the 
# probability graph and the result of the Kolmogorov-Smirnov
	n<-length(T$T)
	abcisse <- log(sort(T$T));
	ordinate <- log(-log(abs(1-seq(1:n)/n)));
	reg <- lm(ordinate [1:n-1]~abcisse[1:n-1]);
	if(toplot)
	{
		plot(abcisse,ordinate );
		abline(reg);
	}
	kstest<-ks.test(T$T,"pexp");
	
	return (list(LinearRegression=reg,KolmogorvSmirnovTest=kstest))
}

