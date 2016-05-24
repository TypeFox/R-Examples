lambdaMax <-function(X,y,Epis=FALSE)
{	
	N = nrow(X);
	K = ncol(X);
	output<-.C("ProjectCorr",
			N 	= as.integer(N),
			P  	= as.integer(K),
			y0 		= as.double(y),
			BASIS 	= as.double(X),
			lmax 	= as.double(0),
			epis 	= as.integer(Epis),
	PACKAGE 	="EBglmnet");
	lambda_Max = output$lmax
	return(lambda_Max);
}

