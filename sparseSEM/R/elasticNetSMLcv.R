
elasticNetSMLcv <- function(Y,X,Missing,B,alpha_factors = seq(1,0.05, -0.05), lambda_factors =10^seq(-0.2,-4,-0.2) ,Verbose = 0){
	M 					= nrow(Y);
	N 					= ncol(Y);
	cat("\telastic net SML version_1;",M, "Genes, ", N , "samples; Verbose: ", Verbose, "\n\n")
	f 					= matrix(1,M,1);
	stat 				= rep(0,6);
#------------------------------------------------------R_package parameter
	nAlpha 	= length(alpha_factors);
	nLambda = length(lambda_factors);
	mseStd 	= rep(0,nLambda*2);
#------------------------------------------------------R_package parameter

	#dyn.load("elasticSMLv1.dll")
	tStart 				= proc.time();
	output<-.C("mainSML_adaENcv",
				Y 		= as.double(Y),
				X 		= as.double(X),
				M  		= as.integer(M),
				N 		= as.integer(N),			
				Missing 	= as.integer(Missing),
				B 		= as.double(B),
				f 		= as.double(f),
				stat 	= as.double(stat),
				alpha 	= as.double(alpha_factors),
				nAlpha 	= as.integer(nAlpha),
				lambda 	= as.double(lambda_factors),
				nLambda = as.integer(nLambda),
				mseStd 	= as.double(mseStd),
				verbose = as.integer(Verbose),
				package = "sparseSEM"); 

	tEnd 				= proc.time();
	simTime 			= tEnd - tStart;
	#dyn.unload("elasticSMLv1.dll")
	cat("\t computation time:", simTime[1], "sec\n");

	Bout = matrix(output$B,nrow= M, ncol = M, byrow = F);
	fout = matrix(output$f,nrow= M, ncol = 1, byrow = F);
	stat = matrix(output$stat,nrow = 6,ncol = 1, byrow = F);
#------------------------------------------------------R_package parameter
	mseStd = matrix(output$mseStd,nrow= nLambda, ncol = 2, byrow = F);
	mseStd = cbind(lambda_factors,mseStd);
	colnames(mseStd)<-c("lambda","mean Error", "Std");
#------------------------------------------------------R_package parameter


	SMLresult 			<- list(Bout,fout,stat,simTime[1],mseStd);
	names(SMLresult)	<-c("weight","F","statistics","simTime","residuals")
	return(SMLresult)
}

