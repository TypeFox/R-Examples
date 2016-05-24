elasticNetSML <-
function(Y,X,Missing,B,Verbose = 0){
	M = nrow(Y);
	N = ncol(Y);
	cat("\telastic net SML version_1;",M, "Genes, ", N , "samples; Verbose: ", Verbose, "\n\n")
	f = matrix(1,M,1);
	stat = rep(0,6);

	#dyn.load("elasticSMLv1.dll")
	tStart 	= proc.time();
	output<-.C("mainSML_adaEN",
				Y 	= as.double(Y),
				X 	= as.double(X),
				M  	= as.integer(M),
				N 		= as.integer(N),			
				Missing 	= as.integer(Missing),
				B 	= as.double(B),
				f = as.double(f),
				stat = as.double(stat),
				verbose = as.integer(Verbose),
				package = "sparseSEM"); 

	tEnd = proc.time();
	simTime = tEnd - tStart;
	#dyn.unload("elasticSMLv1.dll")
	cat("\t computation time:", simTime[1], "sec\n");

	Bout = matrix(output$B,nrow= M, ncol = M, byrow = F);
	fout = matrix(output$f,nrow= M, ncol = 1, byrow = F);
	stat = matrix(output$stat,nrow = 6,ncol = 1, byrow = F);
	stat

	SMLresult 			<- list(Bout,fout,stat,simTime[1]);
	names(SMLresult)	<-c("weight","F","statistics","simTime")
	return(SMLresult)
}
