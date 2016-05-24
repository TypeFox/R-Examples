EBelasticNet.Gaussian <-
function(BASIS,Target,lambda,alpha,Epis = "no",verbose = 0 ){
	N 				= nrow(BASIS);
	K 				= ncol(BASIS);
	if (verbose>0) cat("EBEN Gaussian Model, Epis: ",Epis,"\n");
	if(Epis == "yes"){
		N_effect 		= (K+1)*K/2;
#-----------------------------------------
		blupCol 		=5;
		Beta 			= rep(0,N_effect *blupCol);
#-----------------------------------------	


#dyn.load("fEBLinearNeFull.dll")

		output<-.C("elasticNetLinearNeEpisEff",
			BASIS 	= as.double(BASIS),
			Target 	= as.double(Target),
			lamda 	= as.double(lambda),
			alpha 	= as.double(alpha),
			Beta 		= as.double(Beta),
			WaldScore 	= as.double(0),
			Intercept 	= as.double(0),
			N 		= as.integer(N),
			K 		= as.integer(K),
			verbose = as.integer(verbose),
			residual = as.double(0),
			PACKAGE="EBEN");
	
#dyn.unload("elasticNetLinearNeFull2.dll")	

	}else {
		N_effect 		= K;
		blupCol =4;
		Beta 			= rep(0,N_effect *blupCol);

#		dyn.load("fEBLinearNeMainEff.so")

		output<-.C("elasticNetLinearNeMainEff",
			BASIS 	= as.double(BASIS),
			Target 	= as.double(Target),
			lamda 	= as.double(lambda),
			alpha 	= as.double(alpha),
			Beta 		= as.double(Beta),
			WaldScore 	= as.double(0),
			Intercept 	= as.double(0),
			N 		= as.integer(N),
			K 		= as.integer(K),
			verbose = as.integer(verbose),
			residual = as.double(0),
			PACKAGE="EBEN");
#		dyn.unload("fEBLinearNeMainEff.so")

	}		
#-------------------------------------------------------------------
	if(Epis == "yes"){	
		result 			= matrix(output$Beta,N_effect,blupCol); #5th column: Used.
		ToKeep 			= which(result[,5]!=0);	
	}else
	{
		result 			= matrix(output$Beta,N_effect,blupCol);
		ToKeep 			= which(result[,3]!=0);
	}
	if(length(ToKeep)==0) { Blup = matrix(0,1,blupCol); 
	}else
	{
		nEff 	= length(ToKeep);
		Blup 		= matrix(result[ToKeep,],nEff,blupCol);#doesn't matter if blupCol=4: in R, if columns is small, only keep the first blupCol columns.
	}
	if(Epis == "yes"){
		blupMain 		= Blup[Blup[,1] ==Blup[,2],];
		nMain 			= length(blupMain)/blupCol;
		blupMain 		= matrix(blupMain,nMain,blupCol);
		#
		blupEpis 		= Blup[Blup[,1] !=Blup[,2],];
		nEpis 			= length(blupEpis)/blupCol;
		blupEpis 		= matrix(blupEpis,nEpis,blupCol);
		
		order1 			= order(blupMain[,1]);
		order2 			= order(blupEpis[,1]);
		Blup 			= rbind(blupMain[order1,],blupEpis[order2,]);	
	}
	#t-test:
	Blup 			= Blup[,1:4,drop = FALSE]; 			#will not report 5th column of Epis model;
	t 				= abs(Blup[,3])/(sqrt(Blup[,4])+ 1e-20);
	pvalue 			= 2*(1- pt(t,df=(N-1)));
	Blup 			= cbind(Blup,t,pvalue); 			#M x 6
	#col1: index1
	#col2: index2
	#col3: beta
	#col4: variance
	#col5: t-value
	#col6: p-value
	
	
	fEBresult 			<- list(Blup,output$WaldScore,output$Intercept,output$residual,lambda,alpha);
	rm(list= "output")	
	names(fEBresult)		<-c("weight","WaldScore","Intercept","residVar","lambda","alpha")
	return(fEBresult)
	
}
