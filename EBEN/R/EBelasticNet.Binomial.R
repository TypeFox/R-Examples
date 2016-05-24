EBelasticNet.Binomial <-
function(BASIS,Target,lambda,alpha,Epis = "no",verbose = 0){
	N 				= nrow(BASIS);
	K 				= ncol(BASIS);
	if (verbose>0) 	cat("EBEN Logistic Model, NE prior,Epis: ",Epis,"\n");
	if(Epis == "yes"){
		N_effect 		= 2*K;
		Beta 			= rep(0,N_effect *4);

		#dyn.load("fEBBinaryNeFull.dll")

		output<-.C("ElasticNetBinaryNEfull",
			BASIS 	= as.double(BASIS),
			Target 	= as.double(Target),
			lambda 	= as.double(lambda),
			alpha 	= as.double(alpha),
			logLikelihood = as.double(0),
			Beta 		= as.double(Beta),
			WaldScore 	= as.double(0),
			Intercept 	= as.double(rep(0,2)),
			N 		= as.integer(N),
			K 		= as.integer(K),
			verbose = as.integer(verbose),
			bMax 		= as.integer(N_effect),
			PACKAGE="EBEN");
		#dyn.unload("fEBBinaryNeFull.dll")
	} else {
		N_effect 		= K;
		Beta 			= rep(0,N_effect *4);

		#dyn.load("fEBBinaryNEmainEff.dll")
#dyn.load("ElasticNetBinaryNEmainEff.dll")
		output<-.C("ElasticNetBinaryNEmainEff",
			BASIS 	= as.double(BASIS),
			Target 	= as.double(Target),
			lamda 	= as.double(lambda),
			alpha 	= as.double(alpha),
			logLikelihood = as.double(0),
			Beta 		= as.double(Beta),
			WaldScore 	= as.double(0),
			Intercept 	= as.double(rep(0,2)),
			N 		= as.integer(N),
			K 		= as.integer(K),
			verbose = as.integer(verbose),
			bMax 		= as.integer(N_effect),
			PACKAGE="EBEN");
		#dyn.unload("fEBBinaryNEmainEff.dll")
	#dyn.unload("ElasticNetBinaryNEmainEff.dll")
	}
	result 			= matrix(output$Beta,N_effect,4);

	ToKeep 			= which(result[,3]!=0);
	if(length(ToKeep)==0) {  Blup = matrix(0,1,4)
	}else	Blup 	 = result[ToKeep,,drop= FALSE];
	if(Epis == "yes"){
		blupMain 		= Blup[Blup[,1] ==Blup[,2],];
		nMain 			= length(blupMain)/4;
		blupMain 		= matrix(blupMain,nMain,4);
		#
		blupEpis 		= Blup[Blup[,1] !=Blup[,2],];
		nEpis 			= length(blupEpis)/4;
		blupEpis 		= matrix(blupEpis,nEpis,4);
		
		order1 			= order(blupMain[,1]);
		order2 			= order(blupEpis[,1]);
		Blup 			= rbind(blupMain[order1,],blupEpis[order2,]);	
	}
	
	#t-test:
	t 				= abs(Blup[,3])/(sqrt(Blup[,4])+ 1e-20);
	pvalue 			= 2*(1- pt(t,df=(N-1)));
	Blup 			= cbind(Blup,t,pvalue); 			#M x 6
	#col1: index1
	#col2: index2
	#col3: beta
	#col4: variance
	#col5: t-value
	#col6: p-value
	
	fEBresult 			<- list(Blup,output$logLikelihood,output$WaldScore,output$Intercept,lambda,alpha);
	rm(list= "output")	
	names(fEBresult)		<-c("weight","logLikelihood","WaldScore","Intercept","lambda","alpha")
	return(fEBresult)
	
}
