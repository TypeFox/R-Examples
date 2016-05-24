EBlassoNEG.Binomial <-function(BASIS,Target,a_gamma,b_gamma,Epis = FALSE,verbose = 0,group = FALSE){
	N 				= nrow(BASIS);
	K 				= ncol(BASIS);
	if (verbose>0)  cat("EBLASSO Logistic Model, NEG prior, Epis: ",Epis,"\n");
	if(Epis){
		N_effect 		= (K+1)*K/2;
		#N_effect 		= 2*K;
		Beta 			= rep(0,N_effect *4);

		#dyn.load("fEBBinaryFull.dll")

		output<-.C("fEBBinaryFull",
			BASIS 		= as.double(BASIS),
			Target 		= as.double(Target),
			a_gamma 	= as.double(a_gamma),
			b_gamma 	= as.double(b_gamma),
			logLikelihood = as.double(0),
			Beta 		= as.double(Beta),
			WaldScore 	= as.double(0),
			Intercept 	= as.double(rep(0,2)),
			N 			= as.integer(N),
			K 			= as.integer(K),
			verbose 	= as.integer(verbose),
			bMax 		= as.integer(N_effect),
			group 		= as.integer(group),
			PACKAGE="EBglmnet");
		#dyn.unload("fEBBinaryFull.dll")
	} else {
		N_effect 		= K;
		Beta 			= rep(0,N_effect *4);

		#dyn.load("fEBBinaryMainEff.dll")

		output<-.C("fEBBinaryMainEff",
			BASIS 	= as.double(BASIS),
			Target 	= as.double(Target),
			a_gamma 	= as.double(a_gamma),
			b_gamma 	= as.double(b_gamma),
			logLikelihood = as.double(0),
			Beta 		= as.double(Beta),
			WaldScore 	= as.double(0),
			Intercept 	= as.double(rep(0,2)),
			N 		= as.integer(N),
			K 		= as.integer(K),
			verbose 	= as.integer(verbose),			
			PACKAGE="EBglmnet");
		#dyn.unload("fEBBinaryMainEff.dll")
	}

	result 			= matrix(output$Beta,N_effect,4);

	ToKeep 			= which(result[,3]!=0);
	if(length(ToKeep)==0) { Blup = matrix(0,1,4)
	}else
	{
		nEff 	= length(ToKeep);
		#Blup 		= matrix(result[ToKeep,],nEff,4);
		Blup 		= result[ToKeep,,drop=FALSE];
	}
	if(Epis){
		blupMain 		= Blup[Blup[,1] ==Blup[,2],,drop = FALSE];
		#
		blupEpis 		= Blup[Blup[,1] !=Blup[,2],,drop = FALSE];
	
		order1 			= order(blupMain[,1]);
		order2 			= order(blupEpis[,1]);
		Blup 			= rbind(blupMain[order1,],blupEpis[order2,]);	
	}
	#t- test:
	t 				= abs(Blup[,3])/(sqrt(Blup[,4])+ 1e-20);
pvalue 			= 2*(1- pt(t,df=(N-1)));
	Blup 			= cbind(Blup,t,pvalue); 			#M x 6
	colnames(Blup) = c("locus1","locus2","beta","posterior variance","t-value","p-value");	
	#col1: index1
	#col2: index2
	#col3: beta
	#col4: variance
	#col5: t-value
	#col6: p-value
	hyperparameters = c(a_gamma, b_gamma);
	names(hyperparameters) = c("a", "b");
	fEBresult 			<- list(Blup,output$logLikelihood,output$WaldScore,output$Intercept[1],hyperparameters);
	rm(list= "output")	
	names(fEBresult)		<-c("fit","logLikelihood","WaldScore","Intercept","hyperparameters")
	return(fEBresult)
	
}
