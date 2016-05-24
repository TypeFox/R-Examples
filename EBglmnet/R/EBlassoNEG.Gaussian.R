EBlassoNEG.Gaussian <-
function(BASIS,Target,a_gamma,b_gamma,Epis = FALSE,verbose = 0,group = FALSE){
	N 					= nrow(BASIS);
	K 					= ncol(BASIS);
	if (verbose>0) cat("EBLASSO Gaussian Model, NEG prior, N: ",N,",K: ",K,", Epis: ",Epis,"\n");
	if(Epis){
		N_effect 		= (K+1)*K/2;
		#N_effect 		= 2*K;
		Beta 			= rep(0,N_effect *4);

		#dyn.load("fEBLinearFullFloat.so")

		output<-.C("fEBLinearEpisEff",
			BASIS 		= as.double(BASIS),
			Target 		= as.double(Target),
			a_gamma 	= as.double(a_gamma),
			b_gamma 	= as.double(b_gamma),
			Beta 		= as.double(Beta),
			WaldScore 	= as.double(0),
			Intercept 	= as.double(0),
			N 			= as.integer(N),
			K 			= as.integer(K),
			ver 		= as.integer(verbose),
			bMax 		= as.integer(N_effect),
			residual 	= as.double(0),			
			group 		= as.integer(group),
			PACKAGE 	="EBglmnet");
		#dyn.unload("fEBLinearFullFloat.so")
	}else {
		N_effect 		= K;
		Beta 			= rep(0,N_effect *4);
		#dyn.load("fEBLinearMainEff.so")

		output<-.C("fEBLinearMainEff",
			BASIS 		= as.double(BASIS),
			Target 		= as.double(Target),
			a_gamma 	= as.double(a_gamma),
			b_gamma 	= as.double(b_gamma),
			Beta 		= as.double(Beta),
			WaldScore 	= as.double(0),
			Intercept 	= as.double(0),
			N 			= as.integer(N),
			K 			= as.integer(K),
			ver 		= as.integer(verbose),
			residual 	= as.double(0),
			PACKAGE		="EBglmnet");
#		dyn.unload("fEBLinearMainEff.so")
	}	
	
	result 				= matrix(output$Beta,N_effect,4);
	ToKeep 				= which(result[,3]!=0);
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
	fEBresult 			<- list(Blup,output$WaldScore,output$Intercept,output$residual,hyperparameters);
	rm(list= "output")	
	names(fEBresult)	<-c("fit","WaldScore","Intercept","residual variance","hyperparameters")
	return(fEBresult)
	
}
