EBlassoNEG.GaussianCV <-function(BASIS,Target,nFolds,foldId,Epis=FALSE, verbose = 0, group = FALSE){
	cat("EBLASSO Linear Model, NEG prior,Epis: ",Epis, ";", nFolds, "fold cross-validation\n");
	N 				= nrow(BASIS);
	#set.seed(proc.time())
	if (missing(foldId)) 
	{
		if(N%%nFolds!=0){
			foldId 		= sample(c(rep(1:nFolds,floor(N/nFolds)),1:(N%%nFolds)),N);
		}else{
			foldId 		= sample(rep(1:nFolds,floor(N/nFolds)),N);
		}
	}
	a_r1 			= c(0.01, 0.05, 0.1,0.5,1);
	b_r1 			= a_r1;
	N_step1 		= length(a_r1);
	a_r2 			= c(1, 0.5, 0.1, 0.05, 0.01, -0.01,-0.1,-0.2,-0.3,-0.4,-0.5, -0.6, -0.7, -0.8, -0.9);
	b_r2 			= c(0.01, 0.05, 0.1,0.5,1);
	N_step2 		= length(a_r2) -1;
	N_step3 		= length(b_r2);	
	N_step  		= N_step1 + N_step2 + N_step3;	
	#Likelihood 		= mat.or.vec(N_step,4);
	#logL 			= mat.or.vec(nFolds,1);
	MeanSqErr 					= mat.or.vec(N_step,4);
	SSE 						= matrix(rep(0,nFolds),nFolds,1);
	stp 			 = 1;
	
	nLogL = rep(0,4);
	pr = "lasso"; #1LassoNEG; 2: lasso; 3EN
	model = "gaussian";#0linear; 1 binomial
	
	#------------------------------------------ step one ----------------------------------
	for (i_s1 in 1:N_step1){		
		a_gamma 			= a_r1[i_s1];
		b_gamma 			= b_r1[i_s1];
	if(verbose >=0) cat("Testing step", stp, "\t\ta: ",a_gamma, "b: ", b_gamma,"\t")
			
			hyperpara = c(a_gamma, b_gamma);
			logL = CVonePair(BASIS,Target,nFolds, foldId,hyperpara,Epis,pr,model,verbose,group);
			
		if(verbose >=0) cat("sum squre error",logL[3],"\n");
		MeanSqErr[stp,] 		= logL;
		stp 					= stp + 1;
	}	
index		 				= which.min(MeanSqErr[1:N_step1,3]);
previousL 					= MeanSqErr[index,3] + MeanSqErr[index,4];
previousMin 				= MeanSqErr[index,3];		
b_gamma 					= b_r1[index];
index 						= which(a_r2>=b_gamma);
a_rS2 						= a_r2[-index]; 			# starts at smaller a
N_step2 					= length(a_rS2)	
	#------------------------------------------ step two ----------------------------------	
	for(i_s2 in 1:N_step2){
		a_gamma 			= a_rS2[i_s2];
	if(verbose >=0) cat("Testing step", stp, "\t\ta: ",a_gamma, "b: ", b_gamma,"\t")
			hyperpara = c(a_gamma, b_gamma);
			logL = CVonePair(BASIS,Target,nFolds, foldId,hyperpara,Epis,pr,model,verbose,group);
			
		if(verbose >=0) cat("sum squre error",logL[3],"\n");
		MeanSqErr[stp,] 		= logL;
		currentL				= MeanSqErr[stp,3];
		stp 					= stp + 1;	
		# break out of 2nd step
	if((currentL - previousL)>0 && a_gamma <0){break;}
		if(currentL < previousMin)
		{
			preStp 				= stp -1;
			previousL 			= MeanSqErr[preStp,3] + MeanSqErr[preStp,4];	
			previousMin  		= MeanSqErr[preStp,3];
		}	
	}
	nStep 						= stp - 1;
	index 						= which.min(MeanSqErr[1:nStep,3]);
	a_gamma 					= MeanSqErr[index,1];
	previousL 					= MeanSqErr[index,3] + MeanSqErr[index,4];	
	previousMin 				= MeanSqErr[index,3];	
	bstep2 						= MeanSqErr[index,2];
	b_rS2 						= b_r2;
	index 						= which(b_r2==bstep2)
	Nbcut 						= length(index);	
	b_rS2 					= b_r2[-index];
	N_step3 				= N_step3 - Nbcut;

	#------------------------------------------ step three ----------------------------------
	for(i_s3 in 1:N_step3){
		b_gamma 			= b_rS2[i_s3];
	if(verbose >=0) cat("Testing step", stp, "\t\ta: ",a_gamma, "b: ", b_gamma,"\t")
					hyperpara = c(a_gamma, b_gamma);
			logL = CVonePair(BASIS,Target,nFolds, foldId,hyperpara,Epis,pr,model,verbose,group);
			
		if(verbose >=0) cat("sum squre error",logL[3],"\n");
	MeanSqErr[stp,] 		= logL;
	#currentL				= Likelihood[stp,3] + Likelihood[stp,4];
	currentL				= MeanSqErr[stp,3];
	stp 					= stp + 1;			
		# break out of 3rd step
		if((currentL - previousL)>0){break;}
		if(currentL < previousMin)
		{
			preStp 				= stp -1;
			previousL 			= MeanSqErr[preStp,3] + MeanSqErr[preStp,4];
			previousMin  		= MeanSqErr[preStp,3];			
		}		
	}
	nStep = stp - 1;
	#index 					= which.max(Likelihood[1:nStep,3]);
	index 						= which.min(MeanSqErr[1:nStep,3]);
	a_gamma 					= MeanSqErr[index,1];	
	b_gamma 					= MeanSqErr[index,2];
	opt_para 				= c(a_gamma,b_gamma);
	names(opt_para) 		= c("a_optimal","b_optimal");
	colnames(MeanSqErr) = c("a","b","Mean Square Error","standard error");	
	result 					<- list(MeanSqErr[1:nStep,],opt_para);
	names(result)			<-c("CrossValidation","optimal hyperparameter");
	return(result);
}
