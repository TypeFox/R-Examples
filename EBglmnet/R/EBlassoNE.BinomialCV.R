EBlassoNE.BinomialCV <- function(BASIS,Target,nFolds,foldId,Epis=FALSE, verbose = 0)
{
nStep= 19
	cat("EBLASSO Logistic Model, NE prior,Epis: ",Epis, ";", nFolds, "fold cross-validation\n");
	N 					= nrow(BASIS);
	K 					= ncol(BASIS);
	#set.seed(proc.time())
	if (missing(foldId)) 
	{
		if(N%%nFolds!=0){
			foldId 			= sample(c(rep(1:nFolds,floor(N/nFolds)),1:(N%%nFolds)),N);
		}else{
			foldId 			= sample(rep(1:nFolds,floor(N/nFolds)),N);
		}
	}
	lambda_Max = lambdaMax(BASIS,Target,Epis);
	lambda_Min 			= log(0.0001*lambda_Max);
	step 				= (log(lambda_Max) - lambda_Min)/nStep;
	Lambda 				= exp(seq(from = log(lambda_Max),to=lambda_Min,by= -step))
	N_step 				= length(Lambda);

	step 				= 1;
	alpha 				= 1;
	nAlpha 				= 1;
		
	Likelihood 			= mat.or.vec((N_step*nAlpha),4);
	#logL 				= mat.or.vec(nFolds,1);
	nLogL = rep(0,4);
	pr = "lasso"; #1LassoNEG; 2: lasso; 3EN
	model = "binomial";#0linear; 1 binomial
	group = 0;
	
	for (i_s in 1:N_step){			
			lambda 		= Lambda[i_s];

	if(verbose >=0) cat("\tTesting step", step, "\t\tlambda: ",lambda,"\t")
				
			hyperpara = c(alpha, lambda);
			logL = CVonePair(BASIS,Target,nFolds, foldId,hyperpara,Epis,pr,model,verbose,group);	
			
			logL[3] = -logL[3]; #C produces negative logL;
			Likelihood[step,] = logL;
			if(verbose >=0) cat("log Likelihood",logL[3],"\n");
			step 			= step + 1;
	}
		colnames(Likelihood) = c("alpha","lambda","logLikelihood","standard error");	

	index 				= which.max(Likelihood[,3]);
	Res.lambda			= Likelihood[index,2];
	Res.alpha 			= Likelihood[index,1];
	result 				<- list(Likelihood,Res.lambda);
	names(result)		<-c("CrossValidation","optimal hyperparameter");
	return(result);
}
