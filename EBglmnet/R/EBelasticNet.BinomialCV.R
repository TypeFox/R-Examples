EBelasticNet.BinomialCV <- function(BASIS,Target,nFolds,foldId,Epis=FALSE,verbose = 0)
{
	nStep = 19;
	cat("EBEN Logistic Model,Epis: ",Epis, ";", nFolds, "fold cross-validation\n");
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

	lambda_Min 			= log(0.001*lambda_Max);
	step 				= (log(lambda_Max) - lambda_Min)/nStep;
	Lambda 				= exp(seq(from = log(lambda_Max),to=lambda_Min,by= -step))
	N_step 				= length(Lambda);

	step 				= 1;
	Alpha 				=  seq(from = 0.9, to = 0.1, by = -0.1)
	nAlpha 				= length(Alpha);
		
	Likelihood 			= mat.or.vec((N_step*nAlpha),4);
	#logL 				= mat.or.vec(nFolds,1);
	logL1alpha			= matrix(0,N_step,2);# temp matrix to keep MSE + std in each step

	nLogL = rep(0,4);
	pr = "elastic net"; #1LassoNEG; 2: lasso; 3EN
	model = "binomial";#0linear; 1 binomial
	group = 0;
	for(i_alpha in 1:nAlpha){
		alpha 			= Alpha[i_alpha];
		if(verbose >=0) cat("Testing alpha", i_alpha, "/",nAlpha,":\t\talpha: ",alpha,"\t")
		for (i_s in 1:N_step){
			
			lambda 		= Lambda[i_s];
			
			hyperpara = c(alpha, lambda);
			logL = CVonePair(BASIS,Target,nFolds, foldId,hyperpara,Epis,pr,model,verbose,group);
			logL[3] = -logL[3]; #C produces negative logL;
			Likelihood[step,] = logL;
			logL1alpha[i_s,] 	= logL[3:4];

			step 			= step + 1;
		}
		index = which.max(logL1alpha[,1]);
		lambda= Lambda[index];
		if(verbose >=0) cat("lambda:",lambda,"\t max log Likelihood",logL1alpha[index,1],"\n");		
	}
	colnames(Likelihood) = c("alpha","lambda","logLikelihood","standard error");
	index 				= which.max(Likelihood[,3]);
	Res.lambda			= Likelihood[index,2];
	Res.alpha 			= Likelihood[index,1];
	opt_para 				= c(Res.alpha,Res.lambda);	
	result 				<- list(Likelihood,opt_para);
	names(result)		<-c("CrossValidation","optimal hyperparameter");
	return(result);
}
