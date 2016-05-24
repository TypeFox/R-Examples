EBelasticNet.GaussianCV <-
function(BASIS,Target,nFolds,foldId,Epis=FALSE, verbose = 0)
{

	nStep = 19;
	#early stop: for each alpha, if next lambda > SSEmin, then stop.
	cat("EB-Elastic Net Linear Model, Epis: ",Epis, ";", nFolds, "fold cross-validation\n");
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
	#Alpha 				= seq(from = 1, to = 0.1, by = -0.1)
	Alpha 				= seq(from = 0.9, to = 0.1, by = -0.1)
	nAlpha 				= length(Alpha);
		
	MSEcv 				= mat.or.vec((N_step*nAlpha),4);
	MSEeachAlpha		= mat.or.vec(nAlpha,4); # minimum MSE for each alpha
	MeanSqErr 			= mat.or.vec(nFolds,1);
	SSE1Alpha			= matrix(1e10,N_step,2);# temp matrix to keep MSE + std in each step

	nLogL = rep(0,4);
	pr = "elastic net"; #1LassoNEG; 2: lasso; 3EN
	model = "gaussian";#0linear; 1 binomial
	group = 0;
	
	for(i_alpha in 1:nAlpha){
		alpha 					= Alpha[i_alpha];
		SSE1Alpha				= matrix(1e10,N_step,2);# temp matrix to keep MSE + std in each step

		if(verbose >=0) cat("Testing alpha", i_alpha, "/",nAlpha,":\t\talpha: ",alpha,"\t")
		for (i_s in 1:N_step){			
			lambda 				= Lambda[i_s];
			min_index 			= which.min(SSE1Alpha[1:(i_s -1),1]);
			previousL 			= SSE1Alpha[min_index,1] + SSE1Alpha[min_index,2];
			#cat("\tTesting step", step, "\t\tlambda: ",lambda,"\t")

			hyperpara = c(alpha, lambda);
			logL = CVonePair(BASIS,Target,nFolds, foldId,hyperpara,Epis,pr,model,verbose,group);
			SSE1Alpha[i_s,] 	= logL[3:4];
			#cat("sum squre error",mean(MeanSqErr),"\n");
			MSEcv[step,]		= logL;
			currentL			= MSEcv[step,3];
			step 				= step + 1;
			# break out of 2nd for loop
			if((currentL - previousL)>0){break;}
		}
		index 					= which.min(SSE1Alpha[,1]);
		lambda 					= Lambda[index];
		if(verbose >=0) cat("lambda: ",lambda, "\t minimum square error: ",SSE1Alpha[index,1],"\n");
		MSEeachAlpha[i_alpha,] 	= c(alpha,lambda, SSE1Alpha[index,]);
	}
	colnames(MSEeachAlpha) = c("alpha","lambda","Mean Square Error","standard error");
	colnames(MSEcv) = c("alpha","lambda","Mean Square Error","standard error");
	index 						= which.min(MSEeachAlpha[,3]);
	Res.lambda					= MSEeachAlpha[index,2];
	Res.alpha 					= MSEeachAlpha[index,1];
	opt_para 				= c(Res.alpha,Res.lambda);
	result 						<- list(MSEcv,opt_para);
	names(result)				<-c("CrossValidation","optimal hyperparameter");
	return(result);
}
