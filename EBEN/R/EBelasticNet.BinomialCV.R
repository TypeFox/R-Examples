EBelasticNet.BinomialCV <- function(BASIS,Target,nFolds,Epis="no",foldId = 0)
{

	cat("EB-Elastic Net Linear Model, Epis: ",Epis, ";", nFolds, "fold cross-validation\n");
	N 					= nrow(BASIS);
	K 					= ncol(BASIS);
	set.seed(1)
	if(length(foldId)!=N)
	{
		if(N%%nFolds!=0){
			foldId 			= sample(c(rep(1:nFolds,floor(N/nFolds)),1:(N%%nFolds)),N);
		}else{
			foldId 			= sample(rep(1:nFolds,floor(N/nFolds)),N);
		}
	}
	lambda_Max			= log(1.1);
	response 			= Target-mean(Target);
	response 			= response/sqrt(sum(response*response));
	for(i_b in 1:K){
		basis 			= BASIS[,i_b];
		basis 			= basis/sum(basis*basis);
		corBy 			= basis%*%response;
		if(corBy>lambda_Max) lambda_Max = corBy;
	}
	
	if(Epis == "yes"){
		Alpha 				= seq(from = 1, to = 0.05, by = -0.1)
		for(i_b in 1:(K-1)){
			for(i_bj in (i_b + 1):K){
				basis 	= BASIS[,i_b]*BASIS[,i_bj];
				basis 	= basis/sum(basis*basis);
				corBy 	= basis%*%(Target-mean(Target));
				if(corBy>lambda_Max) lambda_Max = corBy;
			}
		}		
	}
	lambda_Max 			= lambda_Max*10;
	lambda_Min 			= log(0.001*lambda_Max);
	step 				= (log(lambda_Max) - lambda_Min)/19;
	Lambda 				= exp(seq(from = log(lambda_Max),to=lambda_Min,by= -step))
	N_step 				= length(Lambda);


	step 				= 1;
	if(Epis == "no"){
		Alpha 				= seq(from = 1, to = 0.05, by = -0.05)
	}
	nAlpha 				= length(Alpha);
		
	Likelihood 			= mat.or.vec((N_step*nAlpha),4);
	MSEeachAlpha		= mat.or.vec(nAlpha,4); # minimum MSE for each alpha
	logL 				= mat.or.vec(nFolds,1);
	SSE1Alpha			= matrix(1e10,N_step,2);# temp matrix to keep MSE + std in each step
	for(i_alpha in 1:nAlpha){
		alpha 			= Alpha[i_alpha];
		cat("Testing alpha", i_alpha, "/",nAlpha,":\t\talpha: ",alpha,"\n")
		SSE1Alpha				= matrix(1e10,N_step,2);# temp matrix to keep MSE + std in each step
		
		for (i_s in 1:N_step){
			
			lambda 		= Lambda[i_s];
			min_index 			= which.min(SSE1Alpha[1:(i_s -1),1]);
			previousL 			= SSE1Alpha[min_index,1] + SSE1Alpha[min_index,2];
	cat("\tTesting step", step, "\t\tlambda: ",lambda,"\n")
			for(i in 1:nFolds){
				index  			= which(foldId!=i);
				Basis.Train 	= BASIS[index,];
				Target.Train 	= Target[index];
				index  			= which(foldId == i);
				Basis.Test  	= BASIS[index,];
				Target.Test 	= Target[index];
				SimF2fEB 		<-EBelasticNet.Binomial(Basis.Train,Target.Train,lambda,alpha,Epis,verbose=0);
				M				= length(SimF2fEB$weight)/6;
				Betas 			<- matrix(SimF2fEB$weight,nrow= M,ncol =6, byrow= FALSE);
				Mu  			= Betas[,3];
				Mu0 			= SimF2fEB$Intercept[1];
				
				rm(list="SimF2fEB");
				ntest 			= nrow(Basis.Test);
				#M 		= nrow(Betas);
				if(M==1 && Betas[1,1]== 0)
				{
					logL[i]  	= 0;
				}else
				{
					basisTest 	= matrix(rep(0,ntest*M),ntest,M);
					for(i_basis in 1:M){
						loc1 = Betas[i_basis,1];
						loc2 = Betas[i_basis,2];
						if(loc1==loc2){ 	basisTest[,i_basis] =  Basis.Test[,loc1];}
						else{			basisTest[,i_basis] =  Basis.Test[,loc1]* Basis.Test[,loc2];}
					}
					
					temp 		= exp(Mu0 + basisTest%*%Mu);
					if(max(temp)>1e10) temp[which(temp>1e10)] = 1e5;
					if(min(temp)<1e-10) temp[which(temp<1e-10)] = 1e-5;
					logL[i] 	= mean(Target.Test*log(temp/(1+temp)) + (1-Target.Test)*log(1/(1+temp)));
				}
				
			}#i
			Likelihood[step,] = c(alpha, lambda,mean(logL),sd(logL));
			SSE1Alpha[i_s,] 	= c(-mean(logL),sd(logL)/sqrt(nFolds));#negative logL: minimize
			currentL			= -Likelihood[step,3];
			step 				= step + 1;
			# break out of 2nd for loop
			#if((currentL - previousL)>0){break;}
		}#i_s
		index 					= which.min(SSE1Alpha[,1]);
		lambda 					= Lambda[index];
		MSEeachAlpha[i_alpha,] 	= c(alpha,lambda, -SSE1Alpha[index,]);
	}#i_alpha
	index 				= which.max(MSEeachAlpha[,3]);
	Res.lambda			= MSEeachAlpha[index,2];
	Res.alpha 			= MSEeachAlpha[index,1];
	result 				<- list(Likelihood,Res.alpha,Res.lambda);
	names(result)		<-c("CrossValidation","Alpha_optimal","Lambda_optimal");
	return(result);
}
