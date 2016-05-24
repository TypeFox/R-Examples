#Bayes factors and posterior probabilities with BayesVarSel
#should be a model in the list of models which is nested in all the others
#(taken as the one with largest SSE)
BayesFactor<- function(models, data, prior.betas="Robust", prior.models="Constant", priorprobs=NULL){
	#cat("---------\n")
	#cat("Important note: for these results to make sense, the simplest model should be nested in all the others\n")
	#N is the number of models:
	N<- length(models)
	#n is the sample size
	n<- dim(data)[1]
	#SSE is a vector with SSE's for each model; Dim with the dimension (number of regressors in each)
	SSE<- rep(0,N); Dim<- rep(0,N)
	BFi0<- rep(0,N); PostProbi<- rep(0,N)
	#prior for betas:
	pfb<- substr(tolower(prior.betas),1,1)
	#check if the selected option exists
	if (pfb!="g" && pfb!="r" && pfb!="z" && pfb!="l" && pfb!="f") stop("I am very sorry: prior for betas no valid\n")

	#The .C to be used:
	if (pfb=="g") method<- "gBF"
	if (pfb=="r") method<- "RobustBF"
	if (pfb=="z") method<- "ZSBF"
	if (pfb=="l") method<- "LiangBF"
	if (pfb=="f") method<- "flsBF"
		

	#prior for model space:
	pfms<- substr(tolower(prior.models),1,1)
	if (pfms!="c" && pfms!="u") stop("I am very sorry: prior for models not supported\n")
		if (pfms=="u" && is.null(priorprobs)){stop("A valid vector of prior probabilities must be provided\n")}
		if (pfms=="u" && length(priorprobs)!=N){stop("Vector of prior probabilities with incorrect length\n")}
		if (pfms=="u" && sum(priorprobs<0)>0){stop("Prior probabilities must be positive\n")}
	
	#Prior probabilities of models:
	PriorModels<- rep(0,N)
	if (prior.models=="Constant"){PriorModels<- rep(1,N)}
	if (prior.models=="User"){
		#should coincide with the length of prior.models
		for (i in 1:N){PriorModels[i]<- priorprobs[[names(models)[i]]]}
	}
	
	#list that contains the names of the covariates in each model
	covar.list<- list()
	for (i in 1:N){
		temp<- lm(formula=as.formula(models[[i]]), data=data, y=TRUE, x=TRUE)
		SSE[i]<- sum(temp$residuals^2)
		Dim[i]<- length(temp$coefficients)
		covar.list[[i]]<- dimnames(temp$x)[[2]]
	}
	ordered.SSE<- sort(SSE, index.return=TRUE, decreasing=TRUE)
	#Which acts as null model:
	nullmodel<- ordered.SSE$ix[1]
	
	if (pfb!="f"){			
	 for (i in (1:N)[-nullmodel]){
		#check if the "null" model is nested in all the others
		if (sum(covar.list[[nullmodel]]%in%covar.list[[i]])<Dim[nullmodel]){stop("There is no a model nested in all the others\n")}
		Qi0<- SSE[i]/SSE[nullmodel]
		BFi0[i]<- .C(method, as.integer(n), as.integer(Dim[i]), as.integer(Dim[nullmodel]), as.double(Qi0), as.double(0.0))[5][[1]]		
	 }
    }
	
	if (pfb=="f"){
		p<- max(Dim)-min(Dim)
	 for (i in (1:N)[-nullmodel]){
		#check if the "null" model is nested in all the others
		if (sum(covar.list[[nullmodel]]%in%covar.list[[i]])<Dim[nullmodel]){stop("There is no a model nested in all the others\n")}
		Qi0<- SSE[i]/SSE[nullmodel]
		BFi0[i]<- .C(method, as.integer(p), as.integer(n), as.integer(Dim[i]), as.integer(Dim[nullmodel]), as.double(Qi0), as.double(0.0))[6][[1]]		
	 }
    }
	
	
	BFi0[nullmodel]<- 1
	names(BFi0)<- paste(names(models),".to.",names(models)[nullmodel],sep="")
	PostProbi<- BFi0*PriorModels/sum(BFi0*PriorModels)
	names(PostProbi)<- names(models)
	result<- list()
	result$BFi0<- BFi0
	result$PostProbi<- PostProbi
	result$models<- models
	result$nullmodel<- nullmodel
	class(result)<- "BayesFactor"
	result	
}


