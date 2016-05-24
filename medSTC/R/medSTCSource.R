medSTC <-
function (documents, mlabels, ntopics, initial_c=0.5, lambda=1, rho=0.01, delta_ell=3600, supervised=TRUE, 
primal_svm=1, var_max_iter=20, convergence=1e-4, em_max_iter=100, em_convergence=1e-4, 
svm_alg_type=2, output_dir=".") 
{

	if (initial_c<=0 | lambda <=0 | rho <=0 | delta_ell <=0 | ntopics <=0 | var_max_iter<=0 | em_max_iter<=0 | convergence<=0 | em_convergence<=0) {
			cat("Numeric variables must have positive values.\n")
			return(NULL)
	}
	if(!svm_alg_type %in% c(0,2)){
		cat("svm_alg_type should be 0 or 2.\n")
		return(NULL)
	} 
		
	integerLabels<-as.integer(factor(mlabels))-1L 
    model<-list()
    class_num = length(unique(integerLabels))
    nfolds=1
    model$state<- structure(.Call("medSTCTrain", documents, integerLabels, as.integer(ntopics), as.integer(class_num), as.double(initial_c),
    							as.double(lambda),as.double(rho),as.integer(nfolds), as.double(delta_ell),
    							as.logical(supervised),as.logical(primal_svm),as.integer(var_max_iter), 
    							as.double(convergence), as.integer(em_max_iter), as.double(em_convergence),
    							as.integer(svm_alg_type), output_dir),names=c('doublePramaters','integerParameters',
    							'LogProbabilityOfWordsForTopics','Eta','Mu')
    					)
    names(model$state[[1]])<-c("DeltaEll","Lambda", "Rho", "Gamma","C","Logloss","B", "PoisOffset","Svm_Primalobj")
    names(model$state[[2]])<-c("NumberOfTopics","NumberOfLabels", "NumberOfTerms", "NumberOfDocuments")
    model$ntopics=ntopics 
    model$class_num=class_num
    model$initial_c=initial_c
    model$lambda=lambda
    model$rho=rho
    model$nfolds=nfolds
    model$delta_ell=delta_ell
    model$supervised=supervised
	model$primal_svm=primal_svm
	model$var_max_iter=var_max_iter
	model$convergence=convergence
	model$em_max_iter=em_max_iter
	model$em_convergence=em_convergence 
	model$svm_alg_type=svm_alg_type
	model$output_dir=output_dir
	model$labels = sort(unique(mlabels))
	class(model)="medSTC"
    model
}

predict.medSTC <-
function (object,documents,...) 
{
	model<-object
	result<-list()
    integerLabels<-as.integer(factor(sample(model$labels,length(documents), replace=TRUE)))-1L 
    retval<-structure(.Call("medSTCTest",model$state, documents, integerLabels, as.integer(model$ntopics), 
    	as.integer(model$class_num), as.double(model$initial_c),as.double(model$lambda),as.double(model$rho), 
    	as.integer(model$nfolds), as.double(model$delta_ell),
		as.logical(model$supervised),as.logical(model$primal_svm),as.integer(model$var_max_iter), 
		as.double(model$convergence), as.integer(model$em_max_iter), as.double(model$em_convergence),
		as.integer(model$svm_alg_type), model$output_dir))
		colnames(retval)<-model$labels
		result$scores<-retval
		result$assignments<-colnames(retval)[apply(retval,1,which.max)]
		result
}

