cv.tune <-
function(train, numgrid = 20, classifier = "snn"){
	# classifier can be "knn", "bnn", "ownn", or "snn"
	# tune knn, bnn, ownn by minimizing CV error
	# tune snn w.r.t lambda by a two-stage procedure: step 1: lower 10% quantile of error; step 2: min CIS 
	
	n = dim(train)[1] 
	d = dim(train)[2]-1
	nfolds = 5                                 
	set.seed(1)
	case.folds <- sample(rep(1:nfolds,length.out=n))

	if(classifier == "knn" | classifier == "ownn"){
	
		kmin = 1																	#minimal k = 1 for knn
		kmax = floor(n/2) 															#maximal k = n/2 for knn
		kvec = floor(seq(kmin,kmax,length=numgrid))									#interval of k for knn
		error_knn_cv = matrix(0,numgrid,nfolds)
		
		for(kk in 1:numgrid){
			for(fold in 1:nfolds){		
				traincv <- train[case.folds!=fold,]
				testcv <- train[case.folds==fold,]
			
				KNN= myknn(traincv,testcv[,1:d],kvec[kk])
				error_knn_cv[kk,fold]= (length(which(KNN!=testcv[,d+1])))/length(testcv[,d+1])	
			}
		}

		kopt = floor((kvec[which.min(apply(error_knn_cv,1,mean))[1]])* (5/4)^(4/(d + 4)))   #rescale optimal k for knn w.r.t. sample size, see P2750 of S12.
		parameter.opt = kopt 
		parameter.list = kvec
		
	}else if(classifier == "bnn"){
	
		kmin = 1																	#minimal k = 1 for knn
		kmax = floor(n/2) 															#maximal k = n/2 for knn
		kvec = floor(seq(kmin,kmax,length=numgrid))									#interval of k for knn
		error_knn_cv = matrix(0,numgrid,nfolds)
	
		for(kk in 1:numgrid){
			for(fold in 1:nfolds){		
				traincv <- train[case.folds!=fold,]
				testcv <- train[case.folds==fold,]
			
				KNN= myknn(traincv,testcv[,1:d],kvec[kk])
				error_knn_cv[kk,fold]= (length(which(KNN!=testcv[,d+1])))/length(testcv[,d+1])	
			}
		}

		kopt = floor((kvec[which.min(apply(error_knn_cv,1,mean))[1]])* (5/4)^(4/(d + 4)))   #rescale optimal k for knn w.r.t. sample size, see P2750 of S12.
		qopt = 2^(d/(d+4))*(gamma(2+2/d))^(2*d/(d+4))/kopt          					 	#optimal resampling ratio for bagged nn, see P2746 of S12.
		parameter.opt = qopt 
		parameter.list = 2^(d/(d+4))*(gamma(2+2/d))^(2*d/(d+4)) / kvec

	}else if(classifier == "snn"){
	
		kmin = 1																	#minimal k = 1 for knn
		kmax = floor(n/2) 															#maximal k = n/2 for knn
		kmin_snn = floor(kmin*(5/4)^(4/(d+4))*((2*d+8)/(d+2))^(d/(d+4)))   			#kmin/kmax of snn equals kmin/kmax of ownn
		kmax_snn = floor(kmax*(5/4)^(4/(d+4))*((2*d+8)/(d+2))^(d/(d+4)))
		lammin = (kmin_snn*((2*d+4)/(d*(d+4)))^(d/(d+4))*n^(-4/(d+4)))^((d+4)/d)    					
		lammax = (kmax_snn*((2*d+4)/(d*(d+4)))^(d/(d+4))*n^(-4/(d+4)))^((d+4)/d)							
		Lambda = (seq(lammin^(d/(d+4)),lammax^(d/(d+4)),length=numgrid))^((d+4)/d)  #interval of lambda, guarantee equally space in k.
		
		cis_snn_cv = matrix(0,numgrid,nfolds)
		error_snn_cv = matrix(0,numgrid,nfolds)
		
		for(kk in 1:numgrid){
			for(fold in 1:nfolds){		
				traincv <- train[case.folds!=fold,]
				testcv <- train[case.folds==fold,]
			
				SNN= mysnn(traincv,testcv[,1:d],Lambda[kk])
				error_snn_cv[kk,fold]= (length(which(SNN!=testcv[,d+1])))/length(testcv[,d+1])
			
				set.seed(fold)
				index1 = sample(1:nrow(traincv),floor(nrow(traincv)/2))
				SNN1 = mysnn(traincv[index1,],testcv[,1:d],Lambda[kk])
				SNN2 = mysnn(traincv[-index1,],testcv[,1:d],Lambda[kk])
				cis_snn_cv[kk,fold]= (length(which(SNN1!=SNN2)))/length(testcv[,d+1])
			}
		}

		error_snn = apply(error_snn_cv,1,mean)
		cis_snn = apply(cis_snn_cv,1,mean)
		active_set = which(round(error_snn,6) <= round(quantile(error_snn,0.1),6))				# step 1: lower 10% quantile of error; step 2: min CIS
		lamopt = Lambda[active_set[which.min(cis_snn[active_set])]]
		parameter.opt = lamopt 
		parameter.list = Lambda
	
	}else{

		stop("The argument classifier can only be one of: knn, bnn, ownn, snn!")
	
	}

	return(list(parameter.opt = parameter.opt, parameter.list = parameter.list))

}
