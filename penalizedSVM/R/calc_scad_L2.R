
################################	globaler!
.calc.scad_L2<- function(parms, x.svm, y.svm,  maxIter=700, class.weights=NULL,verbose= TRUE, inner.val.method, cross.inner=5,
											parms.coding, seed=123, ... ){
	# SCAD	
	# input parameters are (default): log2 lambda1 and log2 lambda2 !!!!!
		
	# 1. decode the parameters ############################################################
	
	if (verbose) print("parms.coding")
	if (verbose) print(parms.coding)
	
	
	if (parms.coding=="log2") parms<- 2^(parms)
	if (parms.coding=="none") parms<-parms
		
	# decode the parameters
	lambda1<- parms[1]
	lambda2<- parms[2]
	# round 
	lambda1<-signif(lambda1,3)
	lambda2<-signif(lambda2,3)
	names(lambda1)<-names(lambda2)<- NULL
	if (verbose) {
			print(paste("lambda1=",lambda1))
		 	print(paste("lambda2=",lambda2))
	}
	print(paste("maxIter=",maxIter))
	
	# 2. fit  model #######################################################################
	 fit<- scad_L2.svc(lambda1 = lambda1, lambda2=lambda2, x.svm, y.svm, a = 3.7, maxIter=maxIter, verbose=verbose)
	
	# 3. calculate validation measure #####################################################
	if (inner.val.method == "gacv"){ 
		# if the model is empty --> gacv: = Inf, set to   10^16
    gacv.f<- ifelse (length(fit)>1, findgacv.scad (y.svm, model=fit), 10^16)
    gacv.f<-max(gacv.f, 0 )
    if (length(fit)>1) {
    	fit$q.val <- gacv.f
    	fit$inner.val.method =inner.val.method
  	}
  	q.val<- gacv.f
  }
	 
	if (inner.val.method == "cv"){ 
		# if the model is empty --> gacv: = Inf, set to   10^16
    if (length(fit)>1){
    	# parms.coding="none" because Lambda1 is already transformed value!
    	
    	# if the mode for the whole data exists and we get troubles by cv (reduced data), change seed and try again
    	try.i<-1
    	
    	if (try.i <=5){
    		if (verbose) print(paste("try.i=", try.i))
    		if (exists("cv.list")) rm(cv.list)
    		try( cv.list<-.run.cv (x=x.svm, y=y.svm, fs.method="scad+L2", cross.outer=cross.inner ,lambda1.set=c(lambda1,lambda2), 
    										class.weights=class.weights,
    										parms.coding="none", seed=seed+1,maxIter=maxIter, verbose=verbose))
    		if (!exists("cv.list")) {
    			if (verbose) print("error by cross validation of Elastic Net, try another speed")
    			try.i<-try.i + 1; seed<-seed+1
    		}else try.i<-10
    	}									
  	
   		if (exists("cv.list")) {cv.f<- cv.list$cv.error}
   		else {cv.f<- 10^16 }
    		
    } else {
    	cv.f<- 10^16
    }
    
     
    fit$q.val <- cv.f
    q.val <- cv.f
  }
	
	# skip xqx
	fit<-fit[- which(names(fit) == "xqx")]
	
	if (inner.val.method == "none") q.val<- NA 
	
	
	ret<-list(q.val=q.val, model=fit)
	
	class(ret) <- "penSVM"
	
	return(ret)
}