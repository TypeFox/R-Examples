
################################	globaler!
.calc.scad<- function(parms, x.svm, y.svm,  maxIter=700, class.weights=NULL,verbose= TRUE, inner.val.method, cross.inner=5,
											parms.coding, seed=123, ... ){
	# SCAD	
	# input parameters are (default): log2 lambda1 and log2 lambda2 !!!!!
		
	# 1. decode the parameters ############################################################
	
	if (verbose) print("parms.coding")
	if (verbose)  print(parms.coding)
	
	
	if (parms.coding=="log2") lambda1<- 2^(parms[1])
	if (parms.coding=="none") lambda1<-parms[1]
		
	# round  (8 - precision of the server!!!)
	lambda1<-round(lambda1,3)
	names(lambda1)<- NULL
	
	if (verbose) print(paste("lambda1=",lambda1))
	if (verbose) print(paste("maxIter=",maxIter))
	
	# 2. fit  model #######################################################################
	fit<- scadsvc(lambda1 = lambda1,  x=x.svm, y=y.svm, class.weights= class.weights, maxIter=maxIter, verbose=verbose) 
	
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
    	cv.list<-.run.cv (x=x.svm, y=y.svm, fs.method="scad", cross.outer=cross.inner ,lambda1.set=lambda1, class.weights=class.weights,
    										parms.coding="none", seed=seed+1, maxIter=maxIter, verbose=verbose)
   		cv.f<- unlist(cv.list$cv.error)
    } else {
    	cv.f<- 10^16
    } 
        
    fit$q.val <- unlist( cv.f)
    q.val <- unlist( cv.f)
  }
	
	# skip xqx
	if ("xqx" %in% names(fit) ) fit<-fit[- which(names(fit) == "xqx")]
	
	if (inner.val.method == "none") q.val<- NA 
	fit$inner.val.method<- inner.val.method
	
	#return(fit$gacv.f)
	ret<-list(q.val=q.val, model=fit)
	
	class(ret) <- "penSVM"
	
	return(ret)
}
		