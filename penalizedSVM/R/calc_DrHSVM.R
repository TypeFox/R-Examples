
################################	globaler!
.calc.DrHSVM<- function(parms, x.svm, y.svm,  maxIter=700, class.weights=NULL,verbose= TRUE, type="lasso", 
											inner.val.method, delta=2, cross.inner=5,
											parms.coding, seed=123, ... ){
											
										
	# SCAD	
	# input parameters are (default): log2 lambda1 and log2 lambda2 !!!!!
			
	# 1. decode the parameters ############################################################
	
	if (verbose) print("parms.coding")
	if (verbose) print(parms.coding)
	
	
	if (parms.coding=="log2") lambda2<- 2^(parms[1])
	if (parms.coding=="none") lambda2<-parms[1]
		
	# round  (8 - precision of the server!!!)
	lambda2<-round(lambda2,3)
	names(lambda2)<- NULL
	
	if (verbose) print(paste("lambda2=",lambda2))
	print(paste("maxIter=",maxIter))
	
	
	# 2. fit  model #######################################################################
	# for EPSGO, Problem = DrHSVM , output gacv = misclassification error
	
	fit <- DrHSVM(x=x.svm, y=y.svm, lambda=lambda2, type=type, trace=FALSE,  max.steps=maxIter  )
	names(fit)[names(fit) == "lambda1"] <- "lambda1.set" 			
	 							  
	# 3. calculate validation measure #####################################################
	if (inner.val.method == "gacv"){ 
		 # find model with min GACV 
   gacv.f<- min(fit$GACV)
   
   # optimal model
   fit$q.val <- gacv.f
   fit$w<- fit$beta[ fit$GACV == gacv.f, ]
   names(fit$w)<-colnames(x.svm)
   # delete w =0
   fit$w<-fit$w[fit$w!=0 ]
   fit$xind<- which(colnames(x.svm) %in% names(fit$w))
    
   fit$b<- fit$beta0[ fit$GACV == gacv.f]
   fit$lambda1<-fit$lambda1.set[ fit$GACV == gacv.f]
   fit$lambda2<- lambda2
   
   fit$q.val <-  gacv.f	
   q.val<- gacv.f
  }
	 
	if (inner.val.method == "cv"){ 
		# if the model is empty --> cv: = Inf, set to   10^16
    if (length(fit)>1){
    	# parms.coding="none" because Lambda1 is already transformed value!
    	cv.list<-.run.cv (x=x.svm, y=y.svm, fs.method="DrHSVM", cross.outer=cross.inner ,lambda1.set=lambda2, class.weights=class.weights,
    										parms.coding="none", seed=seed+1, maxIter=maxIter, verbose=verbose)
    			
    	 
      # TODO optimal lambda1 := the opt lam from cv steps with min. missclassification rate #### ??? Is it true?  
    	opt.lam1<-cv.list$set.opt.lam1[ which.max( unlist(sapply(cv.list$model.info.list,"[", "accurancy" ))   ) ]
    	   			
    	# get new coefs if the opt.lam1 is noth in the knick points of the path
    	if (!( opt.lam1  %in% fit$lambda1)){
    		predict.list<-.DrHSVM.predict(object=fit, newx=x.svm, newy=as.factor(y.svm), newlam=opt.lam1,  eps = 1e-10, verbose=verbose)
    			    		
    		fit$b<-predict.list$newbeta0
    		fit$w<-predict.list$newbeta
    		
      } else {
      	fit$b<-predict.list$newbeta0
    		fit$w<-fit$beta[ fit$lambda1.set == opt.lam1, ]
      } 		
    	
    	names(fit$w)<-colnames(x.svm)
    	# delete w =0
   		fit$w<-fit$w[fit$w!=0 ]
   		fit$xind<- which(colnames(x.svm) %in% names(fit$w))
    
   		fit$lambda1.set<-fit$lambda1
  		fit$lambda1<-opt.lam1
  		fit$lambda2<- lambda2
      						
   		cv.f<- cv.list$cv.error
    } else {
    	cv.f<- 10^16
    } 
    fit$q.val <- cv.f
    q.val <- cv.f
  }
		
	if (inner.val.method == "none") q.val<- NA 
	
	fit$inner.val.method <-inner.val.method
	
	ret<-list(q.val=q.val, model=fit)
	
	class(ret) <- "penSVM"
	
	return(ret)
}
		
		
