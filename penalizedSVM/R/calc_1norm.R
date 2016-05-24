
################################	globaler!
.calc.1norm<- function(parms, x.svm, y.svm,  maxIter=700, class.weights=NULL,verbose= TRUE, inner.val.method="cv", cross.inner=5,
											parms.coding, seed=123,  ... ){
	# 1NORM	:  only cv! 
	# input parameters are (default): log2 lambda1 !!!!
		
	# 1. decode the parameters ############################################################
	
	if (verbose) print("parms.coding")
	if (verbose) print(parms.coding)
		
	if (parms.coding=="log2") lambda1<- 2^(parms[1])
	if (parms.coding=="none") lambda1<-parms[1]
		
	# round  (8 - precision of the server!!!)
	lambda1<-round(lambda1,3)
	names(lambda1)<- NULL
	
	if (verbose) print(paste("lambda1=",lambda1))
	if (verbose) print(paste("maxIter=",maxIter))
			
	catch.error<-function()ifelse ( length(grep("  from Lapack routine",geterrmessage() ) )>0, print("internal error of  Lapack routine 'dgesdd'"),print("something else")  )
	options(show.error.messages=TRUE,  error=catch.error)
		
	# 2. fit  model #######################################################################
	if (inner.val.method == "gacv") {
		 inner.val.method<-"cv"
		 if (verbose) print("gacv is not availible for 1norm SVM. Use 5 fold cv instead.")
	}
	
	if (verbose) print(paste("start", cross.inner, "fold cross validation"))
	#  cross.inner fold cv
	# fit$testCorr = correctness in % [0,100],  rewrite as  1- (testCorr)/100 
	if (exists("fit")) rm(fit)
	try(fit <- lpsvm(A=x.svm, d=y.svm, k=cross.inner, nu=0,output=0, delta=10^-3, epsi=lambda1, seed=seed))
	cv<-ifelse ( (exists("fit") & !is.null(fit$testCorr)), 1- (fit$testCorr)/100 , 10^16 )	
	#if fit does not exists
	if (!exists("fit")) fit<- NULL
	
	# choose the lam1(=epsi) with max testing set correctness:  testCorr
	# revert to default
	options(error = NULL)
	
	ret<-list(q.val=cv, model=fit)
	class(ret) <- "penSVM"
	
	return(ret)
}
		
		
		
