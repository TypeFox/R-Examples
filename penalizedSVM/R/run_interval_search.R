
.run.interval<- function(x,y, bounds, 
																parms.coding="none", # or log2 
																fs.method,
																class.weights, 
																seed=123,
																#number of initial points
																N=NULL,
																maxevals = 2000, 
																# minimal value for gacv =0
																fminlower = 0,
																pdf.name= NULL,  
																maxIter= NULL, 
																inner.val.method="cv",
																verbose=TRUE, 
																# for DIRECT show the points 
																show="final",
																cross.inner=5,
																 ...){
	# for  "SCAD + L2" 
	#  do the interval search by applying EPSGO Algorithm 
	# find several points in parameter space with the same min value,
	# chose the point with min # FS in SVM model
	
	#check bounds 
	
	
	# Function measuring the gacv for Elastic SCAD model
	
	if (fs.method== "1norm") 	Q.func<- ".calc.1norm"
	if (fs.method== "scad") 	Q.func<- ".calc.scad"
		
	if (fs.method== "scad+L2") 	Q.func<- ".calc.scad_L2"
	if (fs.method== "DrHSVM") 	Q.func<- ".calc.DrHSVM"
	
	if (verbose) print("start interval search")
	if (verbose) print(paste("inner validation method:", inner.val.method ))
	fit<-EPSGO(Q.func, bounds=bounds, parms.coding=parms.coding, fminlower=fminlower, show=show, N=N,  maxevals=maxevals, 
						 pdf.name=pdf.name,  seed=seed,  
						 verbose=verbose,
						 # Q.func specific parameters:
						 x.svm=x, y.svm=y, class.weights=class.weights, maxIter=maxIter, inner.val.method=inner.val.method,
						 cross.inner=cross.inner )

	# if we have several models mit minimal gacv value, chose those with min number of features (genes)
	# by equal number of genes --> take the first modell 
		
	# take the calculated final model
	if (fs.method != "1norm"){
		
		if (verbose) print("# FS: ")
		sel.models<-  sapply(fit$model.list, "[", "model") [fit$Ytrain == fit$fmin ]
		
		
		len.w<-sapply( sapply(sel.models, "[", "w") , length )
		
		if (verbose) 	print("chose the model with min num of FS ")
		opt.model<- sel.models[[ which.min(len.w)]]
		
		# rewrite fit$model.list
		fit$model.list<-sel.models
		
	}else{ # for L1 SVM need to calculate the final model, without 5-fold cv!
		if (verbose) print("chose the model with min num of FS ")
		
		catch.error<-function()ifelse ( length(grep("  from Lapack routine",geterrmessage() ) )>0, print("internal error of  Lapack routine 'dgesdd'"),print("something else")  )
		options(show.error.messages=TRUE,  error=catch.error)
	
		# parameter in fit is log2value, --> epsi=2^param !
		sel.epsis<- 2^as.numeric(fit$Xtrain [fit$Ytrain == fit$fmin , ])
				
		sel.final.models<-list()
		for (ep in sel.epsis){
			ep.model<- NULL
			try(ep.model <- lpsvm(A=x, d=y, k=0, nu=0,output=0, delta=10^-3, epsi=ep, seed=seed))
	    
	    sel.final.models<-c(sel.final.models, list(q.val=fit$fmin, model=ep.model))
	  }
		
		# Number  of FS 
		tmp.models<-sel.final.models[which(names(sel.final.models) == "model")]
		len.w<-sapply( sapply( tmp.models , "[", "w") , length )
		
		
		# revert to default
		options(error = NULL)
		
		opt.model<- tmp.models[[ which.min(len.w)]]
		opt.model$lambda1<-opt.model$epsi
		
		# rewrite fit$model.list
		fit$model.list<-sel.final.models
	}
	
	f.final<- list(w=		 			opt.model$w,
								 b= 	 			opt.model$b, 
								 xind=			opt.model$xind, 
								 index= 		opt.model$index,
								 fitted=		opt.model$fitted, 
								 type= 			opt.model$type,
								 lambda1=		opt.model$lambda1, 
								 lambda2= 	opt.model$lambda2,
								 iter = 		opt.model$iter,
								 q.val= 			fit$fmin 
								 ,  fit.info = fit ) # model info, all possible candidates (much space !) 
									#)
	
	return(f.final)
	
}
	