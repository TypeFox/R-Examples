svm.fs <- function (x, ...)
  UseMethod ("svm.fs")



`svm.fs.default` <-
				function(x,y, 
				fs.method = c("scad", "1norm", "scad+L2", "DrHSVM"),
				# chose the search method for lambda1,2: 'interval' or 'discrete'
				grid.search=c("interval","discrete"),
				### tuning parameter settings
				#fixed grid for lambda1, lambda2
				lambda1.set=NULL,  
				lambda2.set=NULL,
				# define range for lambda1,2 for interval search
				bounds=NULL, 
				# parms.coding="none" or "log2"
				parms.coding= c("log2","none"),
				 
				# internal parameter for DIRECT
				maxevals=500, 
				### valuidation settings
				# fot nested validation, 'cross.outer'-fold cv
				#cross.outer= 0,
				# method for the inner validation: cross validation, gacv   
				inner.val.method = c("cv", "gacv"),
				# 'cross.inner'-fold cv
				cross.inner= 5,
				# show plots in Direct?
				show= c("none", "final"),
				### other  settings
				# internal parameter for svm
				calc.class.weights=FALSE,
				class.weights=NULL, 
				#seed
				seed=123, 
				# max Iterations for the feature selection svm method
				maxIter=700, 
				# verbose?
				verbose=TRUE,
				...){


##  Input:
#			x: n-by-d data matrix to train (n chips/patients, d clones/genes, d>>n )
#			y: column vector (or factor vector) of target {-1, 1}'s (for n chips/patiens )
#			lam2.range : lambda2 range for elastic net (DrHSVM)

# feature selection - L1 or 1norm svm

	#require(MCRestimate)
	
	print("grid search")
	print(grid.search)
	
	print("show")
	print(show)
	
	
	require(lhs)#  - Latin Hypercube sampling function
	require(tgp)# for GP 
	require(mlegp)


	
	possible.fs<-c("1norm", "scad", "scad+L2", "DrHSVM")
	nn<-length(y) # number of cases (patients)
	nlevels.class<- nlevels(as.factor(y))
	levels.class <- levels(as.factor(y))
	
	# check labels y
	 if (!all((levels(factor(y)) == c("-1","1"))) ) stop ("labels y should be -1 and 1.")
	
	# check for 1norm
	if (fs.method == "1norm"  & inner.val.method == "gacv" )  stop("gacv is not availible for 1norm SVM. Use k fold cv instead.")
	
	#check bounds
	if (grid.search=="interval"){	
		# for  SCAD or L1 norm
			if ((fs.method %in% c("scad", "1norm" ))  & is.null(bounds)){
				bounds=t(data.frame(log2lambda1=c(-10, 10)))
				colnames(bounds)<-c("lower", "upper")	
			}
			
		# for Elastic Net
			if (fs.method %in% c("DrHSVM") & is.null(bounds)){
				bounds=t(data.frame(log2lambda2=c(-10, 10)))
				colnames(bounds)<-c("lower", "upper")	
			}
			
		
		# for Elastic SCAD
			if (fs.method %in% c("scad+L2") & is.null(bounds)){
				bounds=t(data.frame(log2lambda1=c(-10, 10), log2lambda2=c(-10,10)))
				colnames(bounds)<-c("lower", "upper")	
			}
			
			
		# for scad
		if (calc.class.weights){
			class.weights = 100/ table(y)
		}else class.weights =  NULL
	}
	
	possible.inner.val.method<-c("gacv", "cv")
	if (! (inner.val.method %in% possible.inner.val.method ))  stop(paste("You have to use one of following (inner) validation methods:", 
				paste(possible.inner.val.method, collapse=", ")))

	#checks
	if (! (fs.method %in% possible.fs ))  stop(paste("You have to use one of following fecture selection methods:", 
				paste(possible.fs, collapse=", ")))

	print(paste("feature selection method is", fs.method))
	
#	# if cross.outer>0 use outer cv !!!!!          TODO
#	if(cross.outer>0){
#		print(paste("Apply outer ", cross.outer, "-fold cross validation ", sep="" ))
#		res.cv<- .run.cv(x=x,y=y, lambda1.set=lambda1.set, lambda2.set=lambda2.set, bounds=bounds, cross.outer=cross.outer, class.weights=class.weights, seed=seed)
#		print("outer cross validation done...")
#	}   # end of if cross.outer>0 ###

  ########################################################################################
	# create final model

	# set seed again
	if (!is.null(seed)) set.seed(seed)
	
	if (grid.search=="interval") {	
	
			model<-.run.interval(x=x,y=y, fs.method=fs.method, bounds=bounds,parms.coding=parms.coding,
												 class.weights=class.weights,  maxevals=maxevals, seed=seed, 
													maxIter=maxIter, show=show, inner.val.method=inner.val.method, cross.inner=cross.inner, verbose=verbose)
	
	}		
	
	if (grid.search=="discrete") {
	
		model<-.run.discrete(x=x,y=y, fs.method=fs.method,
												lambda1.set=lambda1.set,
												lambda2.set=lambda2.set,
												parms.coding=parms.coding,
												class.weights=class.weights,  
												maxevals=maxevals, 
												seed=seed, 
												maxIter=maxIter, 
												show=show, 
												inner.val.method=inner.val.method, 
												cross.inner=cross.inner, 
												verbose=verbose)
	
	
	}
	



	
	#######################################################################################
	
	# if no outer cv is done, no correct.prediction results ;-)
#	if (cross.outer == 0) {
#		rv <- list(classes=as.factor(y),
#							sample.names = names(y),
#							class.method=paste("svm",fs.method ),
#							cross.outer=cross.outer,
#							seed = seed,
#							model =model,
#							lambda1.set=lambda1.set,
#							lambda2.set=lambda2.set,
#							bounds=bounds,
#							inner.val.method = inner.val.method,
#          		cross.inner= cross.inner,
#              cv.info = models.info.sort,
#              discrete.fit =discrete.fit   )
#    if (inner.val.method != "cv")   rv$cross.inner =NULL        
#	} else {
#	## all data are collected
#		rv <- list(votes=res.cv$vote.table,
#							classes=as.factor(y),
#							table=res.cv$confusion,
#							#sample.names = rownames(vote.matrix),
#							sample.names = names(y),
#							gene.names = colnames(x),
#							correct.prediction=res.cv$res$correct.prediction,
#							correct.class.vote=res.cv$res$correct.class.vote,
#							class.method=paste("svm",fs.method ),
#							cross.outer=cross.outer,
#							seed = seed,
#							#cross.repeat=cross.repeat,
#							#sample.names=sample.names,
#							model =model,
#							lambda1.set=lambda1.set,
#							lambda2.set=lambda2.set,
#							bounds=bounds,
#							inner.val.method = inner.val.method,
#          		cross.inner= cross.inner
#							,cv.info=res.cv$model.info.list
#							,discrete.fit =discrete.fit  
#							)
#		 if (inner.val.method != "cv")   rv$cross.inner =NULL 					
#	}

rv <- list(classes=as.factor(y),
							sample.names = names(y),
							class.method=paste("svm",fs.method ),
							#cross.outer=cross.outer,
							grid.search=grid.search,
							seed = seed,
							model =model,
							lambda1.set=lambda1.set,
							lambda2.set=lambda2.set,
							bounds=bounds,
							inner.val.method = inner.val.method,
							cross.inner= cross.inner )
 if (inner.val.method != "cv")   rv$cross.inner =NULL    

	#class(rv) <- "MCRestimate"
	class(rv) <- "penSVM"
	return(rv)
	# plot.MCRestimate
	# plot.MCRestimate ( rv, rownames.from.object=TRUE )
	
}

