.run.discrete<- function(x,y, fs.method,
						lambda1.set,
						lambda2.set,
						parms.coding,
						class.weights,  
						maxevals, 
						seed=123, 
						maxIter, 
						show, 
						inner.val.method, 
						cross.inner, 
						verbose=TRUE	) {
						
						
# (grid.search=="discrete") 
							
		### SCAD #####
		if (fs.method == "scad" ){
			# fine grid at the begin of the interval !!!
			#if (is.null(lambda1.set)) lambda1.set <- c (0.001, 0.01,  seq(0.1 ,1, 0.1),5)
			if (is.null(lambda1.set)) lambda1.set <- 2^(-10:5)
		  lambda.points<-data.frame(lambda1.set)
		  Q.func<-".calc.scad"
		}
		
		### L1 norm #####
		if (fs.method == "1norm" ){
			if (is.null(lambda1.set)) lambda1.set <- 2^(-10:5)
		  lambda.points<-data.frame(lambda1.set)
		  Q.func<- ".calc.1norm"
		}
		
		### Elastic Net  #####
		if (fs.method == "DrHSVM" ){
			if (is.null(lambda2.set)) lambda2.set <- 2^(-10:5)
		  lambda.points<-data.frame(lambda2.set)
		  Q.func<- ".calc.DrHSVM"
		}
		
		### Elastic SCAD
		if (fs.method == "scad+L2" ){
			if (is.null(lambda1.set)) lambda1.set <- 2^(-10:5)
			if (is.null(lambda2.set)) lambda2.set <- 2^(-10:5)
			
			# write down all points of the grid
			lambda.points<-data.frame("lam1"=rep(lambda1.set, length(lambda2.set)),
			 													"lam2"=rep(lambda2.set, each=length(lambda1.set))	 )
			Q.func<-".calc.scad_L2" 													
		}		


		# calculate models in each point on the grid
		if (verbose)  print(paste("discrete grid search, calculate models in", nrow(lambda.points), "points")) 
		
		models<-apply(lambda.points,1, eval(Q.func) , 
	                  x.svm=x, y.svm=y, maxIter=maxIter,
	                  class.weights=class.weights,
	                  parms.coding=parms.coding, 
	                  inner.val.method=inner.val.method, 
	                  cross.inner = cross.inner,
	                  verbose=verbose,
	                  seed = seed)
		
		#save(models, file=paste("discrete_models_",fs.method, ".RData", sep=""))
		#print("/////////////////////////////////////////////////")
			
		print(str(models))
	  
	  # delete empty models: model: list() or model : "No variable selected."
	  flag.empty<- (sapply(sapply(models, "[",2 ),"[",1 ) == "No variable selected.") | (sapply(sapply(sapply(models, "[",2 ),"[",1 ),is.null))
	  if (verbose) print(paste("out of ", nrow(lambda.points), "models,", sum(flag.empty), "are empty models. Skip them" ) )
	   
	  models<-models[!flag.empty]
	  
	  # if any model exists...
	  if (length(models)>=1){
		  models.info<-data.frame("lambda1"=as.numeric(sapply(
		                          								sapply(models, "[",2 ), 
		                          								"[",
		                          								which(names(sapply(models[1], "[",2 )[[1]])=="lambda1")
		                          								)), 
		  												"lambda2"=as.numeric(sapply(
		                          								sapply(models, "[",2 ), 
		                          								"[",
		                          								which(names(sapply(models, "[",2 )[[1]])=="lambda2")
		                          								)),
		                          "num.w"= sapply(sapply(
		                          								sapply(models, "[",2 ), 
		                          								"[",
		                          								which(names(sapply(models, "[",2 )[[1]])=="w")
		                          								), length), 
		                          "q.val"=as.numeric(sapply(models, "[",1)), 
		                          inner.val.method = inner.val.method)
		                          
			
			
			### find optimal model
			#
			# sort models.info first by q.val and then by num.w ! 
			models.info.sort<- sortmat( models.info, c(which(colnames(models.info) == "q.val"),which(colnames(models.info) == "num.w")) )
			# example: # sortmat( models.info, c(which(colnames(models.info) == "q.val"),which(colnames(models.info) == "num.w")) )
			#  lambda1 num.w q.val inner.val.method
			#4    0.50    10 0.345               cv
			#2    0.03    64 0.345               cv
			#3    0.40     3 0.350               cv
			#1    0.02    74 0.350               cv
				
			# not the same as   sortmat( models.info, c(which(colnames(models.info) == "num.w"),which(colnames(models.info) == "q.val")) )        !!!!
			#  lambda1 num.w q.val inner.val.method
			#3    0.40     3 0.350               cv
			#4    0.50    10 0.345               cv
			#2    0.03    64 0.345               cv
			#1    0.02    74 0.350               cv
				
			if (verbose) {
			  print("optimal model is:") 
			  print(models.info.sort[1,])
			}  
			
			# if only one non-empty model ...
			if (nrow(models.info.sort) == 1 ) rownames(models.info.sort) <- "1" 
				
			model<- models[[as.numeric(rownames(models.info.sort)[1])]]$model
			print(str(model))
			
			# store all models wih min. q.val
			q.min<- min(models.info.sort$q.val)	
			
			discrete.fit = models[as.numeric(rownames(models.info.sort)[which(models.info.sort$q.val == q.min)] ) ]
			discrete.fit<-sapply(discrete.fit, "[", 2 )
	  } else{ # end of if length(models)>1
    	models.info.sort<- NULL
    	model<- NULL
    	discrete.fit<- NULL
    }
				
		f.final<- model
		f.final$cv.info <- models.info.sort # model info, all possible candidates (much space !) 
		f.final$fit.info<-discrete.fit
								
						
						
return(f.final) 						
					
}
