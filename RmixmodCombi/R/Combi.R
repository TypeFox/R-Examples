require(Rmixmod)

setClass(Class = "MixmodCombi",
		 representation = representation(mixmodOutput = "MixmodCluster",
										 hierarchy = "list",
										 ICLNbCluster = "integer"
										)
		)
		
setClass(Class = "MixmodCombiSol",
		 representation = representation(nbCluster = "integer",
		 								 partition = "integer",
		 								 proba = "matrix",
		 								 combiM = "matrix",
		 								 entropy = "numeric"		 								
		 								)
		)

setMethod(f ="print",
		  signature = c("MixmodCombi"),
		  function(x)
		  {
			if (x@mixmodOutput@dataType == "quantitative")
			{
				cat("\n EM/BIC Solution\n")
				cat(" --------------- \n\n")
				cat("Number of components: ", as.character(x@mixmodOutput@bestResult@nbCluster), "\n", sep = "") 
			
				cat("Model name: ", x@mixmodOutput@bestResult@model, "\n\n", sep="")
				for (K in 1:x@mixmodOutput@bestResult@nbCluster)
					{
						cat("Component n.", as.character(K),": ", "\n", sep="")
						cat("		proportion: ", sprintf(fmt = "%4.2f ", x@mixmodOutput@bestResult@parameters@proportions[K]), "\n", sep="")
						if (x@mixmodOutput@nbVariable == 1) cat("		mean: ", sprintf(fmt = "%4.2f  ", x@mixmodOutput@bestResult@parameters@mean[K]), "\n", sep="") else cat("		mean: ", sprintf(fmt = "%4.2f  ", x@mixmodOutput@bestResult@parameters@mean[K,]), "\n", sep="")
					}
				
				cat("\n Combining steps \n")
				cat(" --------------- \n\n")
			
				cl = paste(rep(" ", max(x@mixmodOutput@bestResult@nbCluster-4,0)), "Classes labels after this step", rep(" ", max(x@mixmodOutput@bestResult@nbCluster-4,0)), sep="")
				
				if (x@mixmodOutput@bestResult@nbCluster>4) for (K in 5:x@mixmodOutput@bestResult@nbCluster) cl = paste(" ", cl, " ", sep="")
				
				cat("  Step | Classes combined at this step | Classes labels after this step", "\n", sep="")
				cat("-------|-------------------------------|-------------------------------", "\n", sep="")
				curr = 1:x@mixmodOutput@bestResult@nbCluster
			
				cat("   0   |              ---              |", sprintf(fmt = "%2d  ", curr), "\n", sep="")
						
				for (K in 1:(x@mixmodOutput@bestResult@nbCluster-1))
					{
						Kp = x@mixmodOutput@bestResult@nbCluster - K + 1
						l1 = which(!x@hierarchy[[Kp-1]]@combiM %*% rep(1,Kp) == 1)
						l2 = (x@hierarchy[[Kp-1]]@combiM %*% curr)[l1] - curr[l1]
			
						nc1 = floor((7-nchar(as.character(K)))/2)
						nc2 = (7-nchar(as.character(K))) - nc1
						nc3 = floor((33-nchar(paste(as.character(c(l1)), "  &  ", as.character(l2))))/2) 
						nc4 = 33-nchar(paste(as.character(c(l1)), "  &  ", as.character(l2))) - nc3
						 
						curr <- x@hierarchy[[Kp-1]]@combiM %*% curr - l2*c(rep(0,(l1-1)),1,rep(0,(Kp-1-l1)))
			
						cat(rep(" ", nc1), as.character(K), rep(" ", nc2), "|", rep(" ", nc3), as.character(l1), "  &  ", as.character(l2), rep(" ", nc4), "|", sprintf(fmt = "%2d  ", curr), "\n", sep="")
					
					}
			
				cat("\n Number of clusters selected with ICL: ", x@ICLNbCluster)
				cat("\n Classification for K classes: mixmodCombiOutput@hierarchy[[K]]@partition or mixmodMap_V2M(mixmodCombiOutput@hierarchy[[K]]@partition) \n")
				cat(" Combining matrix ((K+1) classes -> K classes): mixmodCombiOutput@hierarchy[[K]]@combiM\n\n")
			}
			else if (x@mixmodOutput@dataType == "qualitative")
			{
				cat("\n EM/BIC Solution\n")
				cat(" --------------- \n\n")
				cat("Number of components: ", as.character(x@mixmodOutput@bestResult@nbCluster), "\n", sep = "") 
			
				cat("Model name: ", x@mixmodOutput@bestResult@model, "\n\n", sep="")
				for (K in 1:x@mixmodOutput@bestResult@nbCluster)
					{
						cat("Component n.", as.character(K),": ", "\n", sep="")
						cat("		proportion: ", sprintf(fmt = "%4.2f ", x@mixmodOutput@bestResult@parameters@proportions[K]), "\n", sep="")
						if (x@mixmodOutput@nbVariable == 1) cat("		center: ", sprintf(fmt = "%4.2f  ", x@mixmodOutput@bestResult@parameters@center[K]), "\n", sep="") else cat("		center: ", sprintf(fmt = "%4.2f  ", x@mixmodOutput@bestResult@parameters@center[K,]), "\n", sep="")
					}
				
				cat("\n Combining steps \n")
				cat(" --------------- \n\n")
			
				cl = paste(rep(" ", max(x@mixmodOutput@bestResult@nbCluster-4,0)), "Classes labels after this step", rep(" ", max(x@mixmodOutput@bestResult@nbCluster-4,0)), sep="")
				
				if (x@mixmodOutput@bestResult@nbCluster>4) for (K in 5:x@mixmodOutput@bestResult@nbCluster) cl = paste(" ", cl, " ", sep="")
				
				cat("  Step | Classes combined at this step | Classes labels after this step", "\n", sep="")
				cat("-------|-------------------------------|-------------------------------", "\n", sep="")
				curr = 1:x@mixmodOutput@bestResult@nbCluster
			
				cat("   0   |              ---              |", sprintf(fmt = "%2d  ", curr), "\n", sep="")
						
				for (K in 1:(x@mixmodOutput@bestResult@nbCluster-1))
					{
						Kp = x@mixmodOutput@bestResult@nbCluster - K + 1
						l1 = which(!x@hierarchy[[Kp-1]]@combiM %*% rep(1,Kp) == 1)
						l2 = (x@hierarchy[[Kp-1]]@combiM %*% curr)[l1] - curr[l1]
			
						nc1 = floor((7-nchar(as.character(K)))/2)
						nc2 = (7-nchar(as.character(K))) - nc1
						nc3 = floor((33-nchar(paste(as.character(c(l1)), "  &  ", as.character(l2))))/2) 
						nc4 = 33-nchar(paste(as.character(c(l1)), "  &  ", as.character(l2))) - nc3
						 
						curr <- x@hierarchy[[Kp-1]]@combiM %*% curr - l2*c(rep(0,(l1-1)),1,rep(0,(Kp-1-l1)))
			
						cat(rep(" ", nc1), as.character(K), rep(" ", nc2), "|", rep(" ", nc3), as.character(l1), "  &  ", as.character(l2), rep(" ", nc4), "|", sprintf(fmt = "%2d  ", curr), "\n", sep="")
					
					}
			
				cat("\n Number of clusters selected with ICL: ", x@ICLNbCluster)
				cat("\n Classification for K classes: mixmodCombiOutput@hierarchy[[K]]@partition or mixmodMap_V2M(mixmodCombiOutput@hierarchy[[K]]@partition) \n")
				cat(" Combining matrix ((K+1) classes -> K classes): mixmodCombiOutput@hierarchy[[K]]@combiM\n\n")			
			}
		  }
 		
  		 )
  		 
setMethod(f ="show",
		  signature = c("MixmodCombi"),
		  function(object)
		  {
			if (object@mixmodOutput@dataType == "quantitative")
			{
				cat("\n EM/BIC Solution\n")
				cat(" --------------- \n\n")
				cat("Number of components: ", as.character(object@mixmodOutput@bestResult@nbCluster), "\n", sep = "") 
			
				cat("Model name: ", object@mixmodOutput@bestResult@model, "\n\n", sep="")
				for (K in 1:object@mixmodOutput@bestResult@nbCluster)
					{
						cat("Component n.", as.character(K),": ", "\n", sep="")
						cat("		proportion: ", sprintf(fmt = "%4.2f ", object@mixmodOutput@bestResult@parameters@proportions[K]), "\n", sep="")
						if (object@mixmodOutput@nbVariable == 1) cat("		mean: ", sprintf(fmt = "%4.2f  ", object@mixmodOutput@bestResult@parameters@mean[K]), "\n", sep="") else cat("		mean: ", sprintf(fmt = "%4.2f  ", object@mixmodOutput@bestResult@parameters@mean[K,]), "\n", sep="")
					}
				
				cat("\n Combining steps \n")
				cat(" --------------- \n\n")
			
				cl = paste(rep(" ", max(object@mixmodOutput@bestResult@nbCluster-4,0)), "Classes labels after this step", rep(" ", max(object@mixmodOutput@bestResult@nbCluster-4,0)), sep="")
				
				if (object@mixmodOutput@bestResult@nbCluster>4) for (K in 5:object@mixmodOutput@bestResult@nbCluster) cl = paste(" ", cl, " ", sep="")
				
				cat("  Step | Classes combined at this step | Classes labels after this step", "\n", sep="")
				cat("-------|-------------------------------|-------------------------------", "\n", sep="")
				curr = 1:object@mixmodOutput@bestResult@nbCluster
			
				cat("   0   |              ---              |", sprintf(fmt = "%2d  ", curr), "\n", sep="")
						
				for (K in 1:(object@mixmodOutput@bestResult@nbCluster-1))
					{
						Kp = object@mixmodOutput@bestResult@nbCluster - K + 1
						l1 = which(!object@hierarchy[[Kp-1]]@combiM %*% rep(1,Kp) == 1)
						l2 = (object@hierarchy[[Kp-1]]@combiM %*% curr)[l1] - curr[l1]
			
						nc1 = floor((7-nchar(as.character(K)))/2)
						nc2 = (7-nchar(as.character(K))) - nc1
						nc3 = floor((33-nchar(paste(as.character(c(l1)), "  &  ", as.character(l2))))/2) 
						nc4 = 33-nchar(paste(as.character(c(l1)), "  &  ", as.character(l2))) - nc3
						 
						curr <- object@hierarchy[[Kp-1]]@combiM %*% curr - l2*c(rep(0,(l1-1)),1,rep(0,(Kp-1-l1)))
			
						cat(rep(" ", nc1), as.character(K), rep(" ", nc2), "|", rep(" ", nc3), as.character(l1), "  &  ", as.character(l2), rep(" ", nc4), "|", sprintf(fmt = "%2d  ", curr), "\n", sep="")
					
					}
			
				cat("\n Number of clusters selected with ICL: ", object@ICLNbCluster)
				cat("\n Classification for K classes: mixmodCombiOutput@hierarchy[[K]]@partition or mixmodMap_V2M(mixmodCombiOutput@hierarchy[[K]]@partition) \n")
				cat(" Combining matrix ((K+1) classes -> K classes): mixmodCombiOutput@hierarchy[[K]]@combiM\n\n")
			}
			else if (object@mixmodOutput@dataType == "qualitative")
			{
				cat("\n EM/BIC Solution\n")
				cat(" --------------- \n\n")
				cat("Number of components: ", as.character(object@mixmodOutput@bestResult@nbCluster), "\n", sep = "") 
			
				cat("Model name: ", object@mixmodOutput@bestResult@model, "\n\n", sep="")
				for (K in 1:object@mixmodOutput@bestResult@nbCluster)
					{
						cat("Component n.", as.character(K),": ", "\n", sep="")
						cat("		proportion: ", sprintf(fmt = "%4.2f ", object@mixmodOutput@bestResult@parameters@proportions[K]), "\n", sep="")
						if (object@mixmodOutput@nbVariable == 1) cat("		center: ", sprintf(fmt = "%4.2f  ", object@mixmodOutput@bestResult@parameters@center[K]), "\n", sep="") else cat("		center: ", sprintf(fmt = "%4.2f  ", object@mixmodOutput@bestResult@parameters@center[K,]), "\n", sep="")
					}
				
				cat("\n Combining steps \n")
				cat(" --------------- \n\n")
			
				cl = paste(rep(" ", max(object@mixmodOutput@bestResult@nbCluster-4,0)), "Classes labels after this step", rep(" ", max(object@mixmodOutput@bestResult@nbCluster-4,0)), sep="")
				
				if (object@mixmodOutput@bestResult@nbCluster>4) for (K in 5:object@mixmodOutput@bestResult@nbCluster) cl = paste(" ", cl, " ", sep="")
				
				cat("  Step | Classes combined at this step | Classes labels after this step", "\n", sep="")
				cat("-------|-------------------------------|-------------------------------", "\n", sep="")
				curr = 1:object@mixmodOutput@bestResult@nbCluster
			
				cat("   0   |              ---              |", sprintf(fmt = "%2d  ", curr), "\n", sep="")
						
				for (K in 1:(object@mixmodOutput@bestResult@nbCluster-1))
					{
						Kp = object@mixmodOutput@bestResult@nbCluster - K + 1
						l1 = which(!object@hierarchy[[Kp-1]]@combiM %*% rep(1,Kp) == 1)
						l2 = (object@hierarchy[[Kp-1]]@combiM %*% curr)[l1] - curr[l1]
			
						nc1 = floor((7-nchar(as.character(K)))/2)
						nc2 = (7-nchar(as.character(K))) - nc1
						nc3 = floor((33-nchar(paste(as.character(c(l1)), "  &  ", as.character(l2))))/2) 
						nc4 = 33-nchar(paste(as.character(c(l1)), "  &  ", as.character(l2))) - nc3
						 
						curr <- object@hierarchy[[Kp-1]]@combiM %*% curr - l2*c(rep(0,(l1-1)),1,rep(0,(Kp-1-l1)))
			
						cat(rep(" ", nc1), as.character(K), rep(" ", nc2), "|", rep(" ", nc3), as.character(l1), "  &  ", as.character(l2), rep(" ", nc4), "|", sprintf(fmt = "%2d  ", curr), "\n", sep="")
					
					}
			
				cat("\n Number of clusters selected with ICL: ", object@ICLNbCluster)
				cat("\n Classification for K classes: mixmodCombiOutput@hierarchy[[K]]@partition or mixmodMap_V2M(mixmodCombiOutput@hierarchy[[K]]@partition) \n")
				cat(" Combining matrix ((K+1) classes -> K classes): mixmodCombiOutput@hierarchy[[K]]@combiM\n\n")			
			}
		  }
 		
  		 )
  		 
setMethod(f ="show",
		  signature = c("MixmodCombiSol"),
		  function(object)
		  {
				cat("\n   ", as.character(object@nbCluster),"-cluster Combined Solution\n", sep = "")
				cat(" ------------------------------- ")
				cat("\n Corresponding classification: object@partition or mixmodMap_V2M(object@partition) \n")
				cat(" Combining matrix (", as.character(object@nbCluster + 1)," classes -> ", as.character(object@nbCluster)," classes):", "\n\n", sep = "")				  
				print(object@combiM)
		  }
  		 )

setMethod(f ="summary",
		  signature = c("MixmodCombi"),
		  function(object)
		  {
			if (object@mixmodOutput@dataType == "quantitative")
			{
				cat("\n EM/BIC Solution\n")
				cat(" --------------- \n\n")
				cat("Number of components: ", as.character(object@mixmodOutput@bestResult@nbCluster), "\n", sep = "") 
			
				cat("Model name: ", object@mixmodOutput@bestResult@model, "\n\n", sep="")
				for (K in 1:object@mixmodOutput@bestResult@nbCluster)
					{
						cat("Component n.", as.character(K),": ", "\n", sep="")
						cat("		proportion: ", sprintf(fmt = "%4.2f ", object@mixmodOutput@bestResult@parameters@proportions[K]), "\n", sep="")
						if (object@mixmodOutput@nbVariable == 1) cat("		mean: ", sprintf(fmt = "%4.2f  ", object@mixmodOutput@bestResult@parameters@mean[K]), "\n", sep="") else cat("		mean: ", sprintf(fmt = "%4.2f  ", object@mixmodOutput@bestResult@parameters@mean[K,]), "\n", sep="")
					}
				
				cat("\n Combining steps \n")
				cat(" --------------- \n\n")
			
				cl = paste(rep(" ", max(object@mixmodOutput@bestResult@nbCluster-4,0)), "Classes labels after this step", rep(" ", max(object@mixmodOutput@bestResult@nbCluster-4,0)), sep="")
				
				if (object@mixmodOutput@bestResult@nbCluster>4) for (K in 5:object@mixmodOutput@bestResult@nbCluster) cl = paste(" ", cl, " ", sep="")
				
				cat("  Step | Classes combined at this step | Classes labels after this step", "\n", sep="")
				cat("-------|-------------------------------|-------------------------------", "\n", sep="")
				curr = 1:object@mixmodOutput@bestResult@nbCluster
			
				cat("   0   |              ---              |", sprintf(fmt = "%2d  ", curr), "\n", sep="")
						
				for (K in 1:(object@mixmodOutput@bestResult@nbCluster-1))
					{
						Kp = object@mixmodOutput@bestResult@nbCluster - K + 1
						l1 = which(!object@hierarchy[[Kp-1]]@combiM %*% rep(1,Kp) == 1)
						l2 = (object@hierarchy[[Kp-1]]@combiM %*% curr)[l1] - curr[l1]
			
						nc1 = floor((7-nchar(as.character(K)))/2)
						nc2 = (7-nchar(as.character(K))) - nc1
						nc3 = floor((33-nchar(paste(as.character(c(l1)), "  &  ", as.character(l2))))/2) 
						nc4 = 33-nchar(paste(as.character(c(l1)), "  &  ", as.character(l2))) - nc3
						 
						curr <- object@hierarchy[[Kp-1]]@combiM %*% curr - l2*c(rep(0,(l1-1)),1,rep(0,(Kp-1-l1)))
			
						cat(rep(" ", nc1), as.character(K), rep(" ", nc2), "|", rep(" ", nc3), as.character(l1), "  &  ", as.character(l2), rep(" ", nc4), "|", sprintf(fmt = "%2d  ", curr), "\n", sep="")
					
					}
			
				cat("\n Number of clusters selected with ICL: ", object@ICLNbCluster)
				cat("\n Classification for K classes: mixmodCombiOutput@hierarchy[[K]]@partition or mixmodMap_V2M(mixmodCombiOutput@hierarchy[[K]]@partition) \n")
				cat(" Combining matrix ((K+1) classes -> K classes): mixmodCombiOutput@hierarchy[[K]]@combiM\n\n")
			}
			else if (object@mixmodOutput@dataType == "qualitative")
			{
				cat("\n EM/BIC Solution\n")
				cat(" --------------- \n\n")
				cat("Number of components: ", as.character(object@mixmodOutput@bestResult@nbCluster), "\n", sep = "") 
			
				cat("Model name: ", object@mixmodOutput@bestResult@model, "\n\n", sep="")
				for (K in 1:object@mixmodOutput@bestResult@nbCluster)
					{
						cat("Component n.", as.character(K),": ", "\n", sep="")
						cat("		proportion: ", sprintf(fmt = "%4.2f ", object@mixmodOutput@bestResult@parameters@proportions[K]), "\n", sep="")
						if (object@mixmodOutput@nbVariable == 1) cat("		center: ", sprintf(fmt = "%4.2f  ", object@mixmodOutput@bestResult@parameters@center[K]), "\n", sep="") else cat("		center: ", sprintf(fmt = "%4.2f  ", object@mixmodOutput@bestResult@parameters@center[K,]), "\n", sep="")
					}
				
				cat("\n Combining steps \n")
				cat(" --------------- \n\n")
			
				cl = paste(rep(" ", max(object@mixmodOutput@bestResult@nbCluster-4,0)), "Classes labels after this step", rep(" ", max(object@mixmodOutput@bestResult@nbCluster-4,0)), sep="")
				
				if (object@mixmodOutput@bestResult@nbCluster>4) for (K in 5:object@mixmodOutput@bestResult@nbCluster) cl = paste(" ", cl, " ", sep="")
				
				cat("  Step | Classes combined at this step | Classes labels after this step", "\n", sep="")
				cat("-------|-------------------------------|-------------------------------", "\n", sep="")
				curr = 1:object@mixmodOutput@bestResult@nbCluster
			
				cat("   0   |              ---              |", sprintf(fmt = "%2d  ", curr), "\n", sep="")
						
				for (K in 1:(object@mixmodOutput@bestResult@nbCluster-1))
					{
						Kp = object@mixmodOutput@bestResult@nbCluster - K + 1
						l1 = which(!object@hierarchy[[Kp-1]]@combiM %*% rep(1,Kp) == 1)
						l2 = (object@hierarchy[[Kp-1]]@combiM %*% curr)[l1] - curr[l1]
			
						nc1 = floor((7-nchar(as.character(K)))/2)
						nc2 = (7-nchar(as.character(K))) - nc1
						nc3 = floor((33-nchar(paste(as.character(c(l1)), "  &  ", as.character(l2))))/2) 
						nc4 = 33-nchar(paste(as.character(c(l1)), "  &  ", as.character(l2))) - nc3
						 
						curr <- object@hierarchy[[Kp-1]]@combiM %*% curr - l2*c(rep(0,(l1-1)),1,rep(0,(Kp-1-l1)))
			
						cat(rep(" ", nc1), as.character(K), rep(" ", nc2), "|", rep(" ", nc3), as.character(l1), "  &  ", as.character(l2), rep(" ", nc4), "|", sprintf(fmt = "%2d  ", curr), "\n", sep="")
					
					}
			
				cat("\n Number of clusters selected with ICL: ", object@ICLNbCluster)
				cat("\n Classification for K classes: mixmodCombiOutput@hierarchy[[K]]@partition or mixmodMap_V2M(mixmodCombiOutput@hierarchy[[K]]@partition) \n")
				cat(" Combining matrix ((K+1) classes -> K classes): mixmodCombiOutput@hierarchy[[K]]@combiM\n\n")			
			}
		  }
 		
  		 )
  		 
  		 
setMethod(f = "plot",
		  signature = c("MixmodCombi"),
		  function(x, data = x@mixmodOutput@data, dataDim = ncol(data), what = c("classification", "entropy"), reg = c(2), ...)
			{
				oldpar <- par(no.readonly = TRUE)
				on.exit(par(oldpar))
				
				if (any(what == "classification") && is.null(data)) warning("Data not supplied: classification plots can't be displayed")
					
				if (any(what == "classification") && !is.null(data)) 
				{
					## Sort z columns so that one of the two combined column is the last one at each step (prevents the colors and symbols to be mixed as K -> K-1)
						
					curr = 1:x@mixmodOutput@bestResult@nbCluster
					i <- numeric()
					j <- numeric()
						
					for (K in (x@mixmodOutput@bestResult@nbCluster):2)
					{
						l1 = which(!x@hierarchy[[K-1]]@combiM %*% rep(1,K) == 1)
						l2 = (x@hierarchy[[K-1]]@combiM %*% curr)[l1] - curr[l1]
									
						i <- c(curr[l1],i)
						j <- c(l2,j)
			
						curr <- x@hierarchy[[K-1]]@combiM %*% curr - l2*c(rep(0,(l1-1)),1,rep(0,(K-1-l1)))
					}
							
					j <- c(1,j)
					i <- c(0,i)
						
					permutIndices1 <- i
					permutIndices <- j
					permutz <- x@mixmodOutput@bestResult@proba[,permutIndices] 
						
					permutMat <- function(j,K) 
					{
						M <- diag(K)
						M[j,j] <- 0
						M[K,K] <- 0
						M[j,K] <- 1
						M[K,j] <- 1	
								
						return(M)
					}	
						
					combiM = diag(x@mixmodOutput@bestResult@nbCluster)
					currIndices = 1:x@mixmodOutput@bestResult@nbCluster
						
					## Plot
				
					if (x@mixmodOutput@dataType == "quantitative")
					{
						par(mfrow = c(dataDim, dataDim))
							
						for (K in x@mixmodOutput@bestResult@nbCluster:1)
						{
							if (dataDim > 1) for (i in 1 : dataDim^2) plot.new()
							for (i in 1 : dataDim)
							{
								if (dataDim > 1) par(mfg = c(i,i))
								histMixCluster(x, nbCluster = K, variables = colnames(data)[i], permutIndices = permutIndices, combiM = combiM, ...)
							}
							
							
							if (dataDim > 1)
							{
								for (j in 1 : (dataDim-1))
								{
									for (i in (j+1) : dataDim)
									{
										par(mfg = c(i,j))
										data = x@mixmodOutput@data
										combBestResult = x@mixmodOutput@bestResult
										combBestResult@partition = mixmodMap_M2V(mixmodMap(t(combiM %*% t(permutz))))
#										pdf(paste("~/Desktop/R_", as.character(K), ".pdf"))
										plotCluster(x = combBestResult, data = data, variable1 = colnames(data)[j], variable2 = colnames(data)[i], ...)
#										dev.off()
									} 
								}
						    }
										
						  	currTitle <- if (K == x@mixmodOutput@bestResult@nbCluster & !K == x@ICLNbCluster) paste( "BIC solution (",as.character(K)," clusters)", sep = "") else if (K == x@ICLNbCluster & !K == x@mixmodOutput@bestResult@nbCluster) paste( "Combined solution with ",as.character(K)," clusters (Number of clusters selected with ICL)", sep = "") else if (K == x@ICLNbCluster & K == x@mixmodOutput@bestResult@nbCluster) paste( "BIC and ICL solution (",as.character(K)," clusters)", sep = "") else	paste( "Combined solution with ",as.character(K)," clusters", sep = "")
#						  	currTitle <- if (K == x@mixmodOutput@bestResult@nbCluster & !K == x@ICLNbCluster) paste( "BIC solution (",as.character(K)," clusters)\n","Model: ", x@mixmodOutput@bestResult@model, "\n", sep = "") else if (K == x@ICLNbCluster & !K == x@mixmodOutput@bestResult@nbCluster) paste( "Combined solution with ",as.character(K)," clusters\n (Number of clusters selected with ICL)", sep = "") else if (K == x@ICLNbCluster & K == x@mixmodOutput@bestResult@nbCluster) paste( "BIC and ICL solution (",as.character(K)," clusters)\n","Model: ", x@mixmodOutput@bestResult@model, sep = "") else	paste( "Combined solution with ",as.character(K)," clusters", sep = "")
#						  	currTitle <- if (K == x@mixmodOutput@bestResult@nbCluster & !K == x@ICLNbCluster) paste( "BIC solution (",as.character(K)," clusters)", sep = "") else if (K == x@ICLNbCluster & !K == x@mixmodOutput@bestResult@nbCluster) paste( "Combined solution with ",as.character(K)," clusters", sep = "") else if (K == x@ICLNbCluster & K == x@mixmodOutput@bestResult@nbCluster) paste( "BIC and ICL solution (",as.character(K)," clusters)", sep = "") else	paste( "Combined solution with ",as.character(K)," clusters", sep = "")
#							currSubTitle <- if (K == x@mixmodOutput@bestResult@nbCluster & !K == x@ICLNbCluster) paste( "Model: ", x@mixmodOutput@bestResult@model, sep = "") else if (K == x@ICLNbCluster & !K == x@mixmodOutput@bestResult@nbCluster) paste( "(Number of clusters selected with ICL)", sep = "") else if (K == x@ICLNbCluster & K == x@mixmodOutput@bestResult@nbCluster) paste( "Model: ", x@mixmodOutput@bestResult@model, sep = "") else	""		
						  	par(oma=c(0,0,2,0))#, mar = c(5,4,0,2)+0.1)
#							title(main = currTitle, sub = currSubTitle, outer = TRUE)
							title(main = currTitle, outer = TRUE)
						
						   	combiM = combMat(K, which(permutIndices == permutIndices1[K]), K) %*% combiM 
						   	par(ask = TRUE)
						}
					}
					else if (x@mixmodOutput@dataType == "qualitative")
					{
						plot.new()
						
						for (K in x@mixmodOutput@bestResult@nbCluster:1)
						{
							data = x@mixmodOutput@data
							mixmodOutput = x@mixmodOutput
							mixmodOutput@bestResult@partition = mixmodMap_M2V(mixmodMap(t(combiM %*% t(permutz))))
							plot(mixmodOutput)
										
						  	currTitle <- 	if (K == x@mixmodOutput@bestResult@nbCluster & !K == x@ICLNbCluster) paste( "BIC solution (",as.character(K)," clusters). ","Model: ", x@mixmodOutput@bestResult@model, sep = "") else if (K == x@ICLNbCluster & !K == x@mixmodOutput@bestResult@nbCluster) paste( "Combined solution with ",as.character(K)," clusters\n (Number of clusters selected with ICL)", sep = "") else if (K == x@ICLNbCluster & K == x@mixmodOutput@bestResult@nbCluster) paste( "BIC and ICL solution (",as.character(K)," clusters). ","Model: ", x@mixmodOutput@bestResult@model, sep = "") else	paste( "Combined solution with ",as.character(K)," clusters", sep = "")
						  	par(oma=c(0,0,2,0), mar = c(5,4,2,2)+0.1)
							title(currTitle, outer = TRUE)
			
					   		combiM = combMat(K, which(permutIndices == permutIndices1[K]), K) %*% combiM 
					   		par(ask = TRUE)
					   	}
					}
				}
			
				if (any(what == "entropy")) par(mfrow = c(1,1)); mixmodEntPlot(z = x@mixmodOutput@bestResult@proba, combiM = sapply(x@hierarchy, function(x) x@combiM), reg = reg, ICL = x@ICLNbCluster, ...)
			
			}
		 )
		 
setMethod(
		  f="barplot",
		  signature=c("MixmodCombi"),
		  function(height, ...)
		  	{
			    	barplotMixCluster(height, ...)
			}
		)

setMethod(
		  f="hist",
		  signature=c("MixmodCombi"),
		  function(x, ...)
		  	{
			    	histMixCluster(x, ...)
			}
		)

Combi <- function(data, mixmodOutput, n = nrow(data), d = ncol(data))
{
	combiM <- list()
	combiM[[mixmodOutput@bestResult@nbCluster]] <- diag(mixmodOutput@bestResult@nbCluster)
	tau <- list()
	tau[[mixmodOutput@bestResult@nbCluster]] = mixmodOutput@bestResult@proba
	classif <- list()
	classif[[mixmodOutput@bestResult@nbCluster]] = mixmodMap_V2M(mixmodOutput@bestResult@partition, K = mixmodOutput@bestResult@nbCluster)
	entropy <- list()
	entropy[[mixmodOutput@bestResult@nbCluster]] = - sum(xlog( tau[[mixmodOutput@bestResult@nbCluster]] ))

	for (K in mixmodOutput@bestResult@nbCluster:2)
		{
			dEnt <- matrix(0,nrow=K-1, ncol=K)
			preCombiTau <- tau[[K]]
			for (l1 in 1:(K-1))
				{
					for (l2 in (l1+1):K)
						{
							postCombiTau <- t(combMat(K,l1,l2) %*% t(preCombiTau))
							dEnt[l1,l2] <- sum(xlog(postCombiTau[,l1])) - sum(xlog(preCombiTau[,l1])+xlog(preCombiTau[,l2]))
						}	
				}	
			l1=which(dEnt==max(dEnt),arr.ind=TRUE)[1]
			l2=which(dEnt==max(dEnt),arr.ind=TRUE)[2]
			
			combiM[[K-1]] <- combMat(K,l1,l2)
			tau[[K-1]] = t(combiM[[K-1]] %*% t(tau[[K]]))
			classif[[K-1]] = mixmodMap(tau[[K-1]])
			entropy[[K-1]] = - sum(xlog( tau[[K-1]] ))
			
		}
	
	# Compute ICL selection 

	ICL = numeric()
	for (k in 1:length(mixmodOutput@nbCluster)) ICL[k] = mixmodOutput@results[[k]]@criterionValue[2]
	KICL = mixmodOutput@results[[which.min(ICL)]]@nbCluster
	
	hierarchy = list()
	for (K in 1:mixmodOutput@bestResult@nbCluster) hierarchy[[K]] = new("MixmodCombiSol", nbCluster = K, partition = mixmodMap_M2V(classif[[K]]), proba = tau[[K]], combiM = combiM[[K]], entropy = entropy[[K]])
	output <- new("MixmodCombi", mixmodOutput = mixmodOutput, hierarchy = hierarchy, ICLNbCluster = KICL)
	return(output)
}
