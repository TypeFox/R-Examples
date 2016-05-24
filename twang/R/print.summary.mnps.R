# Produces a summary table for ps object 
print.summary.mnps <- function(x, ...)
{
      dots <- list(...)
      
      if(is.null(dots$pairwiseMax)) dots$pairwiseMax <- TRUE


      if(x$estimand == "ATE"){
      	cat("Summary of pairwise comparisons:\n")
      	print(x$comp)
      	cat("\nSample sizes and effective sample sizes:\n")
      	print(x$ess)
      		
      }
      if(x$estimand == "ATT"){
      	nSum <- length(x$summaryList)      		
      	if(!dots$pairwiseMax){
      		cat("Summary of mnps object:\n")

      	for(i in 1:nSum){
      		cat("Summary of observations receiving treatment ", x$levExceptTreatATT[i], " weighted to match the observations receiving treatment ", x$treatATT, ".\n", sep = "")
      	   	if(!is.null(dots$digits)) obj <- round(x$summaryList[[i]][,-c(1,2,3,8)], digits = digits)
      		else obj <- x$summaryList[[i]][,-c(1,2,3,8)]
		    obj <- data.frame(obj)	
#      class(obj) <- "matrix"
	      	names(obj) <- c("ESS", "max.es","mean.es","max.ks","mean.ks","iter")

      print(obj)
      cat("\n")
      }
      }
      else{
      	cat("Summary of pairwise comparisons:\n")
      	obj <- as.data.frame(x$summaryList[[1]][,c("max.es","max.ks")])
      	for(i in 2:nSum){
      		obj$max.es <- apply(cbind(obj$max.es, x$summaryList[[i]][,"max.es"]), 1, max, na.rm = TRUE)
      		obj$max.ks <- apply(cbind(obj$max.ks, x$summaryList[[i]][,"max.ks"]), 1, max, na.rm = TRUE)
      	}
      	
      	obj <- data.frame(obj)
      	names(obj) <- c("max.es","max.ks")
      	print(obj)
      	
      	cat("Sample sizes and effective sample sizes:\n")
		print(x$ess)      	
      	
      	}
      }
      invisible(x)
      }

