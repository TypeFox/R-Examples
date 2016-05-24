
##################################################################
# with function for a nested imputation list
with.NestedImputationList <- function(data, expr,fun,...){
  pf<-parent.frame()
  NB <- data$Nimp["Between"]  
  results0 <- as.list( seq( 1 , NB) )
  data0 <- data
  for (bb in 1:NB){
	  data <- list()
	  data$imputations <- data0$imputations[[bb]]
	  #**************************************
	  # This code is copied from the with.imputationList
	  # function from the mitools package
	  if ( ! is.null(match.call()$expr) ){
		expr <- substitute(expr)
		results <- lapply(data$imputations, function(dataset){
						eval(expr, dataset, enclos=pf)
						} )
				} else {
		results <- lapply(data$imputations, FUN=fun ,...)
	  }

	  if (all(sapply(results, inherits,  what="imputationResult"))){
		class(results)<-"imputationResultList"
		results$call<-sys.call(-1)
	  } else {
		attr(results,"call")<-sys.call(-1)
	  }
	  
		results0[[bb]] <- results
		#*************************************************
		}
    class(results0) <- "NestedImputationResultList"
    return(results0) 
}
####################################################################