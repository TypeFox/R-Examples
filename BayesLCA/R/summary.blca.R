summary.blca <-
function(object, ...){
	sum0<-sum2<-NULL
	sum2$ItemProb<- list(alpha=object$prior$alpha, beta=object$prior$beta)
	sum2$ClassProb<- object$prior$delta
	if(is.matrix(sum2$ItemProb$beta)){
	rownames(sum2$ItemProb$beta)<- rownames(sum2$ItemProb$alpha)<- names(sum2$ClassProb)<- paste("Group", rep(1:length(sum2$ClassProb)))
	} else names(sum2$ClassProb) <- "Group 1"
	
	if(is.matrix(sum2$ItemProb$alpha)){
	  if(is.null(colnames(object$itemprob))){ 
	    colnames(sum2$ItemProb$beta)<- colnames(sum2$ItemProb$alpha)<- paste("Col", rep(1:ncol(object$itemprob)))
	    } else colnames(sum2$ItemProb$beta)<- colnames(sum2$ItemProb$alpha)<- colnames(object$itemprob)
	    } else if(is.null(colnames(object$itemprob))){
	      names(sum2$ItemProb$beta)<- names(sum2$ItemProb$alpha)<- paste("Col", rep(1:ncol(object$itemprob)))
	      } else names(sum2$ItemProb$beta)<- names(sum2$ItemProb$alpha)<- colnames(object$itemprob)
	

	sum0$method<- object$method
	sum0$printnames<- object$printnames
	sum0$sum1<- object$sum1
	sum0$sum2<- sum2
	class(sum0)<-"summary.blca"
	return(sum0)
	}
