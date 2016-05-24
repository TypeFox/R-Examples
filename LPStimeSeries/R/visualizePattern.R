visualizePattern <- function(object, x, which.terminal, orient=c(2,2)) {
  if (is.null(object$forest)) {
    stop("No forest component in ", deparse(substitute(object)))
  }

  if(!is.matrix(x)){
	if(length(x)>0){ #single time series
		x <- t(as.matrix(x))
	}
	else{
		stop("data (x) has 0 rows")
	}   
  } 
  
  if(is.null(which.terminal)){
	  stop("Terminal node info is not provided")
  }
  
  nofterminals <- c(1,apply(object$forest$nodestatus,2, function(x) sum(x==-1)))
  nofterminals <- cumsum(nofterminals)
  
  if(which.terminal > nofterminals[object$ntree+1]){
	  stop(paste("Total number of terminal nodes is",nofterminals[object$ntree+1],"which is less than",which.terminal))
  }

  whichtree <- findInterval(which.terminal,nofterminals)
  terminal <- which.terminal - nofterminals[whichtree] + 1
  mdim <- ncol(x)
  ntest <- nrow(x)
  x <- t(data.matrix(x))

  keepIndex <- c("predictpatterns", "targetpatterns")
  ans <- .C("regForest_pattern",  
			as.double(x),
			as.integer(ntest),
			as.integer(whichtree),
			as.integer(terminal),
			as.double(object$segment.length),
			as.integer(mdim),
			as.integer(object$ntree),
			object$forest$leftDaughter,
			object$forest$rightDaughter,
			object$forest$nodestatus,
			object$forest$nodedepth,
			object$forest$nrnodes,
			object$forest$xbestsplit,
			object$forest$bestvar,
			object$forest$splitType,
			object$forest$ndbigtree,
			as.integer(object$maxdepth),
			as.integer(object$target),
			as.integer(object$target.type),
			predictpatterns = double(ntest * mdim),
			targetpatterns = double(ntest * mdim),
		PACKAGE = "LPStimeSeries")[keepIndex]
	
	ans$targetpatterns[ans$targetpatterns==-999] <- NA		
	resT <- t(matrix(ans$targetpatterns,nrow=mdim))

	ans$predictpatterns[ans$predictpatterns==-999] <- NA		
	res <- t(matrix(ans$predictpatterns,nrow=mdim))

	bu <- apply(res,1,function(x) sum(!is.na(x)))
	ind <- order(-bu)

    lims=max(abs(x))
    par(mfrow=orient)
    for(k in 1:prod(orient)){
		plot(c(1:nrow(x)),c(1:nrow(x)),main=paste('Pattern',k,'from series',ind[k]),col=0,
			 ylim=c(-1.05*lims,1.05*lims),ylab='',xlab='')
		points(resT[ind[k],])
		points(res[ind[k],],col=2,pch=2)
		legend("topleft",c("Target","Predictor"),pch=c(1,2),col=c(1,2))
	}

	out=list(predictor=res,target=resT,tree=whichtree,terminal=terminal)
	invisible(out)
}
