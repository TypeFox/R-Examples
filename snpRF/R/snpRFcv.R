snpRFcv <- function(trainx.autosome=NULL,trainx.xchrom=NULL,trainx.covar=NULL, trainy, cv.fold=5, scale="log", step=0.5,
                 mtry=function(p) max(1, floor(sqrt(p))), recursive=FALSE,
                 ...) {
    
    n <- as.numeric(as.character(unique(c(nrow(trainx.autosome),nrow(trainx.xchrom),nrow(trainx.covar)))))
    if(length(n)==0) stop("all trainx.* data missing")
    if(length(n)>1) stop("all trainx.* data do not agree in dimension")

    ### strip names from trainx.* matrices ###

    if(!is.null(trainx.autosome)) dimnames(trainx.autosome)[[2]]<-paste("A",1:dim(trainx.autosome)[2],sep="")
    if(!is.null(trainx.covar)) dimnames(trainx.covar)[[2]]<-paste("A",1:dim(trainx.covar)[2],sep="")
    if(!is.null(trainx.xchrom)) xchrom.names<-paste("X",1:(dim(trainx.xchrom)[2]/2),sep="")

    p <- sum(c(ncol(trainx.autosome),ncol(trainx.xchrom)/2,ncol(trainx.covar)),na.rm=T) 
    if (scale == "log") {
        k <- floor(log(p, base=1/step))
        n.var <- round(p * step^(0:(k-1)))
        same <- diff(n.var) == 0
        if (any(same)) n.var <- n.var[-which(same)]
        if (! 1 %in% n.var) n.var <- c(n.var, 1)
    } else {
        n.var <- seq(from=p, to=1, by=step)
    }
    k <- length(n.var)
    cv.pred <- vector(k, mode="list")
    for (i in 1:k) cv.pred[[i]] <- trainy
    ## Generate the indices of the splits
    ## Stratify the classes for classification problem.
    f <- trainy
   
    nlvl <- table(f)
    idx <- numeric(n)
    for (i in 1:length(nlvl)) {
        idx[which(f == levels(f)[i])] <-  sample(rep(1:cv.fold, length=nlvl[i]))
    }

    for (i in 1:cv.fold) {
        ## cat(".")

	
        all.rf <- snpRF(x.autosome=trainx.autosome[idx != i, , drop=FALSE],
	       	  	x.xchrom=trainx.xchrom[idx != i, , drop=FALSE],
	       	  	x.covar=trainx.covar[idx != i, , drop=FALSE],
                        y=trainy[idx != i],
                        xtest.autosome=trainx.autosome[idx == i, , drop=FALSE],
                        xtest.xchrom=trainx.xchrom[idx == i, , drop=FALSE],
                        xtest.covar=trainx.covar[idx == i, , drop=FALSE],
                        ytest=trainy[idx == i],
                        mtry=mtry(p), importance=TRUE, ...)

		
        cv.pred[[1]][idx == i] <- all.rf$test$predicted
        impvar <- dimnames(all.rf$importance)[[1]][order(all.rf$importance[,1], decreasing=TRUE)]
        for (j in 2:k) {
            imp.idx <- impvar[1:n.var[j]]
	    if(any(imp.idx %in% dimnames(trainx.autosome)[[2]])){
	    	       sub.autosome<-trainx.autosome[,imp.idx[imp.idx %in%  dimnames(trainx.autosome)[[2]]],drop=F]
            }else{
		       sub.autosome<-NULL
	    }
	    if(any(imp.idx %in% dimnames(trainx.covar)[[2]])){
	    	       sub.covar<-trainx.covar[,imp.idx[imp.idx %in%  dimnames(trainx.covar)[[2]]],drop=F]
            }else{
		       sub.covar<-NULL
	    }
	    if(any(substr(imp.idx,1,1)=="X")){
	        x.tmp<-as.numeric(substr(imp.idx[substr(imp.idx,1,1)=="X"],2,nchar(imp.idx)))
		x.tmp<-c((2*x.tmp)-1,2*x.tmp)
		x.tmp<-x.tmp[order(x.tmp)]
		sub.xchrom<-trainx.xchrom[,x.tmp,drop=F]
	    }else{
		sub.xchrom<-NULL
	    }

	   
            sub.rf <- snpRF(x.autosome=sub.autosome[idx != i, , drop=FALSE],
	       	  	    x.xchrom=sub.xchrom[idx != i, , drop=FALSE],
	       	  	    x.covar=sub.covar[idx != i, , drop=FALSE],
                            y=trainy[idx != i],
                            xtest.autosome=sub.autosome[idx == i, , drop=FALSE],
                            xtest.xchrom=sub.xchrom[idx == i, , drop=FALSE],
                            xtest.covar=sub.covar[idx == i, , drop=FALSE],
                            ytest=trainy[idx == i],
                            mtry=mtry(n.var[j]), importance=recursive, ...)
			    

            cv.pred[[j]][idx == i] <- sub.rf$test$predicted
            ## For recursive selection, use importance measures from the sub-model.
            if (recursive) {
                impvar <-
                    dimnames(all.rf$importance)[[1]][order(sub.rf$importance[,1], decreasing=TRUE)]
      }
      NULL
    }
    NULL
  }
  ## cat("\n")
  error.cv <- sapply(cv.pred, function(x) mean(trainy != x))
  
  names(error.cv) <- names(cv.pred) <- n.var
  list(n.var=n.var, error.cv=error.cv, predicted=cv.pred)
}
