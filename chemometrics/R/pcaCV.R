"pcaCV" <-
function(X,amax,center=TRUE,scale=TRUE,repl=50,
  segments=4,segment.type = c("random", "consecutive", "interleaved"),length.seg,
  trace = FALSE,plot.opt=TRUE, ...)
{
# Repeated Double Cross Validation (PF, November 6, 2007) for obtaining the
# optimal number of PCA components
#    see also function "mvr" and "mvrCv" from library(pls)
#
# INPUT: 
# X ... X data
# amax ... maximum number of components to use, e.g. the number of X-variables
# center ... centring of X (FALSE or TRUE)
# scale ... scaling of X (FALSE or TRUE)
# repl ... number of replication for Cross Validation
# segments ... number of segments used for splitting into training and test data
# segment.type ... "random", "consecutive", "interleaved" splitting into training and test data
# length.seg ... number of parts for training and test data, overwrites segments
# trace ... if TRUE intermediate results are reported
# plot.opt ... make plot with optimum choice of components
# ... additional plotting arguments
#
# OUTPUT:
# MSEP ... MSEPs for all numbers of components
# ExplVar ... explained variance (R^2) for all numbers of components


#require(pls)

    if (missing(amax)) {
        amax <- min(nrow(X) - 1, ncol(X))
    }
    else {
        if (amax < 1 || amax > min(nrow(X) - 1, ncol(X)))
            stop("Invalid number of components, amax")
    }

    ###################################################################################################
    X <- as.matrix(scale(X,center=center,scale=scale))
    amax <- min(amax, nrow(X) - max(sapply(segments, length)) - 1)

    optcomp <- matrix(NA,nrow=segments,ncol=repl)
    MSEP <- matrix(NA,nrow=repl,ncol=amax)
    dimnames(MSEP) <- list(paste("rep",1:repl), paste("PC",1:amax))
    Fit <- matrix(NA,nrow=repl,ncol=amax)
    dimnames(Fit) <- list(paste("rep",1:repl), 1:amax)
    for (i in 1:repl){

      if (missing(length.seg)) {
            segment <- cvsegments(nrow(X), k = segments, type = segment.type)
        }
      else {
            segment <- cvsegments(nrow(X), length.seg = length.seg, type = segment.type)
      }
      if (trace)
        cat(paste("Replication: ",i))

      MSEPj <- matrix(NA,nrow=segments,ncol=amax)
      Fitj <- matrix(NA,nrow=segments,ncol=amax)

      for (n.seg in 1:length(segment)) {
        if (trace)
            cat(n.seg, "")
        seg <- segment[[n.seg]]
        obsuse <- as.numeric(unlist(segment[-n.seg]))
	Xtrain <- X[obsuse,]
        obstest <- as.numeric(unlist(segment[n.seg]))
	Xtest <- X[obstest,]

	# PCA
        if (ncol(Xtrain)>nrow(Xtrain)){
		e <- eigen(Xtrain%*%t(Xtrain))
		Ttrain <- e$vectors %*% diag(sqrt(e$values))
		Ptrain <- t(Xtrain) %*% Ttrain %*% diag(1/e$values)
	}
	else {
		Xtrain_svd <- svd(Xtrain)
		Ptrain <- Xtrain_svd$v
	}

	Ttest <- Xtest %*% Ptrain
	for (j in 1:amax){
	  MSEPj[n.seg,j] <- sum((Xtest - Ttest[,1:j] %*% t(Ptrain[,1:j]))^2)
	}
	Fitj[n.seg,] <- MSEPj[n.seg,]/sum(Xtest^2)
      }

      MSEP[i,] <- apply(MSEPj,2,mean)
      Fit[i,] <- 1-apply(Fitj,2,mean)
    }

    if (plot.opt){
      boxplot(as.data.frame(Fit),ylab="Explained variance",xlab="Number of components",...)
    }

    list(ExplVar=Fit,MSEP=MSEP)
}

