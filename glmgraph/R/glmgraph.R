############## glmgraph #######################

glmgraph <- function(X, Y, L, family=c("gaussian","binomial"), penalty=c("MCP","lasso") ,mcpapproach=c("mmcd","adaptive","original"),
gamma=8, lambda1,nlambda1=100, lambda2=c(0,1e-2 * 2^(0:7)),
eps=1e-3,max.iter=2000,dfmax=round(ncol(X)/2),penalty.factor=rep(1,ncol(X)),standardize=TRUE,warn=FALSE,...)
{

	family <- match.arg(family)
  	penalty <- match.arg(penalty)
  	mcpapproach <- match.arg(mcpapproach)
  	if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty")
  	if (length(penalty.factor)!=ncol(X)) stop("penalty.factor does not match up with X")
  	
  	diagL <- diag(L)
  	diag(L) <- 0
  	L <- -L
  	
  	std <- XX <- X.center <- X.scale <- bt <- NULL
  	if(standardize){
  		std <- standardizeX(X)
  		XX <- std[[1]]
  		X.center <- std[[2]]
  		X.scale <- std[[3]] 
  	}else{
  		XX <- X
  		X.scale <- apply(X,2,sd)
  	}
  	
  	nz <- which(X.scale > 1e-6)  	
  	
  	if (length(nz)==0)  stop("please specify penalty.cols again")
  	if (length(nz) != ncol(XX) ){
  		 XX <- XX[ ,nz, drop=FALSE]
  	}
  	  	
  	YY <- Y <- as.numeric(Y)  
  		  	
  	if(family=="gaussian")  YY <- scale(Y)
  	  	
  	n <- length(Y)
  	p <- ncol(XX)
  	lambda1.min.ratio=ifelse(n>p,1e-4,1e-2)

  	if (missing(lambda1)) {
    	lambda1 <- setupLambda(XX, YY, family,lambda1.min.ratio, nlambda1, penalty.factor[nz])
    	lambda1 <- as.numeric(formatC(lambda1,digits=6, format="f"))
  	}
  	
    nlambda1 <- length(lambda1)
  	nlambda2 <- length(lambda2)
  	
  	## Fit
  	if (family=="gaussian") {
		res <- .Call("cdfit_gaussian", XX, YY, L[nz,nz], diagL[nz], penalty, lambda1, lambda2, gamma, penalty.factor[nz], eps, max.iter, dfmax,standardize,warn)  		
  	}else if (family=="binomial") {
  		res <- .Call("cdfit_binomial", XX, YY, L[nz,nz], diagL[nz], penalty, lambda1, lambda2, gamma, penalty.factor[nz], eps, max.iter, dfmax, mcpapproach, warn)  		
  	}
  	
  	b	<- res[[1]]
    iter <- res[[2]]
       
    lambda2idx <- NULL
    for(i in seq(nlambda2)){
    	if(sum(is.na(iter[[i]]))<nlambda1) {
    		lambda2idx <- c(lambda2idx,i)
    	}else{
    		if (warn) warning(paste( "Algorithm failed to converge for some lambda2: ",lambda2[i],"\n"))
    	}
    }
    if(is.null(lambda2idx)) {
    	warning(paste( "No solution under current lambda1: ",lambda1," and lambda2 combination: ",lambda2,"\n"))
    	return(NULL)
    }
    
    b <- b[lambda2idx]
    iter <- iter[lambda2idx]
    
    nlambda2 <- length(lambda2idx)
    lambda2 <- lambda2[lambda2idx]
    
    
    ## Eliminate saturated lambda1 values for each lambda2, if any
  	betas <- lambda1s <-loglik <- df <- list() 
   
  	## unstandardizeX
  	for(i in 1:nlambda2){
  		ind <- !is.na(iter[[i]])
     	if (any(is.na(drop(iter[[i]][ind]))) && warn) warning(paste( "Algorithm failed to converge for some values of lambda1 at lambda2: ",lambda2[i]))
  		if(length(iter[[i]][ind])>0){
  			if(standardize) bt <- unstandardizeX(b[[i]][,ind, drop=FALSE], X.center[nz], X.scale[nz],Y,family=family)
  			else bt <- b[[i]][,ind, drop=FALSE]*sd(Y)
  			beta <- matrix(0, nrow=(ncol(X)+1), ncol=length(lambda1[ind]))
  			beta[1,] <- bt[1,]
  			beta[nz+1,] <- bt[-1,]
  			varnames <- colnames(X)
  			if (is.null(colnames(X))) varnames <- paste("V",1:ncol(X),sep="")
 			varnames <- c("(Intercept)", varnames)
  			dimnames(beta) <- list(varnames, round(lambda1[ind],digits=4))  
  			betas[[i]] <- beta
  			lambda1s[[i]] <- lambda1[ind]
  		}  		
  		##add logLik,df
  		if(family=="gaussian"){  			
  			eta <- sweep(X %*% beta[-1,,drop=FALSE], 2, beta[1,,drop=FALSE], "+")
  			p1 <- sweep(eta,1,Y,"-")
  			sig2 <- apply(p1^2,2,var)
  			p1 <- apply(p1^2,2,sum)  		
  			loglik[[i]] <- -n/2*log(2*pi*sig2)-p1/2/sig2			
  		}else if(family=="binomial"){
  			eta <- sweep(X %*% beta[-1,,drop=FALSE], 2, beta[1,,drop=FALSE], "+")
  			p1 <- sweep(eta,1,Y,"*")
  			p1 <- apply(p1,2,sum)
  			p2 <- log(1+exp(eta))
  			p2 <- apply(p2,2,sum)
  			loglik[[i]] <- p1-p2
  		}
	}

	df <- lapply(betas,function(beta) apply(beta[-1,,drop=FALSE]!=0,2,sum))
	names(betas) <- names(df) <- names(lambda1s) <- names(loglik) <- lambda2
	
  ## Output
  obj <- list(betas = betas,
                        lambda1s = lambda1s,
                        p=p,
                        n = n,
                        family=family,
						penalty=penalty,
						mcpapproach=mcpapproach,
						gamma=gamma,
						lambda1.min.ratio=lambda1.min.ratio,
						lambda1=lambda1,
						nlambda1=nlambda1,
						lambda2=lambda2,
						nlambda2=nlambda2,
						eps=eps,
						max.iter=max.iter,
						dfmax=dfmax,
						standardize=standardize,
						penalty.factor=penalty.factor,
						loglik=loglik,
						df=df			
				)
 
 
  structure(obj,class="glmgraph")
}




















