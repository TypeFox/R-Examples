
getmin <- function(lambda,cvm,cvsd){
  cvmin <- min(cvm,na.rm=TRUE)
  idmin <- cvm<=cvmin
  lambda.min <- max(lambda[idmin],na.rm=TRUE)
  idmin <- match(lambda.min,lambda)
  semin <- (cvm+cvsd)[idmin]
  idmin <- cvm<=semin
  lambda.1se <- max(lambda[idmin],na.rm=TRUE)
  list(lambda.min=lambda.min,lambda.1se=lambda.1se,cvmin=cvmin,semin=semin)
}



auc <- function(Y,predmat){
   if(is.null(ncol(predmat)))
	  predmat <- as.matrix(predmat)
	rprobmat <- apply(predmat,2,rank)
	n1 <- sum(Y);n0 <- length(Y)-n1	
	if(n1==0 || n0==0 ) {
		warning("No two classes within fold.")
		return(rep(0,length(Y)))
	}
	R1 <- apply(rprobmat[Y==1,,drop=FALSE],2,sum)
	umat <- R1-n1*(n1+1)/2
	umat <- umat/(n1*n0) #ncol(predmat) is length of lambda1 for current lambda2
	auc <- pmax(umat,1-umat)
	return(auc)
}



setcrit <- function(obj,X,Y,crit,whichi,family=family,type.measure=type.measure){
	pred <- predict(obj, X, type="response")
	if(family=="gaussian"){
		if(type.measure=="mse" || type.measure=="deviance"){       #ncol(pred[[i]])<nlambda1
			for(i in seq(obj$nlambda2)) crit[whichi,seq(ncol(pred[[i]])),i] <- ( sweep( pred[[i]],1,Y,FUN="-") )^2
		}else if(type.measure=="mae"){
			for(i in seq(obj$nlambda2)) crit[whichi,seq(ncol(pred[[i]])),i] <- abs( sweep( pred[[i]],1,Y,FUN="-") )
		}
	}else if(family=="binomial"){
		prob_min <- 1e-5
  		prob_max <- 1-prob_min	
		if(type.measure=="mse"){
			for(i in seq(obj$nlambda2)) {
				predmat <- pmin(pmax(pred[[i]],prob_min),prob_max)
				crit[which(whichi)[Y==1],seq(ncol(predmat)),i] <- ( sweep( predmat[Y==1,],1,Y[Y==1],FUN="-") )^2
				crit[which(whichi)[Y==0],seq(ncol(predmat)),i] <- ( sweep( 1-predmat[Y==0,],1,Y[Y==0],FUN="-") )^2
			}
		}else if(type.measure=="mae"){
			for(i in seq(obj$nlambda2)) {
				predmat <- pmin(pmax(pred[[i]],prob_min),prob_max)
				crit[which(whichi)[Y==1],seq(ncol(predmat)),i] <- abs( sweep( predmat[Y==1,],1,Y[Y==1],FUN="-") )
				crit[which(whichi)[Y==0],seq(ncol(predmat)),i] <- abs( sweep( 1-predmat[Y==0,],1,Y[Y==0],FUN="-") )
			}
		}else if(type.measure=="deviance"){
			for(i in seq(obj$nlambda2)) {
				predmat <- pmin(pmax(pred[[i]],prob_min),prob_max)
				crit[which(whichi)[Y==1],seq(ncol(predmat)),i] <-  -2*log(predmat[Y==1, , drop=FALSE])
				crit[which(whichi)[Y==0],seq(ncol(predmat)),i] <-  -2*log(1-predmat[Y==0, , drop=FALSE])
			}
		}else if(type.measure=="auc"){
			for(i in seq(obj$nlambda2)) {
				predmat <- pmin(pmax(pred[[i]],prob_min),prob_max)
				crit[whichi,seq(ncol(predmat)),i] <- predmat
			}
		}
	}
	return(crit)
}



#setcritcv for others except auc
setcritcv <- function(crit,cv.ind,cv.lambda1,cv.lambda2,nlambda1,nlambda2,nfolds,Y,type.measure=type.measure){
		outmat <-  array(NA,c(nfolds,nlambda1,nlambda2))
		good  <- array(0,c(nfolds,nlambda1,nlambda2))
		nlambda2.nfolds <- max(cv.lambda2)	   # max #lambda2 for different fold
		
		crit[is.infinite(crit)] <- NA
    	for(ifold in seq(nfolds)){
    		whichi <- cv.ind==ifold
			for(i in seq(cv.lambda2[ifold])){  # different #lambda2 for different fold
		    	if(type.measure!="auc") outmat[ifold,,i]=apply(crit[whichi,,i,drop=FALSE],2,mean,na.rm=TRUE)
		    	else outmat[ifold,seq(cv.lambda1[[ifold]][i]),i] <- auc(Y[whichi],crit[whichi,seq(cv.lambda1[[ifold]][i]),i])
				good[ifold,seq(cv.lambda1[[ifold]][i]),i] <- 1
			}
		}
		outmat <- outmat[,,seq(nlambda2.nfolds),drop=FALSE] # remove lambda2 without any solution for any fold
		Nmat <- matrix(0, nlambda2.nfolds,nlambda1)
		for(i in seq( nlambda2.nfolds)){
			Nmat[i,] <- apply(good[,,i,drop=FALSE],2,sum)
		}
		return(list(crit=outmat,Nmat=Nmat,nlambda2.nfolds=nlambda2.nfolds))
}




