
cv.glmgraph <- function(X, Y, L,..., type.measure=c("mse","mae","deviance","auc"), nfolds=5,trace=TRUE) {
  
  type.measure <- match.arg(type.measure)
    
  obj <- glmgraph(X=X, Y=Y, L=L, ...)
  if(is.null(obj)) stop("No solution for whole data fitting, please adapt the tuning parameter.")
  
  n <- obj$n
  p <- obj$p
  lambda1 <- obj$lambda1
  nlambda1 <- obj$nlambda1
  lambda2 <- obj$lambda2
  nlambda2 <- obj$nlambda2
  
  crit <- array(NA, c(n,nlambda1, nlambda2), dimnames=list(y=1:n, lambda1=1:nlambda1,lambda2=1:nlambda2))

  if (obj$family=="gaussian") {
    cv.ind <- ceiling(sample(1:n)/n*nfolds)
  } else if (obj$family=="binomial") {
  	 if (min(table(Y)) < nfolds) stop("nfolds is larger than the smaller of 0/1 in the data; decrease nfolds")
  	 if( (n/nfolds <10) && type.measure=="auc"){
    	warning("Too few (< 10) observations per fold for type.measure='auc' in cv.glmgraph; changed to type.measure='deviance'. Alternatively, use smaller value for nfolds",call.=FALSE)
    	type.measure="deviance"
  	 }
 	 if( (n/nfolds <3)){
    	warning("Option grouped=FALSE enforced in cv.glmgraph, since < 3 observations per fold",call.=FALSE)
  	 }
     ind1 <- which(Y==1)
     ind0 <- which(Y==0)
     n1 <- length(ind1)
     n0 <- length(ind0)
     cv.ind1 <- ceiling(sample(1:n1)/n1*nfolds)
     cv.ind0 <- ceiling(sample(1:n0)/n0*nfolds)
     cv.ind <- numeric(n)
     cv.ind[Y==1] <- cv.ind1
     cv.ind[Y==0] <- cv.ind0
  }  
  
  cv.args <- list()
  cv.args$L <- L
  cv.args$family <- obj$family
  cv.args$penalty <- obj$penalty
  cv.args$gamma <- obj$gamma
  cv.args$lambda1.min.ratio <- obj$lambda1.min.ratio
  cv.args$lambda1 <- obj$lambda1
  cv.args$lambda2 <- obj$lambda2
  cv.args$eps <- obj$eps
  cv.args$max.iter <- obj$max.iter
  cv.args$dfmax <- obj$dfmax
  cv.args$penalty.factor <- obj$penalty.factor
  cv.args$standardize <- obj$standardize
  cv.args$mcpapproach <- obj$mcpapproach
  cv.args$warn <- obj$warn

  
  cv.lambda1 <- list()
  cv.lambda2 <- double(nfolds)  
  
  for (ifold in 1:nfolds) {
    if(trace) cat("Starting CV fold #",ifold,sep="","\n")
    whichi <- cv.ind==ifold
    cv.args$X <- X[!whichi, , drop=FALSE]
   	cv.args$Y <- Y[!whichi]
	obj.i <- do.call(glmgraph,cv.args)

	if(is.null(obj.i)) stop("No solution within nfold, please reset your folds again!")
    X2 <- X[whichi, , drop=FALSE]
    Y2 <- Y[whichi]
	cv.lambda1[[ifold]] <- sapply(obj.i$lambda1s,length)
	cv.lambda2[ifold] <- obj.i$nlambda2
	crit <- setcrit(obj.i,X2,Y2,crit,whichi,family=cv.args$family,type.measure=type.measure)
  }
  
  
  critcv <- setcritcv(crit,cv.ind,cv.lambda1,cv.lambda2,nlambda1,nlambda2,nfolds,Y,type.measure)
       

  ## Return
  cvm <- cvsd <- list()
  cvmat <- data.frame(matrix(rep(0,5*nlambda2),nrow=nlambda2))
  colnames(cvmat) <- c("lambda2","lambda1.min","cvmin","semin","lambda1.1se")

  for(i in 1:critcv$nlambda2.nfolds) {
  	ind <- critcv$Nmat[i,] >=nfolds
   	if(all(ind==FALSE)) {
  			warning("No criteria value for all folds for some lambda2, may need to adjust tuning parameters.")
  			next;		
  	}
	cvm[[i]] <- apply(critcv$crit[,ind,i,drop=FALSE],2,mean,na.rm=TRUE)
	cvsd[[i]] <- sqrt( apply(scale(critcv$crit[,ind,i],cvm[[i]],FALSE)^2,2,mean,na.rm=TRUE)/critcv$Nmat[i,ind])
   	if(type.measure=="auc") lmin <- getmin(lambda1[seq(cvm[[i]])],-cvm[[i]],cvsd[[i]])
   	else lmin <- getmin(lambda1[seq(cvm[[i]])],cvm[[i]],cvsd[[i]])
  	cvmat[i,] <- c(lambda2[i],lmin$lambda.min, lmin$cvmin,lmin$semin,lmin$lambda.1se)
  } 
  

  cvm <- cvm[!sapply(cvm, is.null)] 
  cvsd <- cvsd[!sapply(cvsd, is.null)] 
  cvmat <- cvmat[!apply(cvmat,1,sum)==0,]
    
  id.min <- which.min(cvmat$cvmin)
  id.1se <- which.min(cvmat$semin)
  cvmin <- min(cvmat$cvmin)
  cv.1se <-  min(cvmat$semin)
  
  lambda2.min <- cvmat$lambda2[id.min]
  lambda1.min<-  cvmat$lambda1.min[id.min]
  lambda2.1se <- cvmat$lambda2[id.1se]
  lambda1.1se <-  cvmat$lambda1.1se[id.1se]
    
  cv.args$X <- X
  cv.args$Y <- Y
  cv.args$lambda1 <- lambda1.min
  cv.args$lambda2 <- lambda2.min
  cv.args$dfmax	<-	obj$dfmax
  
  obj.min <- do.call(glmgraph,cv.args)
  
  if(is.null(obj.min)){
    flag.break <- FALSE
  	stepsize <- 0
  	if ( (id.min-stepsize)>=1 ){
  		while( (id.min-stepsize)>=1 ){
  				if(is.element(lambda1.min,lambda1)) ind.lambda1 <- which(lambda1==lambda1.min)+1
  				else ind.lambda1 <- 1
  		 		if(ind.lambda1<=nlambda1){ #decrease lambda1
  		 			while(ind.lambda1<=nlambda1){
  		 				cv.args$lambda1 <- lambda1[ind.lambda1]
  		 				obj.min <- do.call(glmgraph,cv.args)
  		 				if(!is.null(obj.min)){
  		 					lambda1.min <- cv.args$lambda1
  		 					flag.break <- TRUE
  		 					break;
  		 				}
  		 				ind.lambda1 <- ind.lambda1+1
  		 			}
  		 		}else{ #decrease lambda2
  					if(is.element(lambda2.min,lambda2)) ind.lambda2 <- which(lambda2==lambda2.min)-1
  					else ind.lambda2 <- nlambda2
  					while(ind.lambda2>0){
  		 				cv.args$lambda2 <- lambda2[ind.lambda2]
  		 				obj.min <- do.call(glmgraph,cv.args)
  		 				if(!is.null(obj.min)){
  		 					lambda2.min <- cv.args$lambda2
  		 					flag.break <- TRUE
  		 					break;
  		 				}
  		 				ind.lambda2 <- ind.lambda2-1
  		 			}
  		 		}
  		 		if(flag.break) break;
  		 		stepsize <- stepsize+1
  		 		cv.args$lambda1 <- cvmat$lambda2[id.min-stepsize] #decrease both lambda1 and lambda2
  		 		cv.args$lambda2 <- cvmat$lambda1[id.min-stepsize]
  		 		lambda1.min <- cv.args$lambda1
  		 		lambda2.min <- cv.args$lambda2
  		 		obj.min <- do.call(glmgraph,cv.args)
  		 		if(!is.null(obj.min)) break;
 		}
  	}
  } 
  if(is.null(obj.min)) warning("No solution for all tuning parameters. Perhaps need to increase max.iter or adapt tuning parameters.")


  cv.args$lambda1 <- lambda1.1se
  cv.args$lambda2 <- lambda2.1se
  obj.1se <- do.call(glmgraph,cv.args)
  
 if(is.null(obj.1se)){
    flag.break <- FALSE
  	stepsize <- 0
  	if ( (id.1se-stepsize)>=1 ){
  		while( (id.1se-stepsize)>=1 ){
  				if(is.element(lambda1.1se,lambda1)) ind.lambda1 <- which(lambda1==lambda1.1se)+1
  				else ind.lambda1 <- 1
  		 		if(ind.lambda1<=nlambda1){ #decrease lambda1
  		 			while(ind.lambda1<=nlambda1){
  		 				cv.args$lambda1 <- lambda1[ind.lambda1]
  		 				obj.1se <- do.call(glmgraph,cv.args)
  		 				if(!is.null(obj.1se)){
  		 					lambda1.1se <- cv.args$lambda1
  		 					flag.break <- TRUE
  		 					break;
  		 				}
  		 				ind.lambda1 <- ind.lambda1+1
  		 			}
  		 		}else{ #decrease lambda2
  					if(is.element(lambda2.1se,lambda2)) ind.lambda2 <- which(lambda2==lambda2.1se)-1
  					else ind.lambda2 <- nlambda2
  					while(ind.lambda2>0){
  		 				cv.args$lambda2 <- lambda2[ind.lambda2]
  		 				obj.1se <- do.call(glmgraph,cv.args)
  		 				if(!is.null(obj.1se)){
  		 					lambda2.1se <- cv.args$lambda2
  		 					flag.break <- TRUE
  		 					break;
  		 				}
  		 				ind.lambda2 <- ind.lambda2-1
  		 			}
  		 		}
  		 		if(flag.break) break;
  		 		stepsize <- stepsize+1
  		 		cv.args$lambda1 <- cvmat$lambda2[id.1se-stepsize] #decrease both lambda1 and lambda2
  		 		cv.args$lambda2 <- cvmat$lambda1[id.1se-stepsize]
  		 		lambda1.1se <- cv.args$lambda1
  		 		lambda2.1se <- cv.args$lambda2
  		 		obj.1se <- do.call(glmgraph,cv.args)
  		 		if(!is.null(obj.1se)) break;
 		}
  	}
  } 
  if(is.null(obj.1se)) warning("No solution for all tuning parameters. Perhaps need to increase max.iter or adapt tuning parameters.")

  if(type.measure=="auc") {
  	cvmat$cvmin <- -cvmat$cvmin
  	cvmat$semin <- -cvmat$semin
  }
  obj.cv <- list(obj=obj,obj.min=obj.min,obj.1se=obj.1se,cvmat=cvmat,cvmin=cvmin,cv.1se=cv.1se,cvm=cvm,cvsd=cvsd,lambda1.min=lambda1.min,lambda2.min=lambda2.min,lambda1.1se=lambda1.1se,lambda2.1se=lambda2.1se,
  beta.min=obj.min$betas[[1]],beta.1se=obj.1se$betas[[1]],type.measure=type.measure)
    
  structure(obj.cv, class="cv.glmgraph")
}







