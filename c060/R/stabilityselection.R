#require(glmnet)
#require(parallel)

stabpath <- function(y,x,size=0.632,steps=100,weakness=1,mc.cores=getOption("mc.cores", 2L),...){
  fit <- glmnet(x,y,...)
  if(class(fit)[1]=="multnet"|class(fit)[1]=="lognet") y <- as.factor(y)
  #if(class(fit)[1]=="lognet") y <- as.logical(y) 
  p <- ncol(x)
  #draw subsets
  subsets <- sapply(1:steps,function(v){sample(1:nrow(x),nrow(x)*size)})
  
  # parallel computing depending on OS
  # UNIX/Mac
  if (.Platform$OS.type!="windows") {
    res <- mclapply(1:steps,mc.cores=mc.cores,glmnet.subset,subsets,x,y,lambda=fit$lambda,weakness,p,...)
  } else {
    # Windows  
    cl  <- makePSOCKcluster(mc.cores)
    clusterExport(cl,c("glmnet","drop0"))
    res <- parLapply(cl, 1:steps,glmnet.subset,subsets,x,y,lambda=fit$lambda,weakness,p,...)
    stopCluster(cl)
  }
  
  #merging
  res <- res[unlist(lapply(lapply(res,dim),function(x) x[2]==dim(res[[1]])[2]))]
  x <- as.matrix(res[[1]])
  qmat <- matrix(ncol=ncol(res[[1]]),nrow=length(res))
  qmat[1,] <- colSums(as.matrix(res[[1]]))
  for(i in 2:length(res)){
    qmat[i,] <- colSums(as.matrix(res[[i]]))
    x <- x + as.matrix(res[[i]])
  }
  x <- x/length(res)
  qs <- colMeans(qmat)
  out <- list(fit=fit,x=x,qs=qs)	
  class(out) <- "stabpath" 
  return(out)
}

#internal function used by lapply 
glmnet.subset <- function(index,subsets,x,y,lambda,weakness,p,...){
  if(length(dim(y))==2|class(y)=="Surv"){
    glmnet(x[subsets[,index],],y[subsets[,index],],lambda=lambda
           ,penalty.factor= 1/runif(p,weakness,1),...)$beta!=0
  }else{
    if(is.factor(y)&length(levels(y))>2){
      temp <- glmnet(x[subsets[,index],],y[subsets[,index]],lambda=lambda
                     ,penalty.factor= 1/runif(p,weakness,1),...)[[2]]
      temp <- lapply(temp,as.matrix)
      Reduce("+",lapply(temp,function(x) x!=0))
      
    }	
    else{
      glmnet(x[subsets[,index],],y[subsets[,index]],lambda=lambda
             ,penalty.factor= 1/runif(p,weakness,1),...)$beta!=0
    }
  }	
}


#performs error control and returns estimated set of stable variables and corresponding lambda
stabsel <- function(x,error=0.05,type=c("pfer","pcer"),pi_thr=0.6){
  if(pi_thr <= 0.5 | pi_thr >= 1) stop("pi_thr needs to be > 0.5 and < 1!")
  if(class(x$fit)[1]=="multnet"){
  p <- dim(x$fit$beta[[1]])[1]
  }else{
	p <- dim(x$fit$beta)[1]
  }
  type <- match.arg(type)
  switch(type,
         "pcer"={
           if(error>=1 | error<=0)stop("pcer needs to be > 0 and < 1!")
           qv <- ceiling(sqrt(error* p * (2*pi_thr-1)*p)) },
         "pfer"={
          qv <- ceiling(sqrt(error * (2*pi_thr-1)*p)) }
         )
  if(x$qs[length(x$qs)]<=qv){ lpos <- length(x$qs)
    }else{
  lpos <- which(x$qs>qv)[1]
  }
	if(!is.na(lpos)){stable <- which(x$x[,lpos]>=pi_thr)}else{
    stable <- NA
	}
	out <- list(stable=stable,lambda=x$fit$lambda[lpos],lpos=lpos,error=error,type=type)
	return(out)
}

print.stabpath <- function(x,...){
  cat(" stabilitypath","\n",
      dim(x$x)[1],"variables","\n",
      dim(x$x)[2],"lambdas","\n")
}

#plot penalization and stability path 
plot.stabpath <- function(x,error=0.05,type=c("pfer","pcer"),pi_thr=0.6,xvar=c("lambda", "norm", "dev")
                          , col.all="black", col.sel="red",...){
  sel <- stabsel(x,error,type,pi_thr)
  if(class(x$fit)[1]=="multnet"){
    beta = as.matrix(Reduce("+",x$fit$beta))
  }else{
    beta = as.matrix(x$fit$beta)
  }  
    p <- dim(beta)[1]
    which = nonzeroCoef(beta)
    nwhich = length(which)
    switch(nwhich + 1, `0` = {
      warning("No plot produced since all coefficients zero")
      return()
    }, `1` = warning("1 or less nonzero coefficients; glmnet plot is not meaningful"))
    xvar = match.arg(xvar)
    switch(xvar, norm = {
      index = apply(abs(beta), 2, sum)
      iname = "L1 Norm"
    }, lambda = {
      index = log(x$fit$lambda)
      iname = expression(paste("log ",lambda))
    }, dev = {
      index = x$fit$dev
      iname = "Fraction Deviance Explained"
    })
  #}
  #stability path
  cols <- rep(col.all,p)
  cols[sel$stable] <- col.sel
  lwds <- rep(1,p)
  lwds[sel$stable] <- 2
  if(!class(x$fit)[1]=="multnet"){
  par(mfrow=c(2,1))
  matplot(y=t(beta), x=index
          ,type="l",col=cols,lwd=lwds,lty=1,ylab=expression(paste(hat(beta)[i]))
          ,xlab=iname,main="Penalization Path",cex.lab=1,cex.axis=1,las=1,...)
  }
  matplot(y=as.matrix(t(x$x)), x=index
          ,type="l",col=cols,lwd=lwds,lty=1,ylab=expression(paste(hat(Pi)))
          ,xlab=iname,main="Stability Path",ylim=c(0,1),cex.lab=1,cex.axis=1,las=1,...)
  abline(h=pi_thr,col="darkred",lwd=1,lty=1)
  abline(v=index[sel$lpos],col="darkred",lwd=1,lty=1)
  #text(x=20,y=0.9,paste(expression(paste(lambda)),"=",paste(round(sel[[2]],digits=3)),sep=""),cex=0.75)
  return(sel)
}

