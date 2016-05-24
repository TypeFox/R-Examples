diag_lin <- function(t,ic,z,n.total,Nit,n.sim,n.plot,seed,minor_included){
  
  ###order data by time
  tmp <- cbind(t,ic,z)
  s.tmp <- tmp[order(tmp[,1]),]
  times <- s.tmp[,1]
  cause <- s.tmp[,2]
  desX <- as.matrix(s.tmp[,-c(1,2)])
  var.name <- colnames(z) 
  
  ###run crr() for initial value
  out.crr <- crr(times,cause,desX)

  ###fetch the results from C++
  if(out.crr$converged==T){
    N <- dim(desX)[1]
    M <- dim(desX)[2]
    if(n.plot>n.sim){
      n.plot <- n.sim
      warning("n.plot has been forced to equal to n.sim")
    }
    NJP <- 0
    TJP <- vector("numeric",length=N+1)
    betaS <- out.crr$coef
    dlamb0 <- vector("numeric",length=N+1)
    Wbeta <- vector("numeric",length=N*M)
    Ut <- vector("numeric",length=(N+1)*M)
    SI <- vector("numeric",length=M*M)
    beta.var <- vector("numeric",length=M*M) 
    
    B <- vector("numeric",length=(n.plot+1)*N*(M+1))
    uniX <- vector("numeric",length=N*(M+1))
    TC <- vector("integer",length=M+1)
    pval <- vector("numeric",length=M+1) 
    mav <- vector("numeric",length=M+1) 
    
    if(!is.null(seed)) set.seed(seed)
    
    #dyn.load("est_tmp.dll")
    out<- .C("diag_lin",
             as.double(times),as.integer(cause),as.double(desX),
             as.integer(N),as.integer(M),as.integer(Nit),
             NJP=as.integer(NJP),TJP=as.double(TJP),beta=as.double(betaS),
             beta.var=as.double(beta.var),dlamb0=as.double(dlamb0),
             uniX=as.double(uniX),TC=as.integer(TC),B=as.double(B),mav=as.double(mav),
             as.integer(n.sim),as.integer(n.plot),pval=as.double(pval),
             as.integer(minor_included))
    #dyn.unload("est_tmp.dll")
    
    beta.var <- matrix(out$beta.var,nrow=M)
    beta.se <- sqrt(diag(beta.var))
    B <- matrix(out$B,nrow=N,ncol=(M+1)*(n.plot+1))
    uniX <- matrix(out$uniX,nrow=N,ncol=M+1)
    ct <- out$NJP+1
    
    lst <- list(test="linearity and link function", varname=var.name,
                used=N,NJP=out$NJP,TJP=out$TJP[1:ct],n.plot=n.plot,
                beta=out$beta,beta.se=beta.se,dlamb0=out$dlamb0[1:ct],
                pval=out$pval,mav=out$mav,B=B,uniX=uniX,TC=out$TC)
    class(lst)<-"diaglin"
    return(lst)
  }else{
    print("crr not converged")
    return(NULL)
  }
}

print.diaglin <- function(x, ...) {
  m <- length(x$pval)
  cat("\n***P-values for model diagnostics***\n\n")
  cat("Linear functional form for a single covariate:\n")
  p.func <- x$pval[-m]
  names(p.func) <- x$varname[-m]
  p.link <- x$pval[m]
  names(p.link) <- NULL
  print(signif(p.func, 4), ...)
  cat("\nLink function:\n")
  print(signif(p.link, 4), ...)
  invisible()
}


plot.diaglin <- function(x,col=c("red","black"),lty=c(1,2),lwd=c(2,1),
                         txt.pos=c(0.85,0.1),lgd.pos="topright",xlim=NULL,ylim=NULL,
                         select=NULL,...) {
  if(is.null(x)) {
    stop("No object to plot!")
    return
  }
  if(class(x)!="diaglin") {
    stop("Wrong object class!")
  }
  
  m <- length(x$pval)
  if(!all(is.element(select,c(1:m,NULL)))) stop("one or more elements in select are invalid")
  if(is.null(select)) var <- 1:m
  else var <- select
  
  for(k in var){
    colind <- seq(k,m*x$n.plot+k,m)
    xx <- x$uniX[1:(x$TC[k]),k]
    if(x$TC[k]<=3&&k<m) warning(paste("plots may not reflect the validaty of linearity for categorical covariate:",x$varname[k]))
    if(x$TC[k]<=3&&k==m) warning(paste("plots may not reflect the validaty of link function."))

    if(is.null(xlim)) myxlim <- range(xx)
    else myxlim <- xlim
    if(is.null(ylim)) myylim <- 1.3*range(x$B[,colind])
    else myylim <- ylim
    
    mymain <- ifelse(k<m,x$varname[k],"link function")
    xlab <- ifelse(k<m,x$varname[k],expression(paste(Z^T,beta)))
    plot(myxlim,myylim,type="n",main=mymain,
         ylab="Cumulative sum of residuals",xlab=xlab,cex.lab=1.2,...)
    for(i in colind){
      mycol <- ifelse(i==k,col[1],col[2])
      mylty <- ifelse(i==k,lty[1],lty[2])
      mylwd <- ifelse(i==k,lwd[1],lwd[2])
      y <- x$B[1:(x$TC[k]),i]
      points(xx,y,type="s",xlim=xlim,ylim=ylim,col=mycol,lty=mylty,lwd=mylwd)
    }
    legend(lgd.pos,c("Observed","Simulated"),col=col,lwd=lwd,lty=lty,bty="n")
    xpos <- txt.pos[1]*(myxlim[2]-myxlim[1])+myxlim[1]
    text(xpos,txt.pos[2]*(myylim[2]-myylim[1])+myylim[1],paste("p-value=",x$pval[k]),font=1,cex=1)
  }
  invisible()
}
