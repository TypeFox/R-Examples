diag_prop <- function(t,ic,z,n.total,Nit,n.sim,n.plot,seed,minor_included){

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
      warning("n.plot has force to n.sim")
    }
    NJP <- 0
    TJP <- vector("numeric",length=N+1)
    betaS <- out.crr$coef
    dlamb0 <- vector("numeric",length=N+1)
    Wbeta <- vector("numeric",length=N*M)
    Ut <- vector("numeric",length=(N+1)*M)
    SI <- vector("numeric",length=M*M)
    beta.var <- vector("numeric",length=M*M) 
    
    Bt <- vector("numeric",length=(n.plot+1)*(N+1)*M)
    pvalt <- vector("numeric",length=M+1)  
    mavt <- vector("numeric",length=M+1) 
    
    if(!is.null(seed)) set.seed(seed)

    #dyn.load("est_tmp.dll")
    out<- .C("diag_prop",
             as.double(times),as.integer(cause),as.double(desX),
             as.integer(N),as.integer(M),as.integer(Nit),
             NJP=as.integer(NJP),TJP=as.double(TJP),beta=as.double(betaS),
             beta.var=as.double(beta.var),dlamb0=as.double(dlamb0),
             Bt=as.double(Bt),mavt=as.double(mavt),
             as.integer(n.sim),as.integer(n.plot),pvalt=as.double(pvalt),
             as.integer(minor_included))
    #dyn.unload("est_tmp.dll")

    beta.var <- matrix(out$beta.var,nrow=M)
    beta.se <- sqrt(diag(beta.var))
    ct <- out$NJP+1
    Bt <- matrix(out$Bt,nrow=N+1,ncol=M*(n.plot+1))[1:ct,]  

    lst <- list(test="proportionality", varname=var.name,
                used=N,NJP=out$NJP,TJP=out$TJP[1:ct],n.plot=n.plot,
                beta=out$beta,beta.se=beta.se,dlamb0=out$dlamb0[1:ct],
                pval=out$pvalt,mav=out$mavt,B=Bt)
    class(lst)<-"diagprop"
    return(lst)
  }else{
    print("crr not converged")
    return(NULL)
  }

}


print.diagprop <- function(x, ...) {
  m <- length(x$pval)
  p.single <- x$pval[-m]
  names(p.single) <- x$varname[-m]
  p.overall <- x$pval[m]
  names(p.overall) <- NULL
  cat("\n***P-values for model diagnostics***\n\n")
  cat("Proportionality for a single covariate:\n")
  print(signif(p.single, 4), ...)
  cat("\nOverall proportionality:\n")
  print(signif(p.overall, 4), ...)
  invisible()
}


plot.diagprop <- function(x,col=c("red","black"),lty=c(1,2),lwd=c(2,1),
                          txt.pos=c(0.85,0.1),lgd.pos="topright",xlim=NULL,ylim=NULL,
                          select=NULL,...) {
  if(is.null(x)) {
    stop("No object to plot!")
  }
  if(class(x)!="diagprop") {
    stop("Wrong object class!")
  }
  
  m <- length(x$beta)
  if(!all(is.element(select,c(1:m,NULL)))) stop("one or more elements in select are invalid")
  if(is.null(select)) var <- 1:m
  else var <- select
  
  for(k in var){
    colind <- seq(k,m*x$n.plot+k,m)
    xx <- x$TJP

    if(is.null(xlim)) myxlim <- range(xx)
    else myxlim <- xlim
    if(is.null(ylim)) myylim <- 1.2*range(x$B[,colind])
    else myylim <- ylim
    
    plot(myxlim,myylim,type="n",main=x$varname[k],
         ylab="Standardized Process",xlab="Follow-up Time",cex.lab=1.2,...)
    for(i in colind){
      mycol <- ifelse(i==k,col[1],col[2])
      mylty <- ifelse(i==k,lty[1],lty[2])
      mylwd <- ifelse(i==k,lwd[1],lwd[2])
      y <- x$B[,i]
      points(xx,y,type="s",xlim=myxlim,ylim=myylim,col=mycol,lty=mylty,lwd=mylwd)
    }
    legend(lgd.pos,c("Observed","Simulated"),col=col,lwd=lwd,lty=lty,bty="n")
    xpos <- txt.pos[1]*(myxlim[2]-myxlim[1])+myxlim[1]
    text(xpos,txt.pos[2]*(myylim[2]-myylim[1])+myylim[1],paste("p-value=",x$pval[k]),font=1,cex=1)
  }
}
