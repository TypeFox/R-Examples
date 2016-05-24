plot.gamsel <-
  function(x,newx, index,which=1:p,rugplot=TRUE,ylims,...){
    gamsel.out=x
    x=newx
    if(missing(index))index=length(gamsel.out$lambdas)
    pJ=dim(gamsel.out$alpha)
    p=pJ[1]
    maxJ=pJ[2]
    nonlin <- getActive(gamsel.out, index,type="nonlinear")[[1]]
    active_linear <- getActive(gamsel.out, index,type="linear")[[1]]
    lin=setdiff(active_linear,nonlin)
    zero=setdiff(1:p,union(nonlin,lin))
    if(which[1]=="nonzero")which=seq(p)[-zero]
    if(which[1]=="nonlinear")which=nonlin
    if(is.null(which)){warning("Nothing to plot with this choice of argument which");return()}
    colvals=rep("blue",p)
    colvals[lin]="green"
    colvals[nonlin]="red"
    termmat=drop(predict.gamsel(gamsel.out,x,index=index,type="terms"))
    lambda=gamsel.out$lambda[index]
    if(missing(ylims))ylims=range(termmat[,which])
    for(ell in which){
      o=order(x[,ell])
      plot(x[o,ell], termmat[o,ell], type='l', ylab=paste("f(v",ell,")",sep=""),xlab=paste("v",ell,sep=""),ylim=ylims,col=colvals[ell],lwd=2,...)
      if(rugplot)rug(x[,ell])
    }
  }
