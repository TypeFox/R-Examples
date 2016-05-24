check.convergence <- function(nmax,M,genXn,argsXn=NULL,mode="p",epsilon=0.05,r=2,nb.sp=10,density=FALSE,densfunc=dnorm,probfunc=pnorm,tinf=-3,tsup=3,plotfunc=plot,...) {


  data <- generate(randomgen=genXn,nmax,M,argsgen=argsXn)$data


  if (mode=="p" | mode =="as") {

    critp <- criterion(data=data,epsilon=epsilon,mode="p")$crit
    critas <- criterion(data=data,epsilon=epsilon,mode="as")$crit

    tt <- p.as.plot(data,critp,critas,epsilon,nb.sp,mode=mode)
    
  }


  if (mode=="L") {

    law.plot3d(data,probfunc,tinf,tsup)
    tt <- law.plot2d(data,density=density,densfunc=densfunc,probfunc=probfunc,tinf,tsup)

  }


  
  if (mode=="r") {

    critr <- criterion(data=data,epsilon=epsilon,mode="r",r=r)$crit

    visualize.crit(critr,plotfunc=plotfunc,main=paste("Convergence in r-th mean?",sep=""),...)
    
    mtext(expression(hat(e)[n~bold(',')~'r']),side=2,line=2,las=2)
    
  }

  if (mode != "r")  return(tt)
  
}
