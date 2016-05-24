"plot.OptimDes" <-function(x,xscale="t",l.type=1:5,
      l.col=c("blue","green","purple","red","dark red"),CMadj=F,...)
{
  if(!identical(xscale,"t") && !identical(xscale,"n"))
    stop("unknown type value")

  num.arm<-as.numeric(x$design["num.arm"])
  num.stage<-as.numeric(x$design["num.stage"])
   
  if(isTRUE(CMadj) && !identical(num.arm,1))    
    stop("The CM adjustment is for single-arm trials only.")
  if(isTRUE(CMadj) && identical(num.stage,3))  
    stop("The CM adjustment is for two-stage trials only.")  
  

  Stime<-x$single.stageTime
  xx<-x$test["x"]
  if(!CMadj){
     n0 <- Stime["n0"]
     da <- Stime["DA"]
     sl <- Stime["SL"]
  }else{
     n0 <- Stime["n0E"]
     da <- Stime["DAE"]
     sl <- Stime["SLE"]
  }

  if (CMadj){
        CMadjfac<-Stime["n0E"]/Stime["n0"]
  }else CMadjfac<-1

  x$all.info <- x$all.info[!is.na(x$all.info[,2]),]
  if(identical(num.stage,2)){

  n.res <- CMadjfac*x$all.info[,"n"]
  if(isTRUE(max(n.res)>floor(sum(unlist(x$accrual["m.init"]))))) 
          warning("Projected patient sample size is below the minimum requirement for normal approximation adjustment")
  
  t1.res <- CMadjfac*x$all.info[,"t1"]
  MDA.res <- CMadjfac*x$all.info[,"MDA"]
  EDA.res <- CMadjfac*x$all.info[,"EDA"]
  ES.res <- CMadjfac*x$all.info[,"ES"]
  eda<- CMadjfac*x$result[1]
  es <- CMadjfac*x$result[3]
  n<-ceiling(CMadjfac*x$n[2])

  ETSL.res<-CMadjfac*(x$all.info[,"t1"]+
            (pnorm(x$all.info[,"C1U"])-pnorm(x$all.info[,"C1L"]))*(x$all.info[,"MTSL"]-
             xx-x$all.info[,"t1"])) 
  ETSL.res<- ETSL.res + (pnorm(x$all.info[,"C1U"])-pnorm(x$all.info[,"C1L"]))*xx
  MTSL.res<-CMadjfac*(x$all.info[,"MTSL"]-xx) + xx

  t1<-x$stageTime["t1"]
  C1L<-x$boundary["C1L"]
  C1U<-x$boundary["C1U"]
  mtsl<-x$stageTime["MTSL"]
  etsl<-CMadjfac*(t1+(pnorm(C1U)-pnorm(C1L))*(mtsl-xx-t1)) + (pnorm(C1U)-pnorm(C1L))*xx
  mtsl<-CMadjfac*(mtsl-xx) + xx
  }else{  ## three-stage

  n.res <- CMadjfac*x$all.info[,"n"]
  t1.res <- CMadjfac*x$all.info[,"t1"]
  t2.res <- CMadjfac*x$all.info[,"t2"]
  MDA.res <- CMadjfac*x$all.info[,"MDA"]
  EDA.res <- CMadjfac*x$all.info[,"EDA"]
  ES.res <- CMadjfac*x$all.info[,"ES"]
  ETSL.res<-CMadjfac*x$all.info[,"ETSL"]
  MTSL.res<-CMadjfac*(x$all.info[,"MTSL"]-xx) + xx  
  eda<- CMadjfac*x$result[1]
  es <- CMadjfac*x$result[3]
  n<-ceiling(CMadjfac*x$n[2])
  mtsl<-x$stageTime["MTSL"]
  etsl<-CMadjfac*x$result[2]
  mtsl<-CMadjfac*(mtsl-xx) + xx  
  
  }

  if(identical(xscale,"t"))
  {
    plot(MTSL.res/sl,ETSL.res/sl,
    xlab=paste("Ratio of single-stage study length(","sl=",round(sl,2),")",sep=" "),
    ylab="Ratio of single-stage study value",main=
    paste("target=",x$target,sep=" "),type="n",
    ylim=c(min(t1.res/sl),1.1*max(EDA.res/da)),...)
    lines(MTSL.res/sl,ETSL.res/sl,lty=l.type[1],col=l.col[1])
    lines(MTSL.res/sl,EDA.res/da,lty=l.type[2],col=l.col[2])
    lines(MTSL.res/sl,ES.res/n0,lty=l.type[3],col=l.col[3])
    lines(MTSL.res/sl,t1.res/sl,lty=l.type[4],col=l.col[4])
    if(identical(num.stage,3)) lines(MTSL.res/sl,t2.res/sl,lty=l.type[5],col=l.col[5])
    if(identical(x$target,"EDA"))
    {  
      points(mtsl/sl,eda/da,pch=18)
      text(mtsl/sl,eda/da+0.05,labels=paste("MTSL=",round(mtsl,2),"\nEDA=",round(eda,2),sep=" "),cex=0.8)
      legy<-min(etsl/sl)-0.1 
    }  
    else
    { 
      if(identical(x$target,"ETSL")){ 
      points(mtsl/sl,etsl/sl,pch=18)
      text(mtsl/sl,etsl/sl-0.05,labels=paste("MTSL=",round(mtsl,2),"\nETSL=",round(etsl,2),sep=" "),cex=0.8)  
      legy<- 1.1*max(EDA.res/da) 
      }
      if(identical(x$target,"ES")){
      points(mtsl/sl,es/n0,pch=18)
      text(mtsl/sl,es/n0+0.05,labels=paste("MTSL=",round(mtsl,2),"\nES=",round(es,2),sep=" "),cex=0.8)    
      legy<-min(etsl/sl)-0.1 
      }
    }  
    if(identical(num.stage,2)) {
    legend(max(MTSL.res/sl)/2+0.5,legy,legend=c("ETSL","EDA","ES","t1"),
           lty=l.type[1:4],col=l.col[1:4],cex=0.8,bty="n")
    }else if(identical(num.stage,3)){
    legend(max(MTSL.res/sl)/2+0.5,legy,legend=c("ETSL","EDA","ES","t1","t2"),
           lty=l.type[1:5],col=l.col[1:5],cex=0.8,bty="n")
    }
  }
  else
  {
    plot(n.res/n0,ETSL.res/sl,xlab=paste("Ratio of single-stage sample size (","n0=",n0,")",sep=" "),ylab="Ratio of single-stage study value",main=paste("target=",x$target,sep=" "),type="n",ylim=c(min(t1.res/sl),1.1*max(EDA.res/da)),...)    
    lines(n.res/n0,ETSL.res/sl,lty=l.type[1],col=l.col[1])
    lines(n.res/n0,EDA.res/da,lty=l.type[2],col=l.col[2])
    lines(n.res/n0,ES.res/n0,lty=l.type[3],col=l.col[3])
    lines(n.res/n0,t1.res/sl,lty=l.type[4],col=l.col[4])
    if(identical(num.stage,3))  lines(n.res/n0,t2.res/sl,lty=l.type[5],col=l.col[5])
    if(identical(x$target,"EDA"))
    {  
      points(n/n0,eda/da,pch=18)
      text(n/n0,eda/da+0.05,labels=paste("n=",n,"\nEDA=",round(eda,2),sep=" "),cex=0.8)   
      legy<-min(etsl/sl)-0.1   
    }  
    else
    { 
      if(identical(x$target,"ETSL")){
      points(n/n0,etsl/sl,pch=18)
      text(n/n0,etsl/sl-0.05,labels=paste("n=",n,"\nETSL=",round(etsl,2),sep=" "),cex=0.8)    
      legy<- 1.1*max(EDA.res/da)
      }
      if(identical(x$target,"ES")){
      points(n/n0,es/n0,pch=18)
      text(n/n0,es/n0+0.05,labels=paste("MTSL=",round(mtsl,2),"\nES=",round(es,2),sep=" "),cex=0.8)    
      legy<-min(etsl/sl)-0.1 
      }      
    }  

    if(identical(num.stage,2)) {
    legend(max(n.res/n0)/2+0.5,legy,legend=c("ETSL","EDA","ES","t1"),
           lty=l.type[1:4],col=l.col[1:4],cex=0.8,bty="n")
    }else if(identical(num.stage,3)){
    legend(max(n.res/n0)/2+0.5,legy,legend=c("ETSL","EDA","ES","t1","t2"),
           lty=l.type[1:5],col=l.col[1:5],cex=0.8,bty="n")
    }
  }
  
}
