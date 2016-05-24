library('fgui')# For R gui
library('tcltk')# For R gui
library('SMPracticals')# for qqexp
##############################################################
# n is the target sample size
# T is the target completion time
# P is the prior certainty (0 <= P <= 1)
# m is the sample observed to date
# tm is the time to date
# The inforamtive accrual is using a closed form 
# nmax is the upper limit on the sample size axis of the graph
# set P to zero for a non-informative prior
# set m and tm to zero if no data has been accumulated yet.


### Simulation of an example data for accrual

set.seed(0)
accrual.data<-round(rexp(300,1/(365.25*3/300)))


### Major function for accrual prediction
accrual.n<- function(n,T,P,m,tm,Tp) {
  prob=(T*P+tm)/(T*P+Tp)
  r=n*P+m
  pred.n.Tp=qnbinom(c(0.025,0.5,0.975),r,prob)+m
  dist.n.Tp=rnbinom(10000,r,prob)+m
  return(list(pred.n.Tp,dist.n.Tp))
}

accrual.T<- function(n,T,P,m,tm,np) {
  qB=qbeta(c(0.025,0.5,0.975), np-m, n*P+m)
  rB=rbeta(10000,np-m, n*P+m)
  qTOT=(T*P+tm)*qB/(1-qB)+tm
  rTOT=(T*P+tm)*rB/(1-rB)+tm
  return(list(qTOT,rTOT))
}


## function for accrual plot
accrual.n.plot <- function(n,T,P,m,tm,Tp){
  tlist=seq(tm+1,Tp,length=100)
  accrual.count <- matrix(NA,nrow=100,ncol=3)
  for (i in 1:100) {
    time=tlist[i]
    accrual.count[i,]=accrual.n(n,T,P,m,tm,Tp=time)[[1]]
                               accural.n.dist=accrual.n(n,T,P,m,tm,Tp=time)[[2]]
    }
 
  
 
## Calcualte the duration for n subjects
    
   accrual.T.n=accrual.T(n,T,P,m,tm,n)[[1]]
   lcln=accrual.count[,1]
    midn=accrual.count[,2]
    ucln=accrual.count[,3]

  layout(matrix(c(2,2,2,1),nrow=1))
  par(mar=c(4.1,0.1,2.1,0.1))
  accrual.hist <- cut(accural.n.dist,
    seq(0,max(ucln)*1.2,length=40))
  barplot(table(accrual.hist),horiz=TRUE,
    axes=FALSE,xlab=" ",ylab=" ",space=0,
    col="white",names.arg=rep(" ",39),
    )

  par(mar=c(4.1,4.1,2.1,0.1))
  plot(c(-0.5,Tp),c(0,max(ucln)*1.2),xlab="Time (Months)",
    ylab="Number of patients",type="n",xaxt = 'n')
  axis(1, at=seq(0,T, 6)) 
   
  legenda=paste("Total targeted subjects:", n)
  legendb=paste("Total finish time (months):",T)
  legendd=paste("Time to date (months) :", tm)
  legendc=paste("Subjects recruited to date:",m)
  legende=paste("Subjects in", T, "months:",round(lcln[100]),"(",round(midn[100]),",",round(ucln[100]),")" )
  legendf=paste("Time for",n, "subjects:", round(accrual.T.n[2],digits=1),"(",round(accrual.T.n[1],digits=1),",",round(accrual.T.n[3],digits=1),")" )

  legend(-0.5, max(ucln)*1.2, legend=c("Input Information:",
         legenda,legendb,legendc,legendd,"------------------------","Summary of Results:",legende,legendf))
  polygon(c(c(tm,tlist),rev(c(tm,tlist))),c(c(m,lcln),rev(c(m,ucln))),
    density=-1,col="gray",border=NA)
  lines(c(tm,tlist),c(m,midn),col="white")
  segments(0,0,tm,m)
  lines(rep(n,100),col="red")
  return(list(paste("2.5%=",round(lcln[100])), paste("50%=", round(midn[100])),paste("97.5%=",round(ucln[100]))))
}


accrual.T.plot=function(n,T,P,m,tm,np) {
  nlist=seq(m+1,np)
  accrual.time=matrix(NA,nrow=length(nlist),ncol=3)
  for (i in 1:length(nlist)){
   npred=nlist[i]
   accrual.time[i,]=accrual.T(n,T,P,m,tm,npred)[[1]]
   accural.T.dist=accrual.T(n,T,P,m,tm,np)[[2]]
  
  }

 
## Calcualte the number of subjects for T time
    accrual.n.T=accrual.n(n,T,P,m,tm,T)[[1]]
    
    lclT=accrual.time[,1]
    midT=accrual.time[,2]
    uclT=accrual.time[,3]


  layout(matrix(c(1,2,2,2)))
  par(mar=c(0.1,4.1,0.1,0.1))
  duration.hist <- cut(accural.T.dist,
    seq(0,max(uclT)*1.2,length=40))
  barplot(table(duration.hist),horiz=FALSE,
    axes=FALSE,xlab=" ",ylab=" ",space=0,
    col="white",names.arg=rep(" ",39))

  par(mar=c(4.1,4.1,0.1,0.1))
  plot(c(0,max(uclT)*1.2),c(0,np),xlab="Time (Months)",
    ylab="Number of patients",xaxt = 'n',type="n")
  axis(1, at=seq(0,max(uclT)*1.2, 6)) 
  
    
  legenda=paste("Total targeted subjects:", n)
  legendb=paste("Total finish time (months):",T)
  legendd=paste("Time to date (months) :", tm)
  legendc=paste("Subjects recruited to date:",m)
  legende=paste("Subjects in", T, "months:",round(accrual.n.T[2]),"(",round(accrual.n.T[1]),",",round(accrual.n.T[3]),")" )
  legendf=paste("Time for",n, "subjects:", round(midT[np-m],digits=1),"(",round(lclT[np-m],digits=1),",",round(uclT[np-m],digits=1),")" )

  legend(max(uclT)*0.6, np*0.4, legend=c("Input Information:",
         legenda,legendb,legendc,legendd,"------------------------","Summary of Results:",legende,legendf))


  polygon(c(c(tm,lclT),rev(c(tm,uclT))),c(m:np,np:m),
    density=-1,col="gray",border=NA)
  lines(c(tm,midT),m:np,col="white")
  abline(v = T, col = "red") 
  segments(0,0,tm,m)
  return(list(paste("2.5%=",lclT[np-m]), paste("50%=",midT[np-m]),paste("97.5%=",uclT[np-m])))

}

## Function for diagonistic plot
  accrual.plots=function(w){
  par(mfrow=c(2,2))
  qqexp(w,line=TRUE)
  hist(w, freq = FALSE,xlab="Waiting time",main="" )
  lines(dexp(min(w):max(w),1/mean(w)),col = "red", lwd = 2)
  CumTime=cumsum(w)
  plot(CumTime,w,type="l",xlab="Cumulative time (days)",ylab="Waiting Time")
  sub=1:length(w)
  plot(CumTime,sub,type="l",xlab="Cumulative time (days)",ylab="Number of Subjects Recruited")

}



accrual.gui <- function(){
## Gui for Acurral analysis                        ######
accrualcallbackplot<- function(Targeted_sample_size,Targeted_finish_time_in_months,Your_confidence,
                       Total_subjects_recruited_to_date,Time_to_date_in_months) {
 n=Targeted_sample_size
 T=Targeted_finish_time_in_months
 P=Your_confidence
 m=Total_subjects_recruited_to_date
 tm=Time_to_date_in_months
 return <-accrual.n.plot(n,T,P,m,tm,T)
}


durationcallbackplot<- function(Targeted_sample_size,Targeted_finish_time_in_months,Your_confidence,
                       Total_subjects_recruited_to_date,Time_to_date_in_months) {
 n=Targeted_sample_size
 T=Targeted_finish_time_in_months
 P=Your_confidence
 m=Total_subjects_recruited_to_date
 tm=Time_to_date_in_months
 return <-accrual.T.plot(n,T,P,m,tm,n)
}


diagexpgui <- function(w,Header) {
  if (Header==TRUE){ w <- as.matrix( read.csv(w),header=TRUE)}
  if (Header==FALSE){w <- as.matrix( read.csv(w),header=FALSE)}
  return(accrual.plots(w))
  
}

diagguiCallback <- function( arg ) {
  if( arg=="w" ) {
    datanames1 <- names( read.csv( guiGetValue("w") ) )
    print( datanames1 )
    guiSet( "datanames1", datanames1 )
    }
  
  }



mgui(accrualcallbackplot,title=c("Menu","How many subjects will you recruit?"),
     argSlider=list(Your_confidence=c(0,1,0.05)))


mgui(durationcallbackplot,title=c("Menu","How long will it take to reach the targeted sample size?"),
     argSlider=list(Your_confidence=c(0,1,0.05)))

mgui(diagexpgui,title=c("Menu","Diagnostic Plots"),
     argFilename=list(w=NULL), callback=diagguiCallback, argOption=list(Header=c("TRUE","FALSE")) )
}





