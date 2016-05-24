
"agelme" <-
function(depup,depdo,bpup,bpdo, use,weights=c(1,rep(0,length(depup)-1)),vspan=1,k=length(depup)-1,m=2, diagnostic=FALSE)
{#Cagedepth.fun version 1.4R, 26.09.05,is written by Einar Heegaard, Bjerknes Centre for Climate Research, University of Bergen, Allegaten 55, 5007 Bergen, Norway
#e.mail: einar.heegaard@bio.uib.no
  #Use of this function is your responsibility
  if(missing(use))use<-TRUE
  data<-data.frame(depthup=depup,depthdo=depdo,cageup=bpup,cagedo=bpdo, use=use)
  depup<-depup[use]
  depdo<-depdo[use]
  bpup<-bpup[use]
  bpdo<-bpdo[use]
  weights<-weights[use]
  n<-length(weights)
  x<-(depup+((depdo-depup)/2))-min(depup)
  y<-((bpup+bpdo)/2)-min(bpup)
  sd<-abs(((bpup+bpdo)/2)-bpup)
  fit.w<-1/sd
  fit.w[weights==1]<-1
  fit.con<-gam(y~s(x,k=k,m=m),quasi(link="identity",variance="constant"),weights=drop(fit.w),scale=-1)
  fit.mu<-gam(y~s(x,k=k,m=m),quasi(link="identity",variance="mu"),weights=drop(fit.w),scale=-1)
  fit.pl<-list(constant=fit.con,mu=fit.mu)
  if(diagnostic){
    par(mfrow=c(3,4))
    v<-unlist(lapply(fit.pl,function(v){v[[2]]}))
    x1<-range(v)
    x2<-range(sqrt(abs(v)))
    for(i in 1:2)
    {x3<-fit.pl[[i]][[3]]
    x4<-fit.pl[[i]][[2]]
    x5<-fit.pl[[i]]$y
    plot(x3,x4,ylim=x1,xlab="Fitted",ylab="Residuals")
    lines(loess.smooth(x3,x4,span=vspan))
    abline(h=0,lty=2)
    plot(x3,sqrt(abs(x4)),ylim=x2,ylab="Sqrt(abs(Residuals))",xlab="Fitted")
    lines(loess.smooth(x3,sqrt(abs(x4)),span=vspan))
    abline(h=0,lty=2)
    plot(x3,x5,ylab="Observed",xlab="Fitted")
    lines(c(0,max(x5)),c(0,max(x5)),lty=2)
    lines(loess.smooth(x3,x5,span=vspan))
    qqnorm(x4)
    qqline(x4,lty=2)
    }
  }
  age.res<-list(
     tdf=c(sum(fit.con$hat),sum(fit.mu$hat)),
     weights=fit.w,
#     Constant=data.frame(Depth=xp,Calage=y+min(bpup),Estage=yp,Lowlim=yp1,Upplim=yp2,tsd=sd2,Csd=sd,Rsd=v1),
#     Muvar=data.frame(Depth=xpm,Calage=y+min(bpup),Estage=ypm,Lowlim=ypm1,Upplim=ypm2,Tsd=sdm2,Csd=sd,Rsd=vm),
     RES=data.frame(Constvar=sum(fit.con[[2]]^2)/1000,Muvar=sum(fit.mu[[2]]^2)/1000),
     Models=fit.pl,
     data=data
  )
  class(age.res)<-"agelme"
  age.res
}


"predict.agelme" <-
function(object,v = 1,depth, ...)
  {#Cagenew.fun version 1.4R,26.09.05,is written by Einar Heegaard, Bjerknes Centre for Climate Research, University of Bergen, Allegaten 55, 5007 Bergen, Norway
  #e.mail: einar.heegaard@bio.uib.no
  #Use of this function is your responibility
  depup<-object$data[object$data$use,1]
  depdo<-object$data[object$data$use,2]
  bpup<-object$data[object$data$use,3]
  bpdo<-object$data[object$data$use,4]
  fit.m<-object$Models[[v]]
  xa<-data.frame(fit.w=object$weights)
  xd<-depth-min(depup)
#  attach(xa)
  yp<-predict(fit.m,newdata=data.frame(x=xd),type="response")
  ypm<-predict(fit.m,se.fit=T)
#  detach(pos=2)
  sd<-abs(((bpup+bpdo)/2)-bpup)
  x<-(depup+((depdo-depup)/2))
  sd2<-sqrt(sd^2+ypm$se.fit^2)
  sdp<-approx(x,sd2,depth,rule=2)
  sdp<-sdp$y
  y1<-yp-(1.96*sdp)
  y2<-yp+(1.96*sdp)
  vt<-data.frame(depth=c(depth),Estage=c(yp+min(bpup)),Lowlim=c(y1+min(bpup)),Upplim=c(y2+min(bpup)),Tsd=c(sdp))
  vt<-list(v=ifelse(v==1, "Constant variance","Constant variance"), fit=vt, data=object$data)
  class(vt)<-"fittedAgelme"
  vt
}
