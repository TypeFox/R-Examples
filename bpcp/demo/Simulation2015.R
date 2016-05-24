
createDataBeta<-function(n,param.a,param.b){
  x<-rbeta(n,param.a,param.b)
  cens<-runif(n)
  time<-pmin(x,cens)
  status<-rep(0,n)
  status[time==x]<-1
  list(time=time,status=status,xi=x,ci=cens)
}

createDataExpl<-function(n,muX,urange){
  x<-rexp(n,1/muX)
  cens<-runif(n,min=urange[1],max=urange[2])
  time<-pmin(x,cens)
  status<-rep(0,n)
  status[time==x]<-1
  list(time=time,status=status,xi=x,ci=cens)
}

#double triplet plot
dtplot<-function(n1.low,n1.high,n2.low,n2.high,method1,method2,distn,param.a,param.b,byrates){


  if(distn==2)
  {epsname<-paste0(dirx,"XBeta Triplet",method1,method2," n is "
               , n ," reps is ",nsim, " param.a is ", param.a, " param.b is ",param.b,".eps")
  }

  if(distn==1)
  {epsname<-paste0(dirx,"XExpl Triplet",method1,method2," n is "
               , n ," reps is ",nsim, " param.a is ", param.a, " param.b is ",param.b,".eps")
  }

  postscript(file=epsname,horizontal=FALSE)

  P<-byrates[,s.c]
  n1.low.b<-byrates[,n1.low]
  n1.high.b<-byrates[,n1.high]
  n2.low.b<-byrates[,n2.low]
  n2.high.b<-byrates[,n2.high]

  par(mfrow=c(3,2),mar=c(1.5,4,3,0.5),oma=c(0,0,4,0)+.5)
  plot(P,n1.low.b+n1.high.b,type="l",lwd=1,ylim=c(0,.10),ylab="Total Error",
   main=paste0("Error for ",method1))
  lines(c(0,1),c(.05,.05),lty=2,lwd=1,col="gray")
  lines(P,n1.low.b+n1.high.b)

  plot(P,n2.low.b+n2.high.b,type="l",lwd=1,ylim=c(0,.10),ylab="Total Error",
  main=paste0("Error for ",method2))
  lines(c(0,1),c(.05,.05),lty=2,lwd=1,col="gray")
  lines(P,n2.low.b+n2.high.b)

  par(mar=c(1.5,4,1,.5))
  plot(P,n1.high.b,type="l",lwd=1,ylim=c(0,.10),ylab="Upper Limit Error")
  lines(c(0,1),c(.025,.025),lty=2,lwd=1,col="gray")
  lines(P,n1.high.b)

  plot(P,n2.high.b,type="l",lwd=1,ylim=c(0,.10),ylab="Upper Limit Error")
  lines(c(0,1),c(.025,.025),lty=2,lwd=1,col="gray")
  lines(P,n2.high.b)

  plot(P,n1.low.b,type="l",lwd=1,ylim=c(0,.10),ylab="Lower Limit Error")
  lines(c(0,1),c(.025,.025),lty=2,lwd=1,col="gray")
  lines(P,n1.low.b)

  plot(P,n2.low.b,type="l",lwd=1,ylim=c(0,.10),ylab="Lower Limit Error")
  lines(c(0,1),c(.025,.025),lty=2,lwd=1,col="gray")
  lines(P,n2.low.b)


  dev.off()   
}
#end dtplot function

#library(bpcp)
# the line below is because this package sits on my Helix
library(bpcp,lib.loc="/home/ebrittain/mylibrary/")

library(survival)


onestep.mult<-function(muX,n,urange,param.a,param.b,distn){
   if (distn==1) {d<-createDataExpl(n,muX,urange)}
                 
   if (distn==2) {d<-createDataBeta(n,param.a,param.b)}

   #this is just a check on non monotonicity
   bp.r.notmono<-bpcp(d$time,d$status,nmc=0,alpha=.05,Delta=0,stype="km",midp=FALSE,nonIncrease=FALSE)
   bp.m.notmono<-bpcp(d$time,d$status,nmc=0,alpha=.05,Delta=0,stype="km",midp=TRUE,nonIncrease=FALSE)
           
   
     max.r.lower.diff<-max((diff(bp.r.notmono$lower)))
     max.r.upper.diff<-max((diff(bp.r.notmono$upper)))
     max.m.lower.diff<-max((diff(bp.m.notmono$lower)))
     max.m.upper.diff<-max((diff(bp.m.notmono$upper)))  
        diff.mono<-c(max.r.lower.diff,max.r.upper.diff,max.m.lower.diff,max.m.upper.diff,distn,n,muX,urange,param.a,param.b)
       #diff.mono is a vector for each rep


   bp.r<-bpcp(d$time,d$status,nmc=0,alpha=.05,Delta=0,stype="km",midp=FALSE)
   bp.m<-bpcp(d$time,d$status,nmc=0,alpha=.05,Delta=0,stype="km",midp=TRUE)
   mod<-survfit(Surv(time,status)~1,data=d,conf.lower='modified')
   bp.b<-kmciBorkowf(d$time,d$status,type="logs",alpha=.05)
  
   if (distn<=1) {minS<-1-pexp(urange[2],1/muX)}
 
   if (distn==2) {minS<-1/300}
   minS<-max(1/300,minS)

   minSint<-ceiling(minS*numS)
   numrows<-(numS-1)-minSint+1

   numout<-14
   outc<-matrix(NA,numrows,numout)
   for (Strue.m in minSint:(numS-1)) {
     #This computes all values of Strue between minS and 1 by increment 1/numS
     Strue<-Strue.m/numS
     if (distn<=1) {ttrue<-qexp(1-Strue,1/muX)}
     if (distn>=2) {ttrue<-qbeta(1-Strue,param.a,param.b)}
     bsci.r<-StCI(bp.r,ttrue)
     bsci.m<-StCI(bp.m,ttrue)
     bsci.b<-StCI(bp.b,ttrue)
     bsci.mod<-StCI(mod,ttrue,afterMax="zeroNoNA")
   
     bpci.r<-as.vector(as.matrix(bsci.r[,3:4]))
     bpci.m<-as.vector(as.matrix(bsci.m[,3:4]))
     bpci.b<-as.vector(as.matrix(bsci.b[,3:4]))
     mod.test<-as.vector(as.matrix(bsci.mod[,3:4]))
     low.r<-0
     high.r<-0
     low.m<-0
     high.m<-0
     low.b<-0
     high.b<-0
     low.mod<-0
     high.mod<-0
     eps<-.Machine$double.eps^0.5
     if (bpci.r[1]>Strue+eps) low.r<-1
     if (bpci.r[2]<Strue-eps) high.r<-1
     if (bpci.m[1]>Strue+eps) low.m<-1
     if (bpci.m[2]<Strue-eps) high.m<-1
     if (bpci.b[1]>Strue+eps) low.b<-1
     if (bpci.b[2]<Strue-eps) high.b<-1
     if (mod.test[1]>Strue+eps) low.mod<-1
     if (mod.test[2]<Strue-eps) high.mod<-1

     #only relevant to possible alternate versions of Norm Modified TG
     zero.events<-0
     if (bsci.mod[3]==1) {zero.events<-1}

     low.mod.alt<-low.mod
     high.mod.alt<-high.mod
     # The following allows possible alternate way to handle no events
     if (zero.events>.99) {
      low.mod.alt<-0
      high.mod.alt<-0
          } 

     if (distn<=1) {Scens<-max(0,ttrue-urange[1])/(urange[2]-urange[1])}
     if (distn>=2) {Scens<-ttrue}   



     v<-(c(low.r,high.r,low.m,high.m,low.mod,high.mod,low.b,high.b,
       low.mod.alt,high.mod.alt,
       ttrue,Strue,Scens,zero.events))
     outc[Strue.m-minSint+1,]<-v

   }
   output<-outc
   out<-list(outc=outc,diff.mono=diff.mono)
}




#end onestep.mult

dosims<-function(nsim,n,muX,urange,param.a,param.b,distn,seed){
set.seed(seed)
   if (distn==1) {minS<-1-pexp(urange[2],1/muX)}
   if (distn==2) {minS<-1/300}
   minS<-max(1/300,minS)

   minSint<-ceiling(minS*numS)
   numrows<-(numS-1)-minSint+1

   OUTa<-matrix(NA,nsim*numrows,16)
   OUTdiff<-matrix(NA,nsim,11)
   for (rep in 1:nsim){
      #these lines below are just to keep the computer awake
        rep100<-rep
        one<-1
        if (round(rep100/100)==rep100/100) {plot(rep100, one)}
    get<-onestep.mult(muX,n,urange,param.a,param.b,distn)
    col.rep<-t(rep(rep,numrows))
    col.rep<-t(col.rep)
    col.nsim<-t(rep(nsim,numrows))
    col.nsim<-t(col.nsim)
    all<-cbind(get$outc,col.nsim,col.rep)
    OUTa[(numrows*(rep-1)+1):(numrows*(rep-1)+numrows),]<-all
    OUTdiff[rep,]<-get$diff.mono
    }

if (distn==1)  {
dcsvname<-paste0("RPBiowulfKMSimResults/OUTdiff Expl n is ", n ," OUTdiff.csv")
}

if (distn==2)  {
dcsvname<-paste0("RPBiowulfKMSimResults/OUTdiff Beta n is ", n ," a is ", param.a ,"b is ", param.b," OUTdiff.csv")
}


write.csv(OUTdiff,dcsvname)


   output<-OUTa
}
#end dosims

time0<-proc.time()

#muX and urange values only relevant for exponential case
n<-30
nsim<-10000
muX<-10
urange<-c(0,5)

# Set distn=1 for Expl and =2 for Beta
distn<-2
if (distn==1) {numS<-1000}
if (distn==2) {numS<-300}
plot.pdf<-1

dirx<-"RPBiowulfKMSimResults/"

# Handles alternate way to handle testing after 0 events
#if no.reject.zero is 1 -- then norm modified and tg will be non rejections
# the default for no.reject.zero is 0 -- since the official CI is (1,1)
no.reject.zero<-0


t.c<-11
s.c<-12
c.c<-13
nsim.c<-s.c+3
all.c<-s.c+4

#finds values significantly different from .025 given nsim
signif.j<-nsim
found<-0
for (jj in 1:nsim){
j<-binom.test(x = jj, n = nsim, p = 0.025, alternative = "greater")
calc<- (j$p.value) + found
if (calc<.025) {
    signif.j<-jj
    found<- 1}
 }
signif.j

signif.j.low<-nsim
found<-0

valuelow<-seq(nsim,0,-1)
for (jj in valuelow){
j<-binom.test(x = jj, n = nsim, p = 0.025, alternative = "less")
calc<- (j$p.value) + found
if (calc<.025) {
    signif.j.low<-jj
    found<- 1}
 }
if (found<1) {signif.j.low<- -1}
signif.j.low


param.set.all<-function(param.a,param.b,muX,urange,distn) {
  get2<-dosims(nsim=nsim,n=n,param.a=param.a,param.b=param.b,muX=muX,urange=urange,distn=distn,seed=890761)


  #get2[,s.c] are the Strue values from minS to 1
  #We want to get rejection rates BY each value of Strue

  ux.s<-sort(unique(get2[,s.c]))
  nx.s<-length(ux.s)
  byrates<-matrix(NA,nx.s,all.c)
  byrates.sd<-matrix(NA,nx.s,all.c)
  byrates.se<-matrix(NA,nx.s,all.c)
  for (i in 1:nx.s){
    byrates[i,]<-apply(get2[get2[,s.c]==ux.s[i],],MARGIN=2,mean)
    byrates.sd[i,]<-apply(get2[get2[,s.c]==ux.s[i],],MARGIN=2,sd)
    byrates.se[i,]<-byrates.sd[i,]/sqrt(get2[1,nsim.c])
  }
  byrates
  byrates.sd
  byrates.se

  time1<-proc.time()
  time.sim<-time1-time0
  time.sim

  if (distn==1) {param.a<-NA
                 param.b<-NA}

  dtplot(1,2,3,4,"BPCP","Mid-p BPCP",distn,param.a,param.b,byrates)
  dtplot(5,6,7,8,"Modified Lower","Borkowf Shrink",distn,param.a,param.b,byrates)
  dtplot(9,10,7,8,"Modified Alternate","Borkowf Shrink",distn,param.a,param.b,byrates)

  return(byrates)

}
#end param.set.all

if (distn==1) {byrates<-param.set.all(param.a,param.b,muX,urange,distn)
              num.curves<-1}

if (distn==2)
{
byrates.1.100000<-param.set.all(param.a=1,param.b=100000,muX,urange,distn)
byrates.1.7<-param.set.all(param.a=1,param.b=7,muX,urange,distn)
byrates.1.50<-param.set.all(param.a=1,param.b=50,muX,urange,distn)
byrates.1.2<-param.set.all(param.a=1,param.b=2,muX,urange,distn)
byrates.1.1<-param.set.all(param.a=1,param.b=1,muX,urange,distn)
byrates.2.1<-param.set.all(param.a=2,param.b=1,muX,urange,distn)
byrates.7.1<-param.set.all(param.a=7,param.b=1,muX,urange,distn)
byrates.50.1<-param.set.all(param.a=50,param.b=1,muX,urange,distn)
byrates.100000.1<-param.set.all(param.a=100000,param.b=1,muX,urange,distn)
byrates.100.100<-param.set.all(param.a=100,param.b=100,muX,urange,distn)
byrates..1..1<-param.set.all(param.a=.1,param.b=.1,muX,urange,distn)


byrates<-rbind(byrates.1.7,byrates.1.50,byrates.1.100000,byrates.1.2,
   byrates.1.1,byrates.2.1,byrates.7.1,byrates.50.1,byrates.100000.1,
   byrates.100.100,byrates..1..1)

num.curves<-11}


num.byrates<-nrow(byrates)


plotdense<-function(num,title,num.curves){
   m.title<-paste(title)

  plot(x=NA,y=NA,ylim=c(1,0),xlim=c(0,1),xlab="G(t)",ylab="S(t)")
  title(m.title)

  for (k in 1:num.byrates) {
    colz <- sqrt((byrates[k,num])/.025)
    val.g <-0
    val.r <-0
  
    if (byrates[k,num]<=10)  {val.g <- max(0,255 - (byrates[k,num]-.025)*3400)
                  val.g <- min(255,val.g)
                  val.b <- max(0,255 - (byrates[k,num]-.025)*7000)
                  val.b <- min(255,val.b)
                 }

    colr <- rgb(red=255,green=val.g,blue=val.b,maxColorValue=255)

    if (byrates[k,num]<=signif.j.low/nsim) {points(x=byrates[k,c.c],y=byrates[k,s.c],col=gray(colz),pch=19)}
    if (byrates[k,num]>signif.j.low/nsim){points(x=byrates[k,c.c],y=byrates[k,s.c],col=gray(.90),cex=.1)} 
    
    if (byrates[k,num]>=signif.j/nsim) {points(x=byrates[k,c.c],y=byrates[k,s.c],col=colr,pch=19)}
    if (byrates[k,num]>.10) {points(x=byrates[k,c.c],y=byrates[k,s.c],col="red",pch=19)}
  }
}
#end plotdense

plot4dense<-function(n1,title1,n2,title2,n3,title3,n4,title4,num.curves,plotnum){

  epsname<-paste0(dirx,type," Layout EPS 19 Dense Plot",plotnum," n is ", n ," reps is ", nsim,".eps")
  postscript(file=epsname,horizontal=FALSE)

  # This creates the odd number of plots display
     mat<-(matrix(c(1,3,5,2,4,5),3,2))
     layout(mat,heights=c(3,3,1),widths=c(3,3),TRUE)
     par(mgp=c(0,1,0))

  op<-par(no.readonly=TRUE)

  plot2.1<-plotdense(n1,title1,num.curves)
  plot2.2<-plotdense(n2,title2,num.curves)
  plot2.3<-plotdense(n3,title3,num.curves)
  plot2.4<-plotdense(n4,title4,num.curves)

  plot(x=NA,y=NA,xlim=c(0,.10),ylim=c(1,1),yaxt='n',xlab="",ylab="")
  for (p in 0:100) {
    pval<-p/1000
    colz <- sqrt(pval/.025)
    gvalr <-0
    bvalr <-0
    b<-1
  
    if (pval<=10)  {val.g <- max(0,255 - (pval-.025)*3400)
                  val.g <- min(255,val.g)
                  val.b <- max(0,255 - (pval-.025)*7000)
                  val.b <- min(255,val.b)
                 }

    colr <- rgb(red=255,green=val.g,blue=val.b,maxColorValue=255)
  
     if (pval<=signif.j.low/nsim) {points(y=b,x=pval,col=gray(colz),pch=19)}
     if (pval>signif.j.low/nsim){points(y=b,x=pval,col=gray(.90),cex=.1)} 
      
     if (pval>=signif.j/nsim) {points(y=b,x=pval,col=colr,pch=19)}
  
  }

dev.off()
}
#end plot4dense

par(mar=c(3,3,3,3))


if (distn==1) {type<-"Expl"}
if (distn==2) {type<-"Beta"}

plot4dense(1,"BPCP Lower",3,"Mid-p BPCP Lower",5,"Modified Lower Lower",7,"Borkowf Shrink Lower",num.curves,1)
plot4dense(2,"BPCP Upper",4,"Mid-p BPCP Upper",6,"Modified Lower Upper",8,"Borkowf Shrink Upper",num.curves,2)
plot4dense(1,"BPCP Lower",3,"Mid-p BPCP Lower",9,"Modified Alt Lower",7,"Borkowf Shrink Lower",num.curves,3)
plot4dense(2,"BPCP Upper",4,"Mid-p BPCP Upper",10,"Modified Alt Upper",8,"Borkowf Shrink Upper",num.curves,4)

if (distn==1)  {
csvname<-paste0("RPBiowulfKMSimResults/Expl n is ", n ," Results.csv")
}

if (distn==2)  {
csvname<-paste0("RPBiowulfKMSimResults/Beta n is ", n ," Results.csv")
}


write.csv(byrates,csvname)
