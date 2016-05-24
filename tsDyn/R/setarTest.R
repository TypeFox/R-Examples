#' @export
#' 
###tests
#SSR of linear
 # check for stationarty
#call selectSETAR: obtain 2 best thresh


#buildXth1Common / buildXth1NoCommon / buildXth2Common / buildXth2NoCommon
 #   a<-lm()
 # is.InUnitCircle(a) to check
 # deviance(a) for SSR

#compute/return stat

###bootstrap
#bootstrap data
#selectSETAR()

#SSR of those models: SSR_1thresh / SSR_2threshCommon / SSR_2threshNoCommon

#replicate()



setarTest <- function (x, m, d = 1, steps = d, series, thDelay = 0, nboot=10, trim=0.1, test=c("1vs", "2vs3"), hpc=c("none", "foreach"), check=FALSE)
{
  test<-match.arg(test)
  hpc<-match.arg(hpc)
  
  include<-"const" #other types not implemented in setar.sim
  if(missing(series))
      series <- deparse(substitute(x))
###setarTest 1: SSR and check linear model
  linear<-linear(x, m, d = 1, steps = d, series)
  SSR<-deviance(linear)
  

###setarTest 2: search best thresholds: call selectSETAR
  search<-selectSETAR(x, m=m,d=1, steps=d, thDelay=thDelay, trace=FALSE, include =include, common="none", model="TAR",nthresh=2,trim=trim,criterion = "SSR",thSteps = 7,  plot=FALSE,max.iter=3) 

  firstBests<-search$firstBests
  bests<-search$bests
  thDelay<-search$bests[1]

  
  ### Obtain infos for the two thresh models
  set1<-setar(x, m, d=d, steps=steps, thDelay=thDelay, th=firstBests["th"], trace=FALSE, nested=FALSE,include = c("const", "trend","none", "both"), common="none", model=c("TAR", "MTAR"), nthresh=1,trim=trim)
  
   set2<-setar(x, m, d=d, steps=steps, thDelay=thDelay, th=bests[c("th1", "th2")], trace=FALSE, nested=FALSE,include = c("const", "trend","none", "both"), common="none", model=c("TAR", "MTAR"), nthresh=2,trim=trim)
  #RESTRICTED currently:
  #common: FALSE
  #nthresh=2
  #max.iter=1
 
  n<-length(na.omit(residuals(set2)))
    ###setarTest 10: SSR of TAR(+) and (2),compute F stat and print
    SSR1thresh<-search$firstBests["SSR"]
    SSR2thresh<-search$bests["SSR"]
    SSRs<-matrix(c(SSR, SSR1thresh, SSR2thresh), ncol=3, dimnames=list("SSR", c("AR", "TAR(1)", "TAR(2)")))
    
    ###F test for original data
    Ftest12<-as.numeric(n*(SSR-SSR1thresh)/SSR1thresh)
    Bound<-(n/(2*log(log(n))))*Ftest12
    Ftest13<-as.numeric(n*(SSR-SSR2thresh)/SSR2thresh)
    Ftest23<-as.numeric(n*(SSR1thresh-SSR2thresh)/SSR2thresh)
    Ftests<-matrix(c(Ftest12, Ftest13, Ftest23),ncol=3, dimnames=list("Ftest", c("1vs2", "1vs3", "2vs3")))
 if(nboot>0){   
  ########
  ###Boot
  ########
  bootHoAr<-function(linearObject, type){
    ### Bootstrap under null: ar model
    bootLin<-setar.sim(setarObject=linearObject, type=type)$serie
    
  ###setarTest 1: SSR and check linear model
    linearBoot<-linear(bootLin, m=m)
    SSR<-deviance(linearBoot)
    
  ###setarTest 2: search best thresholds: call selectSETAR
    searchBoot<-selectSETAR(bootLin, m=m,d=1, steps=d, thDelay=thDelay, trace=FALSE, include =include, common="none", model=c("TAR", "MTAR"),nthresh=2,trim=trim,criterion = c("SSR"),thSteps = 7,  plot=FALSE,max.iter=3) 
    
    firstBests<-searchBoot$firstBests
    bests<-searchBoot$bests
    
  ###setarTest 10: SSR of all boot models
    SSR1thresh<-searchBoot$firstBests["SSR"]
    SSR2thresh<-searchBoot$bests["SSR"]
      
  ###F test for boot data
    Ftest12<-as.numeric(n*(SSR-SSR1thresh)/SSR1thresh)
    Ftest13<-as.numeric(n*(SSR-SSR2thresh)/SSR2thresh)
    return(c(Ftest12, Ftest13))
  }
  
  bootHoSetar1<-function(setarObject, type ){
    ### Bootstrap under null: Setar(1t) model
    bootSet<-setar.sim(setarObject=set1, type=type, nthresh=1)$serie
    if(check)
      print(all(bootSet-set1$str$x<0.00005))
  
  
  ###setarTest 2: search best thresholds: call selectSETAR
    searchBoot<-selectSETAR(bootSet, m=m,d=1, steps=d, thDelay=thDelay, trace=FALSE, include =include, common="none", model=c("TAR", "MTAR"),nthresh=2,trim=trim,criterion = c("SSR"),thSteps = 7,  plot=FALSE,max.iter=3) 
    
    firstBests<-searchBoot$firstBests
    bests<-searchBoot$bests
    
  ###setarTest 10: SSR of all boot models
    SSR1thresh<-searchBoot$firstBests["SSR"]
    SSR2thresh<-searchBoot$bests["SSR"]
      
  ###F test for boot data
    Ftest23<-as.numeric(n*(SSR1thresh-SSR2thresh)/SSR2thresh)
    
  
  return(Ftest23)
  }
  
  
  ### Run the function (boot, search best SSR for lin, set 1 and set2 ) and extract results
  
  
  type<-ifelse(check, "check", "boot")
  
  probs<-c(0.9, 0.95, 0.975,0.99)
  if(test=="1vs"){
    Ftestboot<-if(hpc=="none"){
      replicate(nboot, bootHoAr(linear, type))
    } else {
      foreach(i=1:nboot, .export="bootHoAr", .combine="rbind") %dopar% bootHoAr(linear, type)
    }
    Ftestboot12<-Ftestboot[1,]
    Ftestboot13<-Ftestboot[2,]
    PvalBoot12<-mean(ifelse(Ftestboot12>Ftest12,1,0))
    CriticalValBoot12<-quantile(Ftestboot12, probs=probs)
    PvalBoot13<-mean(ifelse(Ftestboot13>Ftest13,1,0))
    CriticalValBoot13<-quantile(Ftestboot13, probs=probs)
    CriticalValBoot<-matrix(c(CriticalValBoot12,CriticalValBoot13), nrow=2, byrow=TRUE, dimnames=list(c("1vs2", "1vs3"), probs))
    PvalBoot<-c(PvalBoot12,PvalBoot13)
  }
  else{
    Ftestboot<-replicate(nboot, bootHoSetar1(set1, type))
    Ftestboot<-if(hpc=="none"){
      replicate(nboot, bootHoSetar1(set1, type))
    } else {
      foreach(i=1:nboot, .export="bootHoSetar1", .combine="rbind") %dopar% bootHoSetar1(set1, type)
    }
    Ftestboot23<-Ftestboot
    PvalBoot23<-mean(ifelse(Ftestboot23>Ftest23,1,0))
    CriticalValBoot23<-quantile(Ftestboot23, probs=probs)
    CriticalValBoot<-matrix(CriticalValBoot23, nrow=1, dimnames=list("2vs3", probs))
    PvalBoot<-PvalBoot23
  }
}
  else{ #nboot=0
    CriticalValBoot<-NULL
    PvalBoot<-NULL
    Ftestboot<-NULL
  }
  
  args<-list(x=x, m=m, d = d, steps = steps, series=series, thDelay = thDelay, nboot=nboot, trim=trim, test=test, check=FALSE)
###res
  res<-list(Ftests=Ftests, SSRs=SSRs, firstBests=search$firstBests["th"],secBests=search$bests[c("th1", "th2")], CriticalValBoot=CriticalValBoot,PvalBoot=PvalBoot, Ftestboot=Ftestboot, nboot=nboot, args=args, Bound=Bound, updated=NULL)
  class(res)<-"Hansen99Test"
  return(res)
}
  
#' @S3method print Hansen99Test
print.Hansen99Test<-function(x,...){
  if(x$args$test=="1vs"){
    cat("Test of linearity against setar(2) and setar(3)\n\n")
    print(matrix(c(x$Ftests[-3], x$PvalBoot), ncol=2, dimnames=list(c("1vs2", "1vs3"), c("Test", "Pval"))))
  }
   else{
    cat("Test of setar(2) against  setar(3)\n\n")
    print(matrix(c(x$Ftests[3], x$PvalBoot), ncol=2, dimnames=list(c("2vs3"), c("Test", "Pval"))))
   }
}

#' @S3method summary Hansen99Test
summary.Hansen99Test<-function(object, ...){
  print.Hansen99Test(object)
  cat("\nCritical values:\n")
  print(object$CriticalValBoot)
  cat("\nSSR of original series:\n")
  print(matrix(object$SSRs, ncol=1, dimnames=list(c("AR", "SETAR(2)", "SETAR(3)"), "SSR")))
  cat("\nThreshold of original series:\n")
  print(matrix(c(object$firstBests, NA, object$secBests),byrow=TRUE, ncol=2, dimnames=list(c("SETAR(2)", "SETAR(3)"), c("th1", "th2"))))
  cat("\nNumber of bootstrap replications: ", object$nboot, "\n")
  if(object$args$test=="1vs"){
    cat("Asymptotic bound: ", object$Bound, "\n")
  }
    
}

#' @S3method plot Hansen99Test
plot.Hansen99Test<-function(x,show.extended=TRUE, ...){
  m<-x$args$m
  test<-x$args$test
  leg<-c("Asymptotic Chi 2", "Bootstrap", "Test value")
  col<-c(3,1,2)
  updated<-ifelse(!is.null(x$updated)&show.extended, TRUE, FALSE)
  if(updated){
      leg<-c(leg[1:3], "Old Bootstrap")
      col<-c(col,4)
   }
  
  
  if(test=="1vs"){
    layout(c(1,2))
    Ftestboot12<-x$Ftestboot[1,]
    Ftestboot13<-x$Ftestboot[2,]
    Ftest12<-x$Ftests[1]
    Ftest13<-x$Ftests[2]

    #Plot 1vs2
    plot(density(Ftestboot12, from=0), xlab="Ftest12", xlim=c(0,max(Ftest12+1,max(Ftestboot12))),ylim=c(0,max(density(Ftestboot12)$y,dchisq(0:Ftest12, df=1+m))), main="")
    title("Test linear AR vs 1 threshold SETAR")
    abline(v=Ftest12, lty=2, col=2)
    curve(dchisq(x, df=1+m, ncp=0), from=0, n=Ftest12+5, add=TRUE, col=3)
    if(updated)
      lines(density(Ftestboot12[1:x$updated], from=0), col=4)
    legend("topright", legend=leg, col=col, lty=c(1,1,2), bg="white")
  
    #Plot 1vs3
    plot(density(Ftestboot13, from=0), xlab="Ftest13", xlim=c(0,max(Ftest13+1,max(Ftestboot13))),ylim=c(0,max(density(Ftestboot13)$y,dchisq(0:Ftest12, df=2*(1+m)))), main="")
    title("Test linear AR vs 2 thresholds SETAR")
    abline(v=Ftest13, lty=2, col=2)
    curve(dchisq(x, df=2*(1+m), ncp=0), from=0, n=Ftest13+5, add=TRUE, col=3)
    if(updated)
      lines(density(Ftestboot13[1:x$updated], from=0), col=4)
    legend("topright", legend=leg, col=col, lty=c(1,1,2), bg="white")
  }
  #plot 2vs3
  else {
    Ftestboot23<-x$Ftestboot
    Ftest23<-x$Ftests[3]
    
    plot(density(Ftestboot23, from=0), xlab="Ftest23", xlim=c(0,max(Ftest23+1,Ftestboot23)), ylim=c(0,max(density(Ftestboot23)$y)), main="")
    title("Test 1 threshold SETAR vs 2 thresholds SETAR")
    abline(v=Ftest23, lty=2, col=2)
    curve(dchisq(x, df=1+m, ncp=0), from=0, n=Ftest23+5, add=TRUE, col=3)
    if(updated)
      lines(density(Ftestboot23[1:x$updated], from=0), col=4)
    legend("topright", legend=leg, col=col, lty=c(1,1,2), bg="white")
  }
}



#' @export
extendBoot<-function(x, nboot){
  if(class(x)!="Hansen99Test")
    stop("Function only works for setarTest object")
  args<-x$args
  n <- nboot
  newTestRuns <-setarTest(x=args$x, m=args$m, d = args$d, steps = args$steps, series=args$series, thDelay = args$thDelay, nboot=n, trim=args$trim, test=args$test, check=args$check)
  if(any(x$Ftests!=newTestRuns$Ftests))
    stop("Problem..")
  OldVal<-x$Ftestboot
  NewVal<-newTestRuns$Ftestboot
  
  probs<-c(0.9, 0.95, 0.975,0.99)
  if(args$test=="1vs"){
    Ftestboot<-cbind(OldVal, NewVal)
    Ftestboot12<-Ftestboot[1,]
    Ftestboot13<-Ftestboot[2,]
    PvalBoot12<-mean(ifelse(Ftestboot12>x$Ftests[1],1,0))
    CriticalValBoot12<-quantile(Ftestboot12, probs=probs)
    PvalBoot13<-mean(ifelse(Ftestboot13>x$Ftests[2],1,0))
    CriticalValBoot13<-quantile(Ftestboot13, probs=probs)
    CriticalValBoot<-matrix(c(CriticalValBoot12,CriticalValBoot13), nrow=2, byrow=TRUE, dimnames=list(c("1vs2", "1vs3"), probs))
    PvalBoot<-c(PvalBoot12,PvalBoot13)
  }
  else{
    Ftestboot<-c(OldVal, NewVal)
    Ftestboot23<-Ftestboot
    PvalBoot23<-mean(ifelse(Ftestboot23>x$Ftests[3],1,0))
    CriticalValBoot23<-quantile(Ftestboot23, probs=probs)
    CriticalValBoot<-matrix(CriticalValBoot23, nrow=1, dimnames=list("2vs3", probs))
    PvalBoot<-PvalBoot23
  }
  newNboot<-ifelse(args$test=="1vs", ncol(Ftestboot), length(Ftestboot))
  #second check
  
  if(x$nboot+newTestRuns$nboot!=newNboot)
    stop("Problem2")
  
  res<-list(Ftests=x$Ftests, SSRs=x$SSRs, firstBests=x$firstBests,secBests=x$secBests, CriticalValBoot=CriticalValBoot,PvalBoot=PvalBoot,Ftestboot=Ftestboot, nboot=newNboot, args=args, updated=x$nboot)
  
  class(res)<-"Hansen99Test"
  return(res)
}
  
  
  

if(FALSE){  
library(tsDyn)
sun<-(sqrt(sunspot.year+1)-1)*2

###Sunsport test
#Test 1vs2 and 1vs3
environment(setarTest)<-environment(star)
Han1<-setarTest(sun, m=11, thDelay=0:1, nboot=2,  trim=0.1, test="1vs")

print(Han1)
summary(Han1)
plot(Han1)

#Test 2vs3
environment(setarTest)<-environment(star)
Han2<-setarTest(sun, m=11, thDelay=0:1, nboot=10, plot=TRUE, trim=0.1, test="2vs3")

print(Han2)
summary(Han2)
plot(Han2)

###two probs
#- deviance of selectSETAR and setar is nto the same
# selectSETAR does not select the good one...
environment(selectSETAR)<-environment(star)
selectSETAR(lynx, m=11, thDelay=1, nthresh=2, th=list(exact=c(5.3,8)), criterion="SSR")
selectSETAR(sun, m=11, thDelay=1, nthresh=2, th=list(exact=c(5.3,8)), criterion="SSR")

###US IP
IP<-read.table(file="/media/sda5/Mes documents/Uni/Statistique/Time Series/Handbooks/datasets/Hansen/Hansen1999Linearity/Matlab/ipdat.txt")

#transform as in Hansen 1999
dat<-log(IP[,1])
dat2<-diff(dat, 12)*100 # dat=(dat(13:length(dat(:,1)))-dat(1:length(dat(:,1))-12))*100
dat3<-dat2[157:length(dat2)] #dat=dat(157:length(dat(:,1)));
dat4<-ts(dat3, start=c(1960,1), freq=12)

end(dat4)
length(dat4)
#plot
plot(dat4)

#save
IIPUs<-dat4
save(IIPUs, file="IIPUs.rda")

###Test extendBoot
# test with 10 bootstrap replications:
a<-setarTest(sun[1:100], m=1, nboot=10)
plot(a)

#use old results and compue 20 new replications
b<-extendBoot(a, n=20)
#see the different distributions:
plot(b)

plot(b, show.extended=FALSE)
}
