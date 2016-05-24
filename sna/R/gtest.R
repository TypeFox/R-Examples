######################################################################
#
# gtest.R
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 2/27/13
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains routines related to null hypothesis testing.
#
# Contents:
#  cug.test
#  cugtest
#  plot.cug.test
#  plot.cugtest
#  plot.qaptest
#  print.cug.test
#  print.cugtest
#  print.qaptest
#  print.summary.cugtest
#  print.summary.qaptest
#  qaptest
#  summary.cugtest
#  summary.qaptest
#
######################################################################


#Second generation (univariate) CUG test routine
cug.test<-function(dat,FUN,mode=c("digraph","graph"),cmode=c("size","edges","dyad.census"),diag=FALSE,reps=1000,ignore.eval=TRUE,FUN.args=list()){
  if(ignore.eval){   #Handle everything in edgelist form here
    #Pre-process the raw input
    dat<-as.edgelist.sna(dat)
    if(is.list(dat))
      return(sapply(dat,cug.test,FUN=FUN,mode=mode,cmode=cmode,diag=diag, reps=reps,ignore.eval=ignore.eval,FUN.args=FUN.args))
    #End pre-processing
    n<-attr(dat,"n")
    if(!diag){
      dat<-dat[dat[,1]!=dat[,2],,drop=FALSE]
      attr(dat,"n")<-n
      ndiag<-0
    }else
      ndiag<-sum(dat[,1]==dat[,2])
    mode<-match.arg(mode)
    cmode<-match.arg(cmode)
    #Generate conditioning statistics
    if(cmode=="size"){
      m<-NULL
      dc<-NULL
    }else if(cmode=="edges"){
      m<-switch(match.arg(mode),
       graph=(NROW(dat)-ndiag)/2+ndiag,
        digraph=NROW(dat)
      )
      dc<-NULL
    }else if(cmode=="dyad.census"){
      m<-NULL
      dc<-dyad.census(dat)
    }
    #Generate randomization functions
    getstat<-function(d){do.call(fun,c(list(d),FUN.args))}
    drawrep<-switch(cmode,
      size=function(n,...){rgraph(n,1,mode=mode,diag=diag,tprob=0.5, return.as.edgelist=TRUE)},
      edges=function(n,m,...){rgnm(n=1,nv=n,m=m,mode=mode,diag=diag, return.as.edgelist=TRUE)},
      dyad.census=function(n,dc,...){rguman(n=1,nv=n,mut=dc[1],asym=dc[2], null=dc[3],method="exact",return.as.edgelist=TRUE)},
    )
  }else{       #For valued data, we're going to use adjacency matrices (sorry)
    #Pre-process the raw input
    dat<-as.sociomatrix.sna(dat)
    if(is.list(dat))
      return(sapply(dat,cug.test,FUN=FUN,mode=mode,cmode=cmode,diag=diag, reps=reps,ignore.eval=ignore.eval,FUN.args=FUN.args))
    else if(length(dim(dat))>2)
      return(apply(dat,1,cug.test,FUN=FUN,mode=mode,cmode=cmode,diag=diag, reps=reps,ignore.eval=ignore.eval,FUN.args=FUN.args))
    #End pre-processing
    n<-NROW(dat)
    m<-NULL
    dc<-NULL
    mode<-match.arg(mode)
    cmode<-match.arg(cmode)
    getstat<-function(d){do.call(fun,c(list(d),FUN.args))}
    drawrep<-switch(cmode,
      size=function(n,...){rgraph(n,1,mode=mode,diag=diag,tprob=0.5)},
      edges=switch(mode,
        digraph=function(n,...){g<-dat; 
          g[upper.tri(g,diag=diag)|lower.tri(g)]<- sample(g[upper.tri(g,diag=diag)|lower.tri(g)]); g
        },
        graph=function(n,...){g<-dat;
          g[upper.tri(g,diag=diag)]<-sample(g[upper.tri(g,diag=diag)]);
          g[lower.tri(g)]<-t(g)[lower.tri(g)]; g
        }
      ),
      dyad.census=function(n,...){g<-rewire.ud(dat,1)[1,,]; 
        if(diag) diag(g)<-sample(diag(g)); g
      },
    )
  }
  #Set things up
  fun<-match.fun(FUN)
  if("mode"%in%names(formals(fun)))
    callmode<-TRUE
  else
    callmode<-FALSE
  if("diag"%in%names(formals(fun)))
    calldiag<-TRUE
  else
    calldiag<-FALSE
  if(callmode)
    FUN.args$mode<-mode
  if(calldiag)
    FUN.args$diag<-diag
  obs<-getstat(dat)
  #Draw replicate stats
  repstats<-vector()
  for(i in 1:reps){
    repstats[i]<-getstat(drawrep(n=n,m=m,dc=dc))
  }
  out<-list(obs.stat=obs,rep.stat=repstats,mode=mode,diag=diag,cmode=cmode, plteobs=mean(repstats<=obs),pgteobs=mean(repstats>=obs),reps=reps)
  class(out)<-"cug.test"
  out
}


#cugtest - Generate, print, and plot CUG (conditional uniform graph) test 
#objects.
cugtest<-function(dat,FUN,reps=1000,gmode="digraph",cmode="density",diag=FALSE,g1=1,g2=2,...){
   out<-list()
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   #End pre-processing
   #First, find the test value for fun on dat
   fun<-match.fun(FUN)
   out$testval<-fun(dat,g1=g1,g2=g2,...)
   #Next, determine on what we are conditioning
   if(cmode=="density"){
      d<-c(gden(dat,g=g1,mode=gmode,diag=diag), gden(dat,g=g2,mode=gmode,diag=diag))
   }else if(cmode=="ties"){
     if(is.list(dat)){
       tie1<-dat[[g1]]
       tie2<-dat[[g2]]
     }else{
       tie1<-dat[g1,,]
       tie2<-dat[g2,,]
     }
   }else{
      d<-c(0.5,0.5)
   }
   if(is.list(dat)){    #Get the graph sizes
     n1<-dim(dat[[g1]])[2]
     n2<-dim(dat[[g2]])[2]
   }else{
     n1<-dim(dat)[2]
     n2<-dim(dat)[2]
   }
   #Now, perform reps replications on random recreations of the data
   out$dist<-vector(mode="numeric",length=reps)
   for(i in 1:reps){
      if(cmode=="ties"){  #Generate random replicates
        dat1<-rgraph(n1,diag=diag,mode=gmode,tielist=tie1)
        dat2<-rgraph(n2,diag=diag,mode=gmode,tielist=tie2)
      }else{
        dat1<-rgraph(n1,tprob=d[1],diag=diag,mode=gmode)
        dat2<-rgraph(n2,tprob=d[2],diag=diag,mode=gmode)
      }
      if(n1==n2){         #Combine into single structure
        datc<-array(dim=c(2,n1,n1))
        datc[1,,]<-dat1
        datc[2,,]<-dat2
      }else
        datc<-list(dat1,dat2)
      out$dist[i]<-fun(datc,g1=1,g2=2,...)  #Compute replicate stat
   }
   #Find p values
   out$pgreq<-mean(as.numeric(out$dist>=out$testval))
   out$pleeq<-mean(as.numeric(out$dist<=out$testval))
   class(out)<-c("cugtest","cug")
   out
}


#plot.cugtest - Plotting method for cugtest
plot.cug.test<-function(x,main="Univariate CUG Test", sub=paste("Conditioning:",x$cmode,"Reps:",x$reps),...){
  xl<-range(c(x$rep.stat,x$obs.stat))
  hist(x$rep.stat,xlim=xl,xlab="CUG Replicates",prob=TRUE,main=main,sub=sub,...)
  abline(v=x$obs.stat,col=2,lwd=3)
}


#plot.cugtest - Plotting method for cugtest
plot.cugtest<-function(x,mode="density",...){
   if(mode=="density"){
      plot(density(x$dist),main="Estimated Density of CUG Replications",xlab="Test Statistic",...)
   }else{
      hist(x$dist,main="Histogram of CUG Replications",xlab="Test Statistic",...)
   }
   abline(v=x$testval,lty=2)
}


#plot.qaptest - Plotting method for qaptest
plot.qaptest<-function(x,mode="density",...){
   if(mode=="density"){
      plot(density(x$dist),main="Estimated Density of QAP Replications",xlab="Test Statistic",...)
   }else{
      hist(x$dist,main="Histogram of QAP Replications",xlab="Test Statistic",...)
   }
   abline(v=x$testval,lty=2)
}


#print.cug.test - Print method for cug.test
print.cug.test<-function(x,...){
  cat("\nUnivariate Conditional Uniform Graph Test\n\n")
  cat("Conditioning Method:",x$cmode,"\nGraph Type:",x$mode,"\nDiagonal Used:", x$diag,"\nReplications:",x$reps,"\n\n")
  cat("Observed Value:",x$obs.stat,"\n")
  cat("Pr(X>=Obs):",x$pgteobs,"\n")
  cat("Pr(X<=Obs):",x$plteobs,"\n\n")
}


#print.cugtest - Print method for cugtest
print.cugtest<-function(x,...){
   cat("\nCUG Test Results\n\n")
   cat("Estimated p-values:\n")
   cat("\tp(f(rnd) >= f(d)):",x$pgreq,"\n")
   cat("\tp(f(rnd) <= f(d)):",x$pleeq,"\n\n")      
}


#print.qaptest - Print method for qaptest
print.qaptest<-function(x,...){
   cat("\nQAP Test Results\n\n")
   cat("Estimated p-values:\n")
   cat("\tp(f(perm) >= f(d)):",x$pgreq,"\n")
   cat("\tp(f(perm) <= f(d)):",x$pleeq,"\n\n")      
}


#print.summary.cugtest - Print method for summary.cugtest
print.summary.cugtest<-function(x,...){
   cat("\nCUG Test Results\n\n")
   cat("Estimated p-values:\n")
   cat("\tp(f(rnd) >= f(d)):",x$pgreq,"\n")
   cat("\tp(f(rnd) <= f(d)):",x$pleeq,"\n")
   cat("\nTest Diagnostics:\n")
   cat("\tTest Value (f(d)):",x$testval,"\n")
   cat("\tReplications:",length(x$dist),"\n")
   cat("\tDistribution Summary:\n")
   cat("\t\tMin:\t",quantile(x$dist,probs=0,names=FALSE),"\n")
   cat("\t\t1stQ:\t",quantile(x$dist,probs=0.25,names=FALSE),"\n")
   cat("\t\tMed:\t",quantile(x$dist,probs=0.5,names=FALSE),"\n")
   cat("\t\tMean:\t",mean(x$dist),"\n")
   cat("\t\t3rdQ:\t",quantile(x$dist,probs=0.75,names=FALSE),"\n")
   cat("\t\tMax:\t",quantile(x$dist,probs=1,names=FALSE),"\n")
   cat("\n")
}


#print.summary.qaptest - Print method for summary.qaptest
print.summary.qaptest<-function(x,...){
   cat("\nQAP Test Results\n\n")
   cat("Estimated p-values:\n")
   cat("\tp(f(perm) >= f(d)):",x$pgreq,"\n")
   cat("\tp(f(perm) <= f(d)):",x$pleeq,"\n")
   cat("\nTest Diagnostics:\n")
   cat("\tTest Value (f(d)):",x$testval,"\n")
   cat("\tReplications:",length(x$dist),"\n")
   cat("\tDistribution Summary:\n")
   cat("\t\tMin:\t",quantile(x$dist,probs=0,names=FALSE),"\n")
   cat("\t\t1stQ:\t",quantile(x$dist,probs=0.25,names=FALSE),"\n")
   cat("\t\tMed:\t",quantile(x$dist,probs=0.5,names=FALSE),"\n")
   cat("\t\tMean:\t",mean(x$dist),"\n")
   cat("\t\t3rdQ:\t",quantile(x$dist,probs=0.75,names=FALSE),"\n")
   cat("\t\tMax:\t",quantile(x$dist,probs=1,names=FALSE),"\n")
   cat("\n")
}


#qaptest - Generate a QAP test object
qaptest<-function(dat,FUN,reps=1000,...){
   out<-list()
   #First, find the test value for fun on dat
   fun<-match.fun(FUN)
   out$testval<-fun(dat,...)
   #Now, perform reps replications on random permutations of the data
   out$dist<-vector(mode="numeric",length=reps)
   for(i in 1:reps){
      out$dist[i]<-fun(rmperm(dat),...)
   }
   #Find p values
   out$pgreq<-mean(as.numeric(out$dist>=out$testval))
   out$pleeq<-mean(as.numeric(out$dist<=out$testval))
   class(out)<-c("qaptest","qap")
   out
}


#summary.cugtest - Summary method for cugtest
summary.cugtest<-function(object, ...){
   out<-object
   class(out)<-c("summary.cugtest",class(out))
   out
}


#summary.qaptest - Summary method for qaptest
summary.qaptest<-function(object, ...){
   out<-object
   class(out)<-c("summary.qaptest",class(out))
   out
}
