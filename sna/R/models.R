######################################################################
#
# models.R
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 2/27/13
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains routines related to stochastic models.
#
# Contents:
#   bbnam
#   bbnam.actor
#   bbnam.bf
#   bbnam.fixed
#   bbnam.jntlik
#   bbnam.jntlik.slice
#   bbnam.pooled
#   bn
#   bn.nlpl.dyad
#   bn.nlpl.edge
#   bn.nlpl.triad
#   bn.nltl
#   brokerage
#   coef.bn
#   coef.lnam
#   consensus
#   eval.edgeperturbation
#   lnam
#   nacf
#   netcancor
#   netlm
#   netlogit
#   npostpred
#   potscalered.mcmc
#   plot.bbnam
#   plot.bbnam.actor
#   plot.bbnam.fixed
#   plot.bbnam.pooled
#   plot.bn
#   plot.lnam
#   print.bbnam
#   print.bbnam.actor
#   print.bbnam.fixed
#   print.bbnam.pooled
#   print.bn
#   print.lnam
#   print.netcancor
#   print.netlm
#   print.netlogit
#   print.summary.bbnam
#   print.summary.bbnam.actor
#   print.summary.bbnam.fixed
#   print.summary.bbnam.pooled
#   print.summary.bn
#   print.summary.brokerage
#   print.summary.lnam
#   print.summary.netcancor
#   print.summary.netlm
#   print.summary.netlogit
#   pstar
#   se.lnam
#   summary.bbnam
#   summary.bbnam.actor
#   summary.bbnam.fixed
#   summary.bbnam.pooled
#   summary.bn
#   summary.brokerage
#   summary.lnam
#   summary.netcancor
#   summary.netlm
#   summary.netlogit 
#
######################################################################


#bbnam - Draw from Butts' Bayesian Network Accuracy Model.  This version uses a 
#Gibbs' Sampler, and assumes error rates to be drawn from conditionally 
#independent betas for each actor.  Note that dat MUST be an n x n x n array, 
#and that the data in question MUST be dichotomous.  Priors are also assumed to 
#be in the right form (n x 2 matrices of alpha, beta pairs for em and ep, and
#an n x n probability matrix for the network itself), and are not checked; 
#default behavior if no priors are provided is the uninformative case.
#Wrapper function for the various bbnam models
bbnam<-function(dat,model="actor",...){
   if(model=="actor")
      bbnam.actor(dat,...)
   else if(model=="pooled")
      bbnam.pooled(dat,...)
   else if(model=="fixed")
      bbnam.fixed(dat,...)
}


#bbnam.actor - Draw from the error-prob-by-actor model
bbnam.actor<-function(dat,nprior=0.5,emprior=c(1,11),epprior=c(1,11),diag=FALSE, mode="digraph",reps=5,draws=1500,burntime=500,quiet=TRUE,anames=NULL,onames=NULL,compute.sqrtrhat=TRUE){
   dat<-as.sociomatrix.sna(dat,simplify=TRUE)
   if(is.list(dat))
     stop("All bbnam input graphs must be of the same order.")
   if(length(dim(dat))==2)
     dat<-array(dat,dim=c(1,NROW(dat),NCOL(dat)))
   #First, collect some basic model parameters and do other "setup" stuff
   m<-dim(dat)[1]
   n<-dim(dat)[2]
   d<-dat
   slen<-burntime+floor(draws/reps)
   out<-list()
   if((!is.matrix(nprior))||(NROW(nprior)!=n)||(NCOL(nprior)!=n))
     nprior<-matrix(nprior,n,n)
   if((!is.matrix(emprior))||(NROW(emprior)!=n)||(NCOL(emprior)!=2)){
     if(length(emprior)==2)
       emprior<-sapply(emprior,rep,n)
     else
       emprior<-matrix(emprior,n,2)
   }
   if((!is.matrix(epprior))||(NROW(epprior)!=n)||(NCOL(epprior)!=2)){
     if(length(epprior)==2)
       epprior<-sapply(epprior,rep,n)
     else
       epprior<-matrix(epprior,n,2)
   }
   if(is.null(anames))
     anames<-paste("a",1:n,sep="")
   if(is.null(onames))
     onames<-paste("o",1:m,sep="")
   #Remove any data which doesn't count...
   if(mode=="graph")
      d<-upper.tri.remove(d)
   if(!diag)
      d<-diag.remove(d)
   #OK, let's get started.  First, create temp variables to hold draws, and draw
   #initial conditions for the Markov chain
   if(!quiet)
      cat("Creating temporary variables and drawing initial conditions....\n")
   a<-array(dim=c(reps,slen,n,n))
   em<-array(dim=c(reps,slen,m))
   ep<-array(dim=c(reps,slen,m))
   for(k in 1:reps){
      a[k,1,,]<-rgraph(n,1,diag=diag,mode=mode)
      em[k,1,]<-runif(m,0,0.5)
      ep[k,1,]<-runif(m,0,0.5)
   }
   #Let the games begin: draw from the Gibbs' sampler
   for(i in 1:reps){
      for(j in 2:slen){
         if(!quiet)
            cat("Repetition",i,", draw",j,":\n\tDrawing adjacency matrix\n")
         #Create tie probability matrix
         ep.a<-aperm(array(sapply(ep[i,j-1,],rep,n^2),dim=c(n,n,m)),c(3,2,1))
         em.a<-aperm(array(sapply(em[i,j-1,],rep,n^2),dim=c(n,n,m)),c(3,2,1))
         pygt<-apply(d*(1-em.a)+(1-d)*em.a,c(2,3),prod,na.rm=TRUE)
         pygnt<-apply(d*ep.a+(1-d)*(1-ep.a),c(2,3),prod,na.rm=TRUE)
         tieprob<-(nprior*pygt)/(nprior*pygt+(1-nprior)*pygnt)
         #Draw Bernoulli graph
         a[i,j,,]<-rgraph(n,1,tprob=tieprob,mode=mode,diag=diag)
         if(!quiet)
            cat("\tAggregating binomial counts\n")
         cem<-matrix(nrow=m,ncol=2)
         cep<-matrix(nrow=m,ncol=2)
         for(x in 1:m){
               cem[x,1]<-sum((1-d[x,,])*a[i,j,,],na.rm=TRUE)
               cem[x,2]<-sum(d[x,,]*a[i,j,,],na.rm=TRUE)
               cep[x,1]<-sum(d[x,,]*(1-a[i,j,,]),na.rm=TRUE)
               cep[x,2]<-sum((1-d[x,,])*(1-a[i,j,,]),na.rm=TRUE)
         }
         if(!quiet)
            cat("\tDrawing error parameters\n")
         em[i,j,]<-rbeta(m,emprior[,1]+cem[,1],emprior[,2]+cem[,2])
         ep[i,j,]<-rbeta(m,epprior[,1]+cep[,1],epprior[,2]+cep[,2])
      }
   }
   if(!quiet)
      cat("Finished drawing from Markov chain.  Now computing potential scale reduction statistics.\n")
   if(compute.sqrtrhat){
      out$sqrtrhat<-vector()
      for(i in 1:n)
         for(j in 1:n)
            out$sqrtrhat<-c(out$sqrtrhat,potscalered.mcmc(aperm(a,c(2,1,3,4))[,,i,j]))
      for(i in 1:m)
         out$sqrtrhat<-c(out$sqrtrhat,potscalered.mcmc(aperm(em,c(2,1,3))[,,i]),potscalered.mcmc(aperm(ep,c(2,1,3))[,,i]))
      if(!quiet)
         cat("\tMax potential scale reduction (Gelman et al.'s sqrt(Rhat)) for all scalar estimands:",max(out$sqrtrhat[!is.nan(out$sqrtrhat)],na.rm=TRUE),"\n")
   }
   if(!quiet)
      cat("Preparing output.\n")
   #Whew, we're done with the MCMC.  Now, let's get that data together.
   out$net<-array(dim=c(reps*(slen-burntime),n,n))
   for(i in 1:reps)
      for(j in burntime:slen){
         out$net[(i-1)*(slen-burntime)+(j-burntime),,]<-a[i,j,,]
      }
   if(!quiet)
      cat("\tAggregated network variable draws\n")
   out$em<-em[1,(burntime+1):slen,]
   out$ep<-ep[1,(burntime+1):slen,]
   if(reps>=2)
      for(i in 2:reps){
         out$em<-rbind(out$em,em[i,(burntime+1):slen,])
         out$ep<-rbind(out$ep,ep[i,(burntime+1):slen,])
      }
   if(!quiet)
      cat("\tAggregated error parameters\n")
   #Finish off the output and return it.
   out$anames<-anames
   out$onames<-onames
   out$nactors<-n
   out$nobservers<-m
   out$reps<-reps
   out$draws<-dim(out$em)[1]
   out$burntime<-burntime
   out$model<-"actor"
   class(out)<-c("bbnam.actor","bbnam")
   out
}


#bbnam.bf - Estimate Bayes Factors for the Butts Bayesian Network Accuracy 
#Model.  This implementation relies on monte carlo integration to estimate the 
#BFs, and tests the fixed probability, pooled, and pooled by actor models.
bbnam.bf<-function(dat,nprior=0.5,em.fp=0.5,ep.fp=0.5,emprior.pooled=c(1,11),epprior.pooled=c(1,11),emprior.actor=c(1,11),epprior.actor=c(1,11),diag=FALSE, mode="digraph",reps=1000){
   dat<-as.sociomatrix.sna(dat,simplify=TRUE)
   if(is.list(dat))
     stop("All bbnam.bf input graphs must be of the same order.")
   if(length(dim(dat))==2)
     dat<-array(dat,dim=c(1,NROW(dat),NCOL(dat)))
   n<-dim(dat)[1]
   if((!is.matrix(nprior))||(NROW(nprior)!=n)||(NCOL(nprior)!=n))
     nprior<-matrix(nprior,n,n)
   if((!is.matrix(emprior.actor))||(NROW(emprior.actor)!=n)|| (NCOL(emprior.actor)!=2)){
     if(length(emprior.actor)==2)
       emprior.actor<-sapply(emprior.actor,rep,n)
     else
       emprior.actor<-matrix(emprior.actor,n,2)
   }
   if((!is.matrix(epprior.actor))||(NROW(epprior.actor)!=n)|| (NCOL(epprior.actor)!=2)){
     if(length(epprior.actor)==2)
       epprior.actor<-sapply(epprior.actor,rep,n)
     else
       epprior.actor<-matrix(epprior.actor,n,2)
   }
   d<-dat
   if(!diag)
      d<-diag.remove(d)
   if(mode=="graph")
      d<-lower.tri.remove(d)
   pfpv<-vector()
   ppov<-vector()
   pacv<-vector()
   #Draw em, ep, and a values for the various models
   for(i in 1:reps){
      a<-rgraph(n,1,tprob=nprior)
      em.pooled<-eval(call("rbeta",1,emprior.pooled[1],emprior.pooled[2]))
      ep.pooled<-eval(call("rbeta",1,epprior.pooled[1],epprior.pooled[2]))
      em.actor<-eval(call("rbeta",n,emprior.actor[,1],emprior.actor[,2]))
      ep.actor<-eval(call("rbeta",n,epprior.actor[,1],epprior.actor[,2]))
      pfpv[i]<-bbnam.jntlik(d,a=a,em=em.fp,ep=ep.fp,log=TRUE)
      ppov[i]<-bbnam.jntlik(d,a=a,em=em.pooled,ep=ep.pooled,log=TRUE)
      pacv[i]<-bbnam.jntlik(d,a=a,em=em.actor,ep=ep.actor,log=TRUE)
   }
   int.lik<-c(logMean(pfpv),logMean(ppov),logMean(pacv))
   int.lik.std<-sqrt(c(var(pfpv),var(ppov),var(pacv)))
   int.lik.std<-(logSub(c(logMean(2*pfpv),logMean(2*ppov),logMean(2*pacv)), 2*int.lik)-log(reps))/2
   #Find the Bayes Factors
   o<-list()
   o$int.lik<-matrix(nrow=3,ncol=3)
   for(i in 1:3)
      for(j in 1:3){
         if(i!=j)
            o$int.lik[i,j]<-int.lik[i]-int.lik[j]
         else
            o$int.lik[i,i]<-int.lik[i]
      }
   o$int.lik.std<-int.lik.std
   o$reps<-reps
   o$prior.param<-list(nprior,em.fp,ep.fp,emprior.pooled,epprior.pooled,emprior.actor,epprior.actor)
   o$prior.param.names<-c("nprior","em.fp","ep.fp","emprior.pooled","epprior.pooled","emprior.actor","epprior.actor")
   o$model.names<-c("Fixed Error Prob","Pooled Error Prob","Actor Error Prob")
   class(o)<-c("bbnam.bf","bayes.factor")
   o
}


#bbnam.fixed - Draw from the fixed probability error model
bbnam.fixed<-function(dat,nprior=0.5,em=0.25,ep=0.25,diag=FALSE,mode="digraph",draws=1500,outmode="draws",anames=NULL,onames=NULL){
   dat<-as.sociomatrix.sna(dat,simplify=TRUE)
   if(is.list(dat))
     stop("All bbnam input graphs must be of the same order.")
   if(length(dim(dat))==2)
     dat<-array(dat,dim=c(1,NROW(dat),NCOL(dat)))
   #How many actors are involved?
   m<-dim(dat)[1]
   n<-dim(dat)[2]
   if((!is.matrix(nprior))||(NROW(nprior)!=n)||(NCOL(nprior)!=n))
     nprior<-matrix(nprior,n,n)
   if(is.null(anames))
     anames<-paste("a",1:n,sep="")
   if(is.null(onames))
     onames<-paste("o",1:m,sep="")
   #Check to see if we've been given full matrices (or vectors) of error probs...
   if(length(em)==m*n^2)
      em.a<-em
   else if(length(em)==n^2)
      em.a<-apply(em,c(1,2),rep,m)
   else if(length(em)==m)
      em.a<-aperm(array(sapply(em,rep,n^2),dim=c(n,n,m)),c(3,2,1))
   else if(length(em)==1)
      em.a<-array(rep(em,m*n^2),dim=c(m,n,n))
   if(length(ep)==m*n^2)
      ep.a<-ep
   else if(length(ep)==n^2)
      ep.a<-apply(ep,c(1,2),rep,m)
   else if(length(ep)==m)
      ep.a<-aperm(array(sapply(ep,rep,n^2),dim=c(n,n,m)),c(3,2,1))
   else if(length(ep)==1)
      ep.a<-array(rep(ep,m*n^2),dim=c(m,n,n))
   #Find the network posterior
   pygt<-apply(dat*(1-em.a)+(1-dat)*em.a,c(2,3),prod,na.rm=TRUE)
   pygnt<-apply(dat*ep.a+(1-dat)*(1-ep.a),c(2,3),prod,na.rm=TRUE)
   npost<-(nprior*pygt)/(nprior*pygt+(1-nprior)*pygnt)
   #Send the needed output
   if(outmode=="posterior")
      npost
   else{
      o<-list()
      o$net<-rgraph(n,draws,tprob=npost,diag=diag,mode=mode)
      o$anames<-anames
      o$onames<-onames
      o$nactors<-n
      o$nobservers<-m
      o$draws<-draws
      o$model<-"fixed"
      class(o)<-c("bbnam.fixed","bbnam")
      o
   }
}


#bbnam.jntlik - An internal function for bbnam
bbnam.jntlik<-function(dat,log=FALSE,...){
   p<-sum(sapply(1:dim(dat)[1],bbnam.jntlik.slice,dat=dat,log=TRUE,...))
   if(!log)
      exp(p)
   else
      p
}


#bbnam.jntlik.slice - An internal function for bbnam
bbnam.jntlik.slice<-function(s,dat,a,em,ep,log=FALSE){
   if(length(em)>1)
      em.l<-em[s]
   else
      em.l<-em
   if(length(ep)>1)
      ep.l<-ep[s]
   else
      ep.l<-ep
   p<-sum(log((1-a)*(dat[s,,]*ep.l+(1-dat[s,,])*(1-ep.l))+a*(dat[s,,]*(1-em.l)+(1-dat[s,,])*em.l)),na.rm=TRUE)
   if(!log)
      exp(p)
   else
      p
}


#bbnam.pooled - Draw from the pooled error model
bbnam.pooled<-function(dat,nprior=0.5,emprior=c(1,11),epprior=c(1,11), diag=FALSE,mode="digraph",reps=5,draws=1500,burntime=500,quiet=TRUE,anames=NULL,onames=NULL,compute.sqrtrhat=TRUE){
   dat<-as.sociomatrix.sna(dat,simplify=TRUE)
   if(is.list(dat))
     stop("All bbnam input graphs must be of the same order.")
   if(length(dim(dat))==2)
     dat<-array(dat,dim=c(1,NROW(dat),NCOL(dat)))
   #First, collect some basic model parameters and do other "setup" stuff
   m<-dim(dat)[1]
   n<-dim(dat)[2]
   d<-dat
   slen<-burntime+floor(draws/reps)
   out<-list()
   if((!is.matrix(nprior))||(NROW(nprior)!=n)||(NCOL(nprior)!=n))
     nprior<-matrix(nprior,n,n)
   if(is.null(anames))
     anames<-paste("a",1:n,sep="")
   if(is.null(onames))
     onames<-paste("o",1:m,sep="")
   #Remove any data which doesn't count...
   if(mode=="graph")
      d<-upper.tri.remove(d)
   if(!diag)
      d<-diag.remove(d)
   #OK, let's get started.  First, create temp variables to hold draws, and draw
   #initial conditions for the Markov chain
   if(!quiet)
      cat("Creating temporary variables and drawing initial conditions....\n")
   a<-array(dim=c(reps,slen,n,n))
   em<-array(dim=c(reps,slen))
   ep<-array(dim=c(reps,slen))
   for(k in 1:reps){
      a[k,1,,]<-rgraph(n,1,diag=diag,mode=mode)
      em[k,1]<-runif(1,0,0.5)
      ep[k,1]<-runif(1,0,0.5)
   }
   #Let the games begin: draw from the Gibbs' sampler
   for(i in 1:reps){
      for(j in 2:slen){
         if(!quiet)
            cat("Repetition",i,", draw",j,":\n\tDrawing adjacency matrix\n")
         #Create tie probability matrix
         ep.a<-array(rep(ep[i,j-1],m*n^2),dim=c(m,n,n))
         em.a<-array(rep(em[i,j-1],m*n^2),dim=c(m,n,n))
         pygt<-apply(d*(1-em.a)+(1-d)*em.a,c(2,3),prod,na.rm=TRUE)
         pygnt<-apply(d*ep.a+(1-d)*(1-ep.a),c(2,3),prod,na.rm=TRUE)
         tieprob<-(nprior*pygt)/(nprior*pygt+(1-nprior)*pygnt)
         #Draw Bernoulli graph
         a[i,j,,]<-rgraph(n,1,tprob=tieprob,mode=mode,diag=diag)
         if(!quiet)
            cat("\tAggregating binomial counts\n")
         cem<-vector(length=2)
         cep<-vector(length=2)
         a.a<-apply(a[i,j,,],c(1,2),rep,m)
         cem[1]<-sum((1-d)*a.a,na.rm=TRUE)
         cem[2]<-sum(d*a.a,na.rm=TRUE)
         cep[1]<-sum(d*(1-a.a),na.rm=TRUE)
         cep[2]<-sum((1-d)*(1-a.a),na.rm=TRUE)
         #cat("em - alpha",cem[1],"beta",cem[2]," ep - alpha",cep[1],"beta",cep[2],"\n")
         if(!quiet)
            cat("\tDrawing error parameters\n")
         em[i,j]<-rbeta(1,emprior[1]+cem[1],emprior[2]+cem[2])
         ep[i,j]<-rbeta(1,epprior[1]+cep[1],epprior[2]+cep[2])
      }
   }
   if(!quiet)
      cat("Finished drawing from Markov chain.  Now computing potential scale reduction statistics.\n")
   if(compute.sqrtrhat){
      out$sqrtrhat<-vector()
      for(i in 1:n)
         for(j in 1:n)
            out$sqrtrhat<-c(out$sqrtrhat,potscalered.mcmc(aperm(a,c(2,1,3,4))[,,i,j]))
      out$sqrtrhat<-c(out$sqrtrhat,potscalered.mcmc(em),potscalered.mcmc(ep))
   if(!quiet)
      cat("\tMax potential scale reduction (Gelman et al.'s sqrt(Rhat)) for all scalar estimands:",max(out$sqrtrhat[!is.nan(out$sqrtrhat)],na.rm=TRUE),"\n")
   }
   if(!quiet)
      cat("Preparing output.\n")
   #Whew, we're done with the MCMC.  Now, let's get that data together.
   out$net<-array(dim=c(reps*(slen-burntime),n,n))
   for(i in 1:reps)
      for(j in burntime:slen){
         out$net[(i-1)*(slen-burntime)+(j-burntime),,]<-a[i,j,,]
      }
   if(!quiet)
      cat("\tAggregated network variable draws\n")
   out$em<-em[1,(burntime+1):slen]
   out$ep<-ep[1,(burntime+1):slen]
   if(reps>=2)
      for(i in 2:reps){
         out$em<-c(out$em,em[i,(burntime+1):slen])
         out$ep<-c(out$ep,ep[i,(burntime+1):slen])
      }
   if(!quiet)
      cat("\tAggregated error parameters\n")
   #Finish off the output and return it.
   out$anames<-anames
   out$onames<-onames
   out$nactors<-n
   out$nobservers<-m
   out$reps<-reps
   out$draws<-length(out$em)
   out$burntime<-burntime
   out$model<-"pooled"
   class(out)<-c("bbnam.pooled","bbnam")
   out
}


#bbnam.probtie - Probability of a given tie
bbnam.probtie<-function(dat,i,j,npriorij,em,ep){
   num<-npriorij
   denom<-1-npriorij
   num<-num*prod(dat[,i,j]*(1-em)+(1-dat[,i,j])*em,na.rm=TRUE)
   denom<-denom*prod(dat[,i,j]*ep+(1-dat[,i,j])*(1-ep),na.rm=TRUE)
   p<-num/(denom+num)
   p
}


#bn - Fit a biased net model
bn<-function(dat,method=c("mple.triad","mple.dyad","mple.edge","mtle"),param.seed=NULL,param.fixed=NULL,optim.method="BFGS",optim.control=list(),epsilon=1e-5){
  dat<-as.sociomatrix.sna(dat,simplify=FALSE)
  if(is.list(dat))
    return(lapply(dat,bn,method=method,param.seed=param.seed, param.fixed=param.fixed,optim.method=optim.method,optim.contol=optim.control,epsilon=epsilon))
  else if(length(dim(dat))>2)
    return(apply(dat,1,bn,method=method,param.seed=param.seed, param.fixed=param.fixed,optim.method=optim.method,optim.contol=optim.control,epsilon=epsilon))
  n<-NROW(dat)
  #Make sure dat is appropriate
  if(!is.matrix(dat))
    stop("Adjacency matrix required in bn.")
  dat<-dat>0
  #Choose the objective function to use
  nll<-switch(match.arg(method),
    mple.edge=match.fun("bn.nlpl.edge"),
    mple.dyad=match.fun("bn.nlpl.dyad"),
    mple.triad=match.fun("bn.nlpl.triad"),
    mtle=match.fun("bn.nltl")
  )
  #Extract the necessary sufficient statistics
  if(match.arg(method)%in%c("mple.edge","mple.dyad")){  #Use dyad census stats
    stats<-matrix(0,nrow=n-1,ncol=4)
    stats<-matrix(.C("bn_dyadstats_R",as.integer(dat),as.double(n), stats=as.double(stats),PACKAGE="sna")$stats,ncol=4)
    stats<-stats[apply(stats[,2:4],1,sum)>0,]  #Strip uneeded rows
  }else if(match.arg(method)=="mple.triad"){          #Use full dyad stats
    stats<-matrix(0,nrow=n,ncol=n)
    stats<-matrix(.C("bn_triadstats_R",as.integer(dat),as.double(n), stats=as.double(stats),PACKAGE="sna")$stats,nrow=n,ncol=n)
  }else if(match.arg(method)=="mtle"){                #Use triad census stats
    stats<-as.vector(triad.census(dat))  #Obtain triad census
  }
  #Initialize parameters (using crudely reasonable values)
  if(is.null(param.seed))
    param<-c(gden(dat),grecip(dat,measure="edgewise"),gtrans(dat),gtrans(dat))
  else{
    param<-c(gden(dat),grecip(dat,measure="edgewise"),gtrans(dat),gtrans(dat))
    if(!is.null(param.seed$pi))
      param[1]<-param.seed$pi
    if(!is.null(param.seed$sigma))
      param[2]<-param.seed$sigma
    if(!is.null(param.seed$rho))
      param[3]<-param.seed$rho
    if(!is.null(param.seed$d))
      param[4]<-param.seed$d
  }
  if(is.null(param.fixed))   #Do we need to fix certain parameter values?
    fixed<-rep(NA,4)
  else{
    fixed<-rep(NA,4)
    if(!is.null(param.fixed$pi))
      fixed[1]<-param.fixed$pi
    if(!is.null(param.fixed$sigma))
      fixed[2]<-param.fixed$sigma
    if(!is.null(param.fixed$rho))
      fixed[3]<-param.fixed$rho
    if(!is.null(param.fixed$d))
      fixed[4]<-param.fixed$d
  }
  param<-pmax(pmin(param,1-epsilon),epsilon)  #Ensure interior starting vals
  param<-log(param/(1-param))  #Transform to logit scale
  #Fit the model
  fit<-optim(param,nll,method=optim.method,control=optim.control,stats=stats, fixed=fixed,dat=dat)
  fit$par<-1/(1+exp(-fit$par))                 #Untransform
  fit$par[!is.na(fixed)]<-fixed[!is.na(fixed)] #Fix
  #Prepare the results
  out<-list(d=fit$par[4],pi=fit$par[1],sigma=fit$par[2],rho=fit$par[3], method=match.arg(method),G.square=2*fit$value,epsilon=epsilon)
  #Add GOF for triads
  if(match.arg(method)=="mtle")
    out$triads<-stats
  else
    out$triads<-as.vector(triad.census(dat))
  out$triads.pred<-.C("bn_ptriad_R", as.double(out$pi),as.double(out$sigma),as.double(out$rho), as.double(out$d),pt=as.double(rep(0,16)),PACKAGE="sna")$pt
  names(out$triads.pred)<-c("003", "012", "102", "021D", "021U", "021C", "111D", "111U", "030T", "030C", "201", "120D", "120U", "120C", "210", "300")
  names(out$triads)<-names(out$triads.pred)
  #Add GOF for dyads, using triad distribution
  if(match.arg(method)%in%c("mple.edge","mple.dyad"))
    out$dyads<-apply(stats[,-1],2,sum)
  else
    out$dyads<-as.vector(dyad.census(dat))
  out$dyads.pred<-c(sum(out$triads.pred*c(0,0,1,0,0,0,1,1,0,0,2,1,1,1,2,3)), sum(out$triads.pred*c(0,1,0,2,2,2,1,1,3,3,0,2,2,2,1,0)), sum(out$triads.pred*c(3,2,2,1,1,1,1,1,0,0,1,0,0,0,0,0)))*choose(n,3)/choose(n,2)/(n-2)
  names(out$dyads.pred)<-c("Mut","Asym","Null")
  names(out$dyads)<-names(out$dyads.pred)
  #Add GOF for edges, using dyad distribution
  out$edges<-c(2*out$dyads[1]+out$dyads[2],2*out$dyads[3]+out$dyads[2])
  out$edges.pred<-c(2*out$dyads.pred[1]+out$dyads.pred[2], 2*out$dyads.pred[3]+out$dyads.pred[2])/2
  names(out$edges.pred)<-c("Present","Absent")
  names(out$edges)<-names(out$edges.pred)
  #Add predicted structure statistics (crude)
  a<-out$d*(n-1)
  out$ss.pred<-c(1/n,(1-1/n)*(1-exp(-a/n)))
  for(i in 2:(n-1))
    out$ss.pred<-c(out$ss.pred,(1-sum(out$ss.pred[1:i])) * (1-exp(-(a-out$pi-out$sigma*(a-1))*out$ss.pred[i])))
  out$ss.pred<-cumsum(out$ss.pred)
  names(out$ss.pred)<-0:(n-1)
  out$ss<-structure.statistics(dat)
  #Return the result
  class(out)<-"bn"
  out
}


#bn.nlpl.dyad - Compute the dyadic -log pseudolikelihood for a biased net model
bn.nlpl.dyad<-function(p,stats,fixed=rep(NA,4),...){
  p<-1/(1+exp(-p))
  #Correct for any fixed parameters
  p[!is.na(fixed)]<-fixed[!is.na(fixed)]
  #Calculate the pseudolikelihood
  lpl<-0
  lpl<-.C("bn_lpl_dyad_R",as.double(stats),as.double(NROW(stats)), as.double(p[1]),as.double(p[2]),as.double(p[3]),as.double(p[4]), lpl=as.double(lpl),PACKAGE="sna")$lpl
  -lpl
}


#bn.nlpl.edge - Compute the -log pseudolikelihood for a biased net model,
#using the directed edge pseudolikelihood.  Not sure why you'd want to
#do this, except as a check on the other results....
bn.nlpl.edge<-function(p,stats,fixed=rep(NA,4),...){
  p<-1/(1+exp(-p))
  #Correct for any fixed parameters
  p[!is.na(fixed)]<-fixed[!is.na(fixed)]
  #Calculate the pseudolikelihood
  lp<-cbind(1-(1-p[1])*((1-p[3])^stats[,1])*((1-p[2])^stats[,1])*(1-p[4]),
    1-((1-p[2])^stats[,1])*(1-p[4]))
  lp<-log(cbind(lp,1-lp))
  lpl<-cbind(2*stats[,2]*lp[,1],stats[,3]*lp[,2],stats[,3]*lp[,3], 2*stats[,4]*lp[,4])
  lpl[is.nan(lpl)]<-0  #Treat 0 * -Inf as 0
  -sum(lpl)
}


#bn.nlpl.triad - Compute the triadic -log pseudolikelihood for a biased net
#model, using the Skvoretz (2003) working paper method
bn.nlpl.triad<-function(p,dat,stats,fixed=rep(NA,4),...){
  p<-1/(1+exp(-p))
  #Correct for any fixed parameters
  p[!is.na(fixed)]<-fixed[!is.na(fixed)]
  #Calculate the pseudolikelihood
  lpl<-0
  lpl<-.C("bn_lpl_triad_R",as.integer(dat),as.double(stats), as.double(NROW(stats)),as.double(p[1]),as.double(p[2]),as.double(p[3]), as.double(p[4]), lpl=as.double(lpl),PACKAGE="sna")$lpl
  -lpl
}


#bn.nltl - Compute the -log triad likelihood for a biased net model
bn.nltl<-function(p,stats,fixed=rep(NA,4),...){
  p<-1/(1+exp(-p))
  #Correct for any fixed parameters
  p[!is.na(fixed)]<-fixed[!is.na(fixed)]
  #Calculate the triad likelihood
  pt<-rep(0,16)
  triprob<-.C("bn_ptriad_R", as.double(p[1]),as.double(p[2]),as.double(p[3]), as.double(p[4]),pt=as.double(pt),PACKAGE="sna")$pt
  -sum(stats*log(triprob))
}


#brokerage - perform a Gould-Fernandez brokerage analysis
brokerage<-function(g,cl){
  #Pre-process the raw input
  g<-as.edgelist.sna(g)
  if(is.list(g))
    return(lapply(g,brokerage,cl))
  #End pre-processing
  N<-attr(g,"n")
  m<-NROW(g)
  classes<-unique(cl)
  icl<-match(cl,classes)
  #Compute individual brokerage measures
  br<-matrix(0,N,5)
  br<-matrix(.C("brokerage_R",as.double(g),as.integer(N),as.integer(m), as.integer(icl), brok=as.double(br),PACKAGE="sna",NAOK=TRUE)$brok,N,5)
  br<-cbind(br,apply(br,1,sum))
  #Global brokerage measures
  gbr<-apply(br,2,sum)
  #Calculate expectations and such
  d<-m/(N*(N-1))
  clid<-unique(cl)         #Count the class memberships
  n<-vector()
  for(i in clid)
    n<-c(n,sum(cl==i))
  n<-as.double(n)          #This shouldn't be needed, but R will generate
  N<-as.double(N)          #integer overflows unless we coerce to double!
  ebr<-matrix(0,length(clid),6)
  vbr<-matrix(0,length(clid),6)
  for(i in 1:length(clid)){  #Compute moments by broker's class
    #Type 1: Within-group (wI)
    ebr[i,1]<-d^2*(1-d)*(n[i]-1)*(n[i]-2)
    vbr[i,1]<-ebr[i,1]*(1-d^2*(1-d))+2*(n[i]-1)*(n[i]-2)*(n[i]-3)*d^3*(1-d)^3
    #Type 2: Itinerant (WO)
    ebr[i,2]<-d^2*(1-d)*sum(n[-i]*(n[-i]-1))
    vbr[i,2]<-ebr[i,2]*(1-d^2*(1-d))+ 2*sum(n[-i]*(n[-i]-1)*(n[-i]-2))*d^3*(1-d)^3
    #Type 3: Representative (bIO)
    ebr[i,3]<-d^2*(1-d)*(N-n[i])*(n[i]-1)
    vbr[i,3]<-ebr[i,3]*(1-d^2*(1-d))+ 2*((n[i]-1)*choose(N-n[i],2)+(N-n[i])*choose(n[i]-1,2))*d^3*(1-d)^3
    #Type 4: Gatekeeping (bOI)
    ebr[i,4]<-ebr[i,3]
    vbr[i,4]<-vbr[i,3]
    #Type 5: Liason (bO)
    ebr[i,5]<-d^2*(1-d)*(sum((n[-i])%o%(n[-i]))-sum(n[-i]^2))
    vbr[i,5]<-ebr[i,5]*(1-d^2*(1-d))+ 4*sum(n[-i]*choose(N-n[-i]-n[i],2)*d^3*(1-d)^3)
    #Total
    ebr[i,6]<-d^2*(1-d)*(N-1)*(N-2)
    vbr[i,6]<-ebr[i,6]*(1-d^2*(1-d))+2*(N-1)*(N-2)*(N-3)*d^3*(1-d)^3
  }
  br.exp<-vector()
  br.sd<-vector()
  for(i in 1:N){
    temp<-match(cl[i],clid)
    br.exp<-rbind(br.exp,ebr[temp,])
    br.sd<-rbind(br.sd,sqrt(vbr[temp,]))
  }
  br.z<-(br-br.exp)/br.sd
  egbr<-vector()                     #Global expections/variances
  vgbr<-vector()
  #Type 1: Within-group (wI)
  egbr[1]<-d^2*(1-d)*sum(n*(n-1)*(n-2))
  vgbr[1]<-egbr[1]*(1-d^2*(1-d))+ sum(n*(n-1)*(n-2)*(((4*n-10)*d^3*(1-d)^3)-(4*(n-3)*d^4*(1-d)^2)+((n-3)*d^5*(1-d))))
  #Type 2: Itinerant (WO)
  egbr[2]<-d^2*(1-d)*sum(n*(N-n)*(n-1))
  vgbr[2]<-egbr[2]*(1-d^2*(1-d))+ (sum(outer(n,n,function(x,y){x*y*(x-1)*(((2*x+2*y-6)*d^3*(1-d)^3)+((N-x-1)*d^5*(1-d)))})) - sum(n*n*(n-1)*(((4*n-6)*d^3*(1-d)^3)+((N-n-1)*d^5*(1-d)))))
  #Type 3: Representative (bIO)
  egbr[3]<-d^2*(1-d)*sum(n*(N-n)*(n-1))
  vgbr[3]<-egbr[3]*(1-d^2*(1-d))+ sum(n*(N-n)*(n-1)*(((N-3)*d^3*(1-d)^3)+((n-2)*d^5*(1-d))))
  #Type 4: Gatekeeping (bOI)
  egbr[4]<-egbr[3]
  vgbr[4]<-vgbr[3]
  #Type 5: Liason (bO)
  egbr[5]<- d^2*(1-d)*(sum(outer(n,n,function(x,y){x*y*(N-x-y)}))-sum(n*n*(N-2*n)))
  vgbr[5]<-egbr[5]*(1-d^2*(1-d))
  for(i in 1:length(n))
    for(j in 1:length(n))
      for(k in 1:length(n))
        if((i!=j)&&(j!=k)&&(i!=k))
          vgbr[5]<-vgbr[5] + n[i]*n[j]*n[k] * ((4*(N-n[j])-2*(n[i]+n[k]+1))*d^3*(1-d)^3-(4*(N-n[k])-2*(n[i]+n[j]+1))*d^4*(1-d)^2+(N-(n[i]+n[k]+1))*d^5*(1-d))
  #Total
  egbr[6]<-d^2*(1-d)*N*(N-1)*(N-2)
  vgbr[6]<-egbr[6]*(1-d^2*(1-d))+ N*(N-1)*(N-2)*(((4*N-10)*d^3*(1-d)^3)-(4*(N-3)*d^4*(1-d)^2)+((N-3)*d^5*(1-d)))
  
  #Return the results
  br.nam<-c("w_I","w_O","b_IO","b_OI","b_O","t")
  colnames(br)<-br.nam
  rownames(br)<-attr(g,"vnames")
  colnames(br.exp)<-br.nam
  rownames(br.exp)<-attr(g,"vnames")
  colnames(br.sd)<-br.nam
  rownames(br.sd)<-attr(g,"vnames")
  colnames(br.z)<-br.nam
  rownames(br.z)<-attr(g,"vnames")
  names(gbr)<-br.nam
  names(egbr)<-br.nam
  names(vgbr)<-br.nam
  colnames(ebr)<-br.nam
  rownames(ebr)<-clid
  colnames(vbr)<-br.nam
  rownames(vbr)<-clid
  out<-list(raw.nli=br,exp.nli=br.exp,sd.nli=br.sd,z.nli=br.z,raw.gli=gbr, exp.gli=egbr,sd.gli=sqrt(vgbr),z.gli=(gbr-egbr)/sqrt(vgbr),exp.grp=ebr, sd.grp=sqrt(vbr),cl=cl,clid=clid,n=n,N=N)
  class(out)<-"brokerage"
  out
}


#coef.bn - Coefficient method for bn
coef.bn<-function(object, ...){
  coef<-c(object$d,object$pi,object$sigma,object$rho)
  names(coef)<-c("d","pi","sigma","rho")
  coef
}


#coef.lnam - Coefficient method for lnam
coef.lnam<-function(object, ...){
   coefs<-vector()
#   cn<-vector()
   if(!is.null(object$beta)){
      coefs<-c(coefs,object$beta)
#      cn<-c(cn,names(object$beta))
   }
   if(!is.null(object$rho1)){
      coefs<-c(coefs,object$rho1)
#      cn<-c(cn,"rho1")
   }
   if(!is.null(object$rho2)){
      coefs<-c(coefs,object$rho2)
#      cn<-c(cn,"rho2")
   }
#   names(coefs)<-cn
   coefs
}


#consensus - Find a consensus structure, using one of several algorithms.  Note 
#that this is currently experimental, and that the routines are not guaranteed 
#to produce meaningful output
consensus<-function(dat,mode="digraph",diag=FALSE,method="central.graph",tol=1e-6,maxiter=1e3,verbose=TRUE,no.bias=FALSE){
   #First, prepare the data
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     stop("consensus requires graphs of identical order.")
   if(is.matrix(dat))
     m<-1
   else
   m<-dim(dat)[1]
   n<-dim(dat)[2]
   if(m==1)
     dat<-array(dat,dim=c(1,n,n))
   if(mode=="graph")
      d<-upper.tri.remove(dat)
   else
      d<-dat
   if(!diag)
      d<-diag.remove(d)
   #Now proceed by method
   #First, use the central graph if called for
   if(method=="central.graph"){
      cong<-centralgraph(d)
   #Try the iterative reweighting algorithm....
   }else if(method=="iterative.reweight"){
      cong<-centralgraph(d)
      ans<-sweep(d,c(2,3),cong,"==")
      comp<-pmax(apply(ans,1,mean,na.rm=TRUE),0.5)
      if(no.bias)
        bias<-rep(0.5,length(comp))
      else
        bias<-apply(sweep(d,c(2,3),!ans,"*"),1,mean,na.rm=TRUE)
      cdiff<-1+tol
      iter<-1
      while((cdiff>tol)&&(iter<maxiter)){
        ll1<-apply(sweep(d,1,log(comp+(1-comp)*bias),"*") + sweep(1-d,1,log((1-comp)*(1-bias)),"*"), c(2,3), sum, na.rm=TRUE)
        ll0<-apply(sweep(1-d,1,log(comp+(1-comp)*(1-bias)),"*") + sweep(d,1,log((1-comp)*bias),"*"), c(2,3), sum, na.rm=TRUE)
        cong<-ll1>ll0
        ans<-sweep(d,c(2,3),cong,"==")
        ocomp<-comp
        comp<-pmax(apply(ans,1,mean,na.rm=TRUE),0.5)
        bias<-apply(sweep(d,c(2,3),!ans,"*"),1,mean,na.rm=TRUE)
        cdiff<-sum(abs(ocomp-comp))
        iter<-iter+1
      }
      if(verbose){
        cat("Estimated competency scores:\n")
        print(comp)
        cat("Estimated bias parameters:\n")
        print(bias)
      }
   #Perform a single reweighting using mean correlation
   }else if(method=="single.reweight"){
      gc<-gcor(d)
      gc[is.na(gc)]<-0
      diag(gc)<-1
      rwv<-apply(gc,1,sum)
      rwv<-rwv/sum(rwv)
      cong<-apply(d*aperm(array(sapply(rwv,rep,n^2),dim=c(n,n,m)),c(3,2,1)),c(2,3),sum)
   #Perform a single reweighting using first component loadings
   }else if(method=="PCA.reweight"){
      gc<-gcor(d)
      gc[is.na(gc)]<-0
      diag(gc)<-1
      rwv<-abs(eigen(gc)$vector[,1])
      cong<-apply(d*aperm(array(sapply(rwv,rep,n^2),dim=c(n,n,m)), c(3,2,1)),c(2,3),sum)
   #Use the (proper) Romney-Batchelder model
   }else if(method=="romney.batchelder"){
     d<-d[!apply(is.na(d),1,all),,]              #Remove any missing informants
     if(length(dim(d))<3)
       stop("Insufficient informant information.")
     #Create the initial estimates
     drate<-apply(d,1,mean,na.rm=TRUE)
     cong<-apply(d,c(2,3),mean,na.rm=TRUE)>0.5   #Estimate graph
     s1<-mean(cong,na.rm=TRUE)
     s0<-mean(1-cong,na.rm=TRUE)
     correct<-sweep(d,c(2,3),cong,"==")          #Check for correctness
     correct<-apply(correct,1,mean,na.rm=TRUE)
     comp<-pmax(2*correct-1,0)                   #Estimate competencies
     if(no.bias)
       bias<-rep(0.5,length(comp))
     else{
       bias<-pmin(pmax((drate-s1*comp)/(1-comp),0),1)
       bias[comp==1]<-0.5
     }
     #Now, iterate until the system converges
     ocomp<-comp+tol+1
     iter<-1
     while((max(abs(ocomp-comp),na.rm=TRUE)>tol)&&(iter<maxiter)){
       ocomp<-comp
       #Estimate cong, based on the local MLE
       ll1<-apply(sweep(d,1,log(comp+(1-comp)*bias),"*") + sweep(1-d,1,log((1-comp)*(1-bias)),"*"), c(2,3), sum, na.rm=TRUE)
       ll0<-apply(sweep(1-d,1,log(comp+(1-comp)*(1-bias)),"*") + sweep(d,1,log((1-comp)*bias),"*"), c(2,3), sum, na.rm=TRUE)
       cong<-ll1>ll0
       s1<-mean(cong,na.rm=TRUE)
       s0<-mean(1-cong,na.rm=TRUE)
       #Estimate competencies, using EM
       correct<-sweep(d,c(2,3),cong,"==")         #Check for correctness
       correct<-apply(correct,1,mean,na.rm=TRUE)
       comp<-pmax((correct-s1*bias-s0*(1-bias))/(1-s1*bias-s0*(1-bias)),0)
       #Estimate bias (if no.bias==FALSE), using EM
       if(!no.bias)
         bias<-pmin(pmax((drate-s1*comp)/(1-comp),0),1)
     }
     #Possibly, dump fun messages
     if(verbose){
       cat("Estimated competency scores:\n")
       print(comp)
       cat("Estimated bias parameters:\n")
       print(bias)
     }
   #Use the Locally Aggregated Structure
   }else if(method=="LAS.intersection"){
      cong<-matrix(0,n,n)
      for(i in 1:n)
        for(j in 1:n)
          cong[i,j]<-as.numeric(d[i,i,j]&&d[j,i,j])
   }else if(method=="LAS.union"){
      cong<-matrix(0,n,n)
      for(i in 1:n)
        for(j in 1:n)
          cong[i,j]<-as.numeric(d[i,i,j]||d[j,i,j])
   }else if(method=="OR.row"){
      cong<-matrix(0,n,n)
      for(i in 1:n)
         cong[i,]<-d[i,i,]
   }else if(method=="OR.col"){
      cong<-matrix(0,n,n)
      for(i in 1:n)
         cong[,i]<-d[i,,i]
   }
   #Finish off and return the consensus graph
   if(mode=="graph")
      cong[upper.tri(cong)]<-t(cong)[upper.tri(cong)]
   if(!diag)
      diag(cong)<-0
   cong
}


#eval.edgeperturbation - Evaluate a function on a given graph with and without a
#given edge, returning the difference between the results in each case.
eval.edgeperturbation<-function(dat,i,j,FUN,...){
   #Get the function in question
   fun<-match.fun(FUN)
   #Set up the perturbation matrices
   present<-dat
   present[i,j]<-1
   absent<-dat
   absent[i,j]<-0
   #Evaluate the function across the perturbation and return the difference
   fun(present,...)-fun(absent,...)
}


#lnam - Fit a linear network autocorrelation model
#y = r1 * W1 %*% y + X %*% b + e, e = r2 * W2 %*% e + nu
#y =  (I-r1*W1)^-1%*%(X %*% b + e)
#y = (I-r1 W1)^-1 (X %*% b + (I-r2 W2)^-1 nu)
#e = (I-r2 W2)^-1 nu
#e = (I-r1 W1) y - X b
#nu = (I - r2 W2) [ (I-r1 W1) y - X b ]
#nu = (I-r2 W2) e
lnam<-function(y,x=NULL,W1=NULL,W2=NULL,theta.seed=NULL,null.model=c("meanstd","mean","std","none"),method="BFGS",control=list(),tol=1e-10){
   #Define the log-likelihood functions for each case
   agg<-function(a,w){
     m<-length(w)
     n<-dim(a)[2]
     mat<-as.double(matrix(0,n,n))
     matrix(.C("aggarray3d_R",as.double(a),as.double(w),mat=mat,as.integer(m), as.integer(n),PACKAGE="sna",NAOK=TRUE)$mat,n,n)
   }
  #Estimate covariate effects, conditional on autocorrelation parameters
  betahat<-function(y,X,W1a,W2a){
    if(nw1==0){
      if(nw2==0){
        return(qr.solve(t(X)%*%X,t(X)%*%y))
      }else{
        tXtW2aW2a<-t(X)%*%t(W2a)%*%W2a
        return(qr.solve(tXtW2aW2a%*%X,tXtW2aW2a%*%y))
      }
    }else{
      if(nw2==0){
        return(qr.solve(t(X)%*%X,t(X)%*%W1a%*%y))
      }else{
        tXtW2aW2a<-t(X)%*%t(W2a)%*%W2a
        qr.solve(tXtW2aW2a%*%X,tXtW2aW2a%*%W1a%*%y)
      }
    }
  }
  #Estimate predicted means, conditional on other effects
  muhat<-function(y,X,W1a,W2a,betahat){
    if(nx>0)
      Xb<-X%*%betahat
    else
      Xb<-0
    switch((nw1>0)+2*(nw2>0)+1,
      y-Xb,
      W1a%*%y-Xb,
      W2a%*%(y-Xb),
      W2a%*%(W1a%*%y-Xb)
    )
  }
  #Estimate innovation variance, conditional on other effects
  sigmasqhat<-function(muhat){
    t(muhat)%*%muhat/length(muhat)
  }
  #Model deviance (for use with fitting rho | beta, sigma)
  n2ll.rho<-function(rho,beta,sigmasq){
    #Prepare ll elements according to which parameters are present
    if(nw1>0){
      W1a<-diag(n)-agg(W1,rho[1:nw1])
      W1ay<-W1a%*%y
      adetW1a<-abs(det(W1a))
    }else{
      W1ay<-y
      adetW1a<-1
    }
    if(nw2>0){
      W2a<-diag(n)-agg(W2,rho[(nw1+1):(nw1+nw2)])
      tpW2a<-t(W2a)%*%W2a
      adetW2a<-abs(det(W2a))
    }else{
      tpW2a<-diag(n)
      adetW2a<-1
    }
    if(nx>0){
      Xb<-x%*%beta
    }else{
      Xb<-0
    }
    #Compute and return
    n*(log(2*pi)+log(sigmasq)) + t(W1ay-Xb)%*%tpW2a%*%(W1ay-Xb)/sigmasq - 2*(log(adetW1a)+log(adetW2a))
  }
  #Model deviance (general purpose)
  n2ll<-function(W1a,W2a,sigmasqhat){
    switch((nw1>0)+2*(nw2>0)+1,
      n*(1+log(2*pi)+log(sigmasqhat)),
      n*(1+log(2*pi)+log(sigmasqhat))-2*log(abs(det(W1a))),
      n*(1+log(2*pi)+log(sigmasqhat))-2*log(abs(det(W2a))),
      n*(1+log(2*pi)+log(sigmasqhat))- 2*(log(abs(det(W1a)))+log(abs(det(W2a))))
    )
  }
  #Conduct a single iterative refinement of a set of initial parameter estimates
  estimate<-function(parm,final=FALSE){
    #Either aggregate the weight matrices, or NULL them
    if(nw1>0)
      W1a<-diag(n)-agg(W1,parm$rho1)
    else
      W1a<-NULL
    if(nw2>0)
      W2a<-diag(n)-agg(W2,parm$rho2)
    else
      W2a<-NULL
    #If covariates were given, estimate beta | rho
    if(nx>0)
      parm$beta<-betahat(y,x,W1a,W2a)
    #Estimate sigma | beta, rho
    parm$sigmasq<-sigmasqhat(muhat(y,x,W1a,W2a,parm$beta))
    #If networks were given, (and not final) estimate rho | beta, sigma
    if(!(final||(nw1+nw2==0))){
      rho<-c(parm$rho1,parm$rho2)
      temp<-optim(rho,n2ll.rho,method=method,control=control,beta=parm$beta, sigmasq=parm$sigmasq)
      if(nw1>0)
        parm$rho1<-temp$par[1:nw1]
      if(nw2>0)
        parm$rho2<-temp$par[(nw1+1):(nw1+nw2)]
    }
    #Calculate model deviance
    parm$dev<-n2ll(W1a,W2a,parm$sigmasq)
    #Return the parameter list
    parm
  }
  #Calculate the expected Fisher information matrix for a fitted model
  infomat<-function(parm){     #Numerical version (requires numDeriv)
    require(numDeriv)
    locnll<-function(par){
      #Prepare ll elements according to which parameters are present
      if(nw1>0){
        W1a<-diag(n)-agg(W1,par[(nx+1):(nx+nw1)])
        W1ay<-W1a%*%y
        ladetW1a<-log(abs(det(W1a)))
      }else{
        W1ay<-y
        ladetW1a<-0
      }
      if(nw2>0){
        W2a<-diag(n)-agg(W2,par[(nx+nw1+1):(nx+nw1+nw2)])
        tpW2a<-t(W2a)%*%W2a
        ladetW2a<-log(abs(det(W2a)))
      }else{
        tpW2a<-diag(n)
        ladetW2a<-0
      }
      if(nx>0){
        Xb<-x%*%par[1:nx]
      }else{
        Xb<-0
      }
      #Compute and return
      n/2*(log(2*pi)+log(par[m]))+ t(W1ay-Xb)%*%tpW2a%*%(W1ay-Xb)/(2*par[m]) -ladetW1a-ladetW2a
    }
    #Return the information matrix
    hessian(locnll,c(parm$beta,parm$rho1,parm$rho2,parm$sigmasq)) 
  }
  #How many data points are there?
  n<-length(y)
  #Fix x, W1, and W2, if needed, and count predictors
  if(!is.null(x)){
    if(is.vector(x))
      x<-as.matrix(x)
    if(NROW(x)!=n)
      stop("Number of observations in x must match length of y.")
    nx<-NCOL(x)
  }else
    nx<-0
  if(!is.null(W1)){
    W1<-as.sociomatrix.sna(W1)
    if(!(is.matrix(W1)||is.array(W1)))
      stop("All networks supplied in W1 must be of identical order.")
    if(dim(W1)[2]!=n)
      stop("Order of W1 must match length of y.")
    if(length(dim(W1))==2)
      W1<-array(W1,dim=c(1,n,n))
    nw1<-dim(W1)[1]
  }else
    nw1<-0
  if(!is.null(W2)){
    W2<-as.sociomatrix.sna(W2)
    if(!(is.matrix(W2)||is.array(W2)))
      stop("All networks supplied in W2 must be of identical order.")
    if(dim(W2)[2]!=n)
      stop("Order of W2 must match length of y.")
    if(length(dim(W2))==2)
      W2<-array(W2,dim=c(1,n,n))
    nw2<-dim(W2)[1]
  }else
    nw2<-0
   #Determine the computation mode from the x,W1,W2 parameters
   comp.mode<-as.character(as.numeric(1*(nx>0)+10*(nw1>0)+100*(nw2>0)))
   if(comp.mode=="0")
      stop("At least one of x, W1, W2 must be specified.\n")
   #How many predictors?   
   m<-switch(comp.mode,
      "1"=nx+1,
      "10"=nw1+1,
      "100"=nw2+1,
      "11"=nx+nw1+1,
      "101"=nx+nw2+1,
      "110"=nw1+nw2+1,
      "111"=nx+nw1+nw2+1
   )
  #Initialize the parameter list
  parm<-list()
  if(is.null(theta.seed)){
    if(nx>0)
      parm$beta<-rep(0,nx)
    if(nw1>0)
      parm$rho1<-rep(0,nw1)
    if(nw2>0)
      parm$rho2<-rep(0,nw2)
    parm$sigmasq<-1
  }else{
    if(nx>0)
      parm$beta<-theta.seed[1:nx]
    if(nw1>0)
      parm$rho1<-theta.seed[(nx+1):(nx+nw1)]
    if(nw2>0)
      parm$rho2<-theta.seed[(nx+nw1+1):(nx+nw1+nw2)]
    parm$sigmasq<-theta.seed[nx+nw1+nw2+1]
  }
  parm$dev<-Inf
  #Fit the model
  olddev<-Inf
  while(is.na(parm$dev-olddev)||(abs(parm$dev-olddev)>tol)){
    olddev<-parm$dev
    parm<-estimate(parm,final=FALSE)
  }
  parm<-estimate(parm,final=TRUE)  #Final refinement
  #Assemble the result
  o<-list()
  o$y<-y
  o$x<-x
  o$W1<-W1
  o$W2<-W2
  o$model<-comp.mode
  o$infomat<-infomat(parm)
  o$acvm<-qr.solve(o$infomat)
  o$null.model<-match.arg(null.model)
  o$lnlik.null<-switch(match.arg(null.model),  #Fit a null model
    "meanstd"=sum(dnorm(y-mean(y),0,as.numeric(sqrt(var(y))),log=TRUE)),
    "mean"=sum(dnorm(y-mean(y),log=TRUE)),
    "std"=sum(dnorm(y,0,as.numeric(sqrt(var(y))),log=TRUE)),
    "none"=sum(dnorm(y,log=TRUE))
  )
  o$df.null.resid<-switch(match.arg(null.model),  #Find residual null df
    "meanstd"=n-2,
    "mean"=n-1,
    "std"=n-1,
    "none"=n
  )
  o$df.null<-switch(match.arg(null.model),  #Find null df
    "meanstd"=2,
    "mean"=1,
    "std"=1,
    "none"=0
  )
  o$null.param<-switch(match.arg(null.model),  #Find null params, if any
    "meanstd"=c(mean(y),sqrt(var(y))),
    "mean"=mean(y),
    "std"=sqrt(var(y)),
    "none"=NULL
  )
  o$lnlik.model<--parm$dev/2
  o$df.model<-m
  o$df.residual<-n-m
  o$df.total<-n
  o$beta<-parm$beta                       #Extract parameters
  o$rho1<-parm$rho1  
  o$rho2<-parm$rho2  
  o$sigmasq<-parm$sigmasq
  o$sigma<-o$sigmasq^0.5
  temp<-sqrt(diag(o$acvm))                #Get standard errors
  if(nx>0)
    o$beta.se<-temp[1:nx]
  if(nw1>0)
    o$rho1.se<-temp[(nx+1):(nx+nw1)]
  if(nw2>0)
    o$rho2.se<-temp[(nx+nw1+1):(nx+nw1+nw2)]
  o$sigmasq.se<-temp[m]
  o$sigma.se<-o$sigmasq.se^2/(4*o$sigmasq)  #This a delta method approximation
  if(!is.null(o$beta)){                   #Set X names
    if(!is.null(colnames(x))){
       names(o$beta)<-colnames(x)
       names(o$beta.se)<-colnames(x)
    }else{
       names(o$beta)<-paste("X",1:nx,sep="")
       names(o$beta.se)<-paste("X",1:nx,sep="")
    }
  }
  if(!is.null(o$rho1)){                     #Set W1 names
    if((!is.null(dimnames(W1)))&&(!is.null(dimnames(W1)[[1]]))){
       names(o$rho1)<-dimnames(W1)[[1]]
       names(o$rho1.se)<-dimnames(W1)[[1]]
    }else{
       names(o$rho1)<-paste("rho1",1:nw1,sep=".")
       names(o$rho1.se)<-paste("rho1",1:nw1,sep=".")
    }
  }
  if(!is.null(o$rho2)){                     #Set W2 names
    if((!is.null(dimnames(W2)))&&(!is.null(dimnames(W2)[[1]]))){
       names(o$rho2)<-dimnames(W2)[[1]]
       names(o$rho2.se)<-dimnames(W2)[[1]]
    }else{
       names(o$rho2)<-paste("rho2",1:nw2,sep=".")
       names(o$rho2.se)<-paste("rho2",1:nw2,sep=".")
    }
  }
  if(nw1>0)                               #Aggregate W1 weights
     W1ag<-agg(W1,o$rho1)
  if(nw2>0)                               #Aggregate W2 weights
     W2ag<-agg(W2,o$rho2)
  o$disturbances<-as.vector(switch(comp.mode,  #The estimated disturbances
    "1"=y-x%*%o$beta,
    "10"=(diag(n)-W1ag)%*%y,
    "100"=(diag(n)-W2ag)%*%y,
    "11"=(diag(n)-W1ag)%*%y-x%*%o$beta,
    "101"=(diag(n)-W2ag)%*%(y-x%*%o$beta),
    "110"=(diag(n)-W2ag)%*%((diag(n)-W1ag)%*%y),
    "111"=(diag(n)-W2ag)%*%((diag(n)-W1ag)%*%y-x%*%o$beta)
  ))
  o$fitted.values<-as.vector(switch(comp.mode,  #Compute the fitted values
    "1"=x%*%o$beta,
    "10"=rep(0,n),
    "100"=rep(0,n),
    "11"=qr.solve(diag(n)-W1ag,x%*%o$beta),
    "101"=x%*%o$beta,
    "110"=rep(0,n),
    "111"=qr.solve(diag(n)-W1ag,x%*%o$beta)
  ))
  o$residuals<-as.vector(y-o$fitted.values)
  o$call<-match.call()
  class(o)<-c("lnam")
  o
}


#nacf - Network autocorrelation function
nacf<-function(net,y,lag.max=NULL,type=c("correlation","covariance","moran","geary"),neighborhood.type=c("in","out","total"),partial.neighborhood=TRUE,mode="digraph",diag=FALSE,thresh=0,demean=TRUE){
  #Pre-process the raw input
  net<-as.sociomatrix.sna(net)
  if(is.list(net))
    return(lapply(net,nacf,y=y,lag.max=lag.max, neighborhood.type=neighborhood.type,partial.neighborhood=partial.neighborhood,mode=mode,diag=diag,thresh=thresh,demean=demean))
  else if(length(dim(net))>2)
    return(apply(net,1,nacf,y=y,lag.max=lag.max, neighborhood.type=neighborhood.type,partial.neighborhood=partial.neighborhood,mode=mode,diag=diag,thresh=thresh,demean=demean))
  #End pre-processing
  if(length(y)!=NROW(net))
    stop("Network size must match covariate length in nacf.")
  #Process y
  if(demean||(match.arg(type)=="moran"))
    y<-y-mean(y)
  vary<-var(y)
  #Determine maximum lag, if needed
  if(is.null(lag.max))
    lag.max<-NROW(net)-1
  #Get the appropriate neighborhood graphs for dat
  neigh<-neighborhood(net,order=lag.max,neighborhood.type=neighborhood.type, mode=mode,diag=diag,thresh=thresh,return.all=TRUE,partial=partial.neighborhood)
  #Form the coefficients
  v<-switch(match.arg(type),
    "covariance"=t(y)%*%y/NROW(net),
    "correlation"=1,
    "moran"=1,
    "geary"=0,
  )
  for(i in 1:lag.max){
    ec<-sum(neigh[i,,])
    if(ec>0){
      v[i+1]<-switch(match.arg(type),
        "covariance"=(t(y)%*%neigh[i,,]%*%y)/ec,
        "correlation"=((t(y)%*%neigh[i,,]%*%y)/ec)/vary,
        "moran"=NROW(net)/ec*sum((y%o%y)*neigh[i,,])/sum(y^2),
        "geary"=(NROW(net)-1)/(2*ec)*sum(neigh[i,,]*outer(y,y,"-")^2)/ sum((y-mean(y))^2),
      )
    }else
      v[i+1]<-0
  }
  names(v)<-0:lag.max
  #Return the results
  v
}


#netcancor - Canonical correlations for network variables. 
netcancor<-function(y,x,mode="digraph",diag=FALSE,nullhyp="cugtie",reps=1000){
   y<-as.sociomatrix.sna(y)
   x<-as.sociomatrix.sna(x)
   if(is.list(x)|is.list(y))
     stop("netcancor requires graphs of identical order.")
   if(length(dim(y))>2){
      iy<-matrix(nrow=dim(y)[1],ncol=dim(y)[2]*dim(y)[3])
   }else{
      iy<-matrix(nrow=1,ncol=dim(y)[1]*dim(y)[2])
      temp<-y
      y<-array(dim=c(1,dim(temp)[1],dim(temp)[2]))
      y[1,,]<-temp
   }
   if(length(dim(x))>2){
      ix<-matrix(nrow=dim(x)[1],ncol=dim(x)[2]*dim(x)[3])
   }else{
      ix<-matrix(nrow=1,ncol=dim(x)[1]*dim(x)[2])
      temp<-x
      x<-array(dim=c(1,dim(temp)[1],dim(temp)[2]))
      x[1,,]<-temp
   }
   my<-dim(y)[1]
   mx<-dim(x)[1]
   n<-dim(y)[2]
   out<-list()
   out$xdist<-array(dim=c(reps,mx,mx))
   out$ydist<-array(dim=c(reps,my,my))
   #Convert the response first.
   for(i in 1:my){
      d<-y[i,,]
      #if(!diag){
      #   diag(d)<-NA
      #}
      #if(mode!="digraph")
      #   d[lower.tri(d)]<-NA
      iy[i,]<-as.vector(d)
   }
   #Now for the independent variables.
   for(i in 1:mx){
      d<-x[i,,]
      #if(!diag){
      #   diag(d)<-NA
      #}
      #if(mode!="digraph")
      #   d[lower.tri(d)]<-NA
      ix[i,]<-as.vector(d)
   }   
   #Run the initial model fit
   nc<-cancor(t(ix),t(iy))  #Had to take out na.action=na.omit, since it's not supported
   #Now, repeat the whole thing an ungodly number of times.
   out$cdist<-array(dim=c(reps,length(nc$cor)))
   for(i in 1:reps){
      #Clear out the internal structures
      iy<-matrix(nrow=dim(y)[1],ncol=dim(y)[2]*dim(y)[3])
      ix<-matrix(nrow=dim(x)[1],ncol=dim(x)[2]*dim(x)[3])
      #Convert (and mutate) the response first.
      for(j in 1:my){
         d<-switch(nullhyp,
            qap = rmperm(y[j,,]),
            cug = rgraph(n,1,mode=mode,diag=diag),
            cugden = rgraph(n,1,tprob=gden(y[j,,],mode=mode,diag=diag),mode=mode,diag=diag),
            cugtie = rgraph(n,1,mode=mode,diag=diag,tielist=y[j,,])
         )
         #if(!diag){
         #   diag(d)<-NA
         #}
         #if(mode!="digraph")
         #   d[lower.tri(d)]<-NA
         iy[j,]<-as.vector(d)
      }
      #Now for the independent variables.
      for(j in 1:mx){
         d<-switch(nullhyp,
            qap = rmperm(x[j,,]),
            cug = rgraph(n,1,mode=mode,diag=diag),
            cugden = rgraph(n,1,tprob=gden(x[j,,],mode=mode,diag=diag),mode=mode,diag=diag),
            cugtie = rgraph(n,1,mode=mode,diag=diag,tielist=x[j,,])
         )
         #if(!diag){
         #   diag(d)<-NA
         #}
         #if(mode!="digraph")
         #   d[lower.tri(d)]<-NA
         ix[j,]<-as.vector(d)
      }   
      #Finally, fit the test model
      tc<-cancor(t(ix),t(iy))         #Had to take out na.action=na.omit, since it's not supported
      #Gather the coefficients for use later...
      out$cdist[i,]<-tc$cor
      out$xdist[i,,]<-tc$xcoef
      out$ydist[i,,]<-tc$ycoef
   }
   #Find the p-values for our monte carlo null hypothesis tests
   out$cor<-nc$cor
   out$xcoef<-nc$xcoef
   out$ycoef<-nc$ycoef
   out$cpgreq<-vector(length=length(nc$cor))
   out$cpleeq<-vector(length=length(nc$cor))
   for(i in 1:length(nc$cor)){
      out$cpgreq[i]<-mean(out$cdist[,i]>=out$cor[i],na.rm=TRUE)
      out$cpleeq[i]<-mean(out$cdist[,i]<=out$cor[i],na.rm=TRUE)
   }
   out$xpgreq<-matrix(ncol=mx,nrow=mx)
   out$xpleeq<-matrix(ncol=mx,nrow=mx)
   for(i in 1:mx){
      for(j in 1:mx){ 
         out$xpgreq[i,j]<-mean(out$xdist[,i,j]>=out$xcoef[i,j],na.rm=TRUE)
         out$xpleeq[i,j]<-mean(out$xdist[,i,j]<=out$xcoef[i,j],na.rm=TRUE)
      }
   }
   out$ypgreq<-matrix(ncol=my,nrow=my)
   out$ypleeq<-matrix(ncol=my,nrow=my)
   for(i in 1:my){
      for(j in 1:my){ 
         out$ypgreq[i,j]<-mean(out$ydist[,i,j]>=out$ycoef[i,j],na.rm=TRUE)
         out$ypleeq[i,j]<-mean(out$ydist[,i,j]<=out$ycoef[i,j],na.rm=TRUE)
      }
   }
   #Having completed the model fit and MC tests, we gather useful information for
   #the end user.  This is a combination of cancor output and our own stuff.
   out$cnames<-as.vector(paste("cor",1:min(mx,my),sep=""))
   out$xnames<-as.vector(paste("x",1:mx,sep=""))
   out$ynames<-as.vector(paste("y",1:my,sep=""))
   out$xcenter<-nc$xcenter
   out$ycenter<-nc$ycenter
   out$nullhyp<-nullhyp
   class(out)<-c("netcancor")
   out
}


#netlm - OLS network regrssion routine (w/many null hypotheses)
#  #QAP semi-partialling "plus" (Dekker et al.)
#  #For Y ~ b0 + b1 X1 + b2 X2 + ... + bp Xp
#  #for(i in 1:p)
#  #  Fit Xi ~ b0* + b1* X1 + ... + bp* Xp (omit Xi)
#  #  Let ei = resid of above lm
#  #  for(j in 1:reps)
#  #    eij = rmperm (ei)
#  #    Fit Y ~ b0** + b1** X1 + ... + bi** eij + ... + bp** Xp
#  #Use resulting permutation distributions to test coefficients
netlm<-function(y,x,intercept=TRUE,mode="digraph",diag=FALSE,nullhyp=c("qap", "qapspp","qapy","qapx","qapallx","cugtie","cugden","cuguman","classical"),test.statistic=c("t-value","beta"),tol=1e-7, reps=1000){
  #Define an internal routine to perform a QR reduction and get t-values
  gettval<-function(x,y,tol){
    xqr<-qr(x,tol=tol)
    coef<-qr.coef(xqr,y)
    resid<-qr.resid(xqr,y)
    rank<-xqr$rank
    n<-length(y)
    rdf<-n-rank
    resvar<-sum(resid^2)/rdf
    cvm<-chol2inv(xqr$qr)
    se<-sqrt(diag(cvm)*resvar)
    coef/se
  }
  #Define an internal routine to quickly fit linear models to graphs
  gfit<-function(glist,mode,diag,tol,rety,tstat){
    y<-gvectorize(glist[[1]],mode=mode,diag=diag,censor.as.na=TRUE)
    x<-vector()
    for(i in 2:length(glist))
      x<-cbind(x,gvectorize(glist[[i]],mode=mode,diag=diag,censor.as.na=TRUE))
    if(!is.matrix(x))
      x<-matrix(x,ncol=1)
    mis<-is.na(y)|apply(is.na(x),1,any)
    if(!rety){
      if(tstat=="beta")
        qr.solve(x[!mis,],y[!mis],tol=tol)
      else if(tstat=="t-value"){
        gettval(x[!mis,],y[!mis],tol=tol)
      }
    }else{
      list(qr(x[!mis,],tol=tol),y[!mis])
    }
  }
  #Get the data in order
  y<-as.sociomatrix.sna(y)
  x<-as.sociomatrix.sna(x)
  if(is.list(y)||((length(dim(y))>2)&&(dim(y)[1]>1))) 
    stop("y must be a single graph in netlm.")
  if(length(dim(y))>2)
    y<-y[1,,]
  if(is.list(x)||(dim(x)[2]!=dim(y)[2]))
    stop("Homogeneous graph orders required in netlm.")
  nx<-stackcount(x)+intercept    #Get number of predictors
  n<-dim(y)[2]                   #Get graph order
  g<-list(y)                     #Put graphs into a list
  if(intercept)
    g[[2]]<-matrix(1,n,n)
  if(nx-intercept==1)
    g[[2+intercept]]<-x
  else
    for(i in 1:(nx-intercept))
      g[[i+1+intercept]]<-x[i,,]
  if(any(sapply(lapply(g,is.na),any)))
    warning("Missing data supplied to netlm; this may pose problems for certain null hypotheses.  Hope you know what you're doing....")
  #Fit the initial baseline model
  fit.base<-gfit(g,mode=mode,diag=diag,tol=tol,rety=TRUE)
  fit<-list()  #Initialize output
  fit$coefficients<-qr.coef(fit.base[[1]],fit.base[[2]])
  fit$fitted.values<-qr.fitted(fit.base[[1]],fit.base[[2]])
  fit$residuals<-qr.resid(fit.base[[1]],fit.base[[2]])
  fit$qr<-fit.base[[1]]
  fit$rank<-fit.base[[1]]$rank
  fit$n<-length(fit.base[[2]])
  fit$df.residual<-fit$n-fit$rank
  tstat<-match.arg(test.statistic)
  if(tstat=="beta")
    fit$tstat<-fit$coefficients
  else if(tstat=="t-value")
    fit$tstat<-fit$coefficients/ sqrt(diag(chol2inv(fit$qr$qr))*sum(fit$residuals^2)/(fit$n-fit$rank))
  #Proceed based on selected null hypothesis
  nullhyp<-match.arg(nullhyp)
  if((nullhyp%in%c("qap","qapspp"))&&(nx==1))  #No partialling w/one predictor
    nullhyp<-"qapy"
  if(nullhyp=="classical"){
    resvar<-sum(fit$residuals^2)/fit$df.residual
    cvm<-chol2inv(fit$qr$qr)
    se<-sqrt(diag(cvm)*resvar)
    tval<-fit$coefficients/se
    #Prepare output
    fit$dist<-NULL
    fit$pleeq<-pt(tval,fit$df.residual)
    fit$pgreq<-pt(tval,fit$df.residual,lower.tail=FALSE)
    fit$pgreqabs<-2*pt(abs(tval),fit$df.residual,lower.tail=FALSE)
  }else if(nullhyp%in%c("cugtie","cugden","cuguman")){
    #Generate replicates for each predictor
    repdist<-matrix(0,reps,nx)
    for(i in 1:nx){
      gr<-g
      for(j in 1:reps){
        #Modify the focal x
        gr[[i+1]]<-switch(nullhyp,
          cugtie<-rgraph(n,mode=mode,diag=diag,replace=FALSE,tielist=g[[i+1]]),
          cugden<-rgraph(n,tprob=gden(g[[i+1]],mode=mode,diag=diag),mode=mode, diag=diag),
          cuguman<-(function(dc,n){rguman(1,n,mut=x[1],asym=x[2],null=x[3], method="exact")})(dyad.census(g[[i+1]]),n)
        )
        #Fit model with modified x
        repdist[j,i]<-gfit(gr,mode=mode,diag=diag,tol=tol,rety=FALSE, tstat=tstat)[i]
      }
    }
    #Prepare output
    fit$dist<-repdist
    fit$pleeq<-apply(sweep(fit$dist,2,fit$tstat,"<="),2,mean)
    fit$pgreq<-apply(sweep(fit$dist,2,fit$tstat,">="),2,mean)
    fit$pgreqabs<-apply(sweep(abs(fit$dist),2,abs(fit$tstat),">="),2, mean)
  }else if(nullhyp=="qapy"){
    #Generate replicates for each predictor
    repdist<-matrix(0,reps,nx)
    gr<-g
    for(i in 1:reps){
      gr[[1]]<-rmperm(g[[1]])  #Permute y
      #Fit the model under replication
      repdist[i,]<-gfit(gr,mode=mode,diag=diag,tol=tol,rety=FALSE,tstat=tstat)
    }
    #Prepare output
    fit$dist<-repdist
    fit$pleeq<-apply(sweep(fit$dist,2,fit$tstat,"<="),2,mean)
    fit$pgreq<-apply(sweep(fit$dist,2,fit$tstat,">="),2,mean)
    fit$pgreqabs<-apply(sweep(abs(fit$dist),2,abs(fit$tstat),">="),2, mean)
  }else if(nullhyp=="qapx"){
    #Generate replicates for each predictor
    repdist<-matrix(0,reps,nx)
    for(i in 1:nx){
      gr<-g
      for(j in 1:reps){
        gr[[i+1]]<-rmperm(gr[[i+1]]) #Modify the focal x
        #Fit model with modified x
        repdist[j,i]<-gfit(gr,mode=mode,diag=diag,tol=tol,rety=FALSE, tstat=tstat)[i]
      }
    }
    #Prepare output
    fit$dist<-repdist
    fit$pleeq<-apply(sweep(fit$dist,2,fit$tstat,"<="),2,mean)
    fit$pgreq<-apply(sweep(fit$dist,2,fit$tstat,">="),2,mean)
    fit$pgreqabs<-apply(sweep(abs(fit$dist),2,abs(fit$tstat),">="),2, mean)
  }else if(nullhyp=="qapallx"){
    #Generate replicates for each predictor
    repdist<-matrix(0,reps,nx)
    gr<-g
    for(i in 1:reps){
      for(j in 1:nx)
        gr[[1+j]]<-rmperm(g[[1+j]])  #Permute each x
      #Fit the model under replication
      repdist[i,]<-gfit(gr,mode=mode,diag=diag,tol=tol,rety=FALSE, tstat=tstat)
    }
    #Prepare output
    fit$dist<-repdist
    fit$pleeq<-apply(sweep(fit$dist,2,fit$tstat,"<="),2,mean)
    fit$pgreq<-apply(sweep(fit$dist,2,fit$tstat,">="),2,mean)
    fit$pgreqabs<-apply(sweep(abs(fit$dist),2,abs(fit$tstat),">="),2, mean)
  }else if((nullhyp=="qap")||(nullhyp=="qapspp")){
    xsel<-matrix(TRUE,n,n)
    if(!diag)
      diag(xsel)<-FALSE
    if(mode=="graph")
      xsel[upper.tri(xsel)]<-FALSE
    #Generate replicates for each predictor
    repdist<-matrix(0,reps,nx)
    for(i in 1:nx){
      #Regress x_i on other x's
      xfit<-gfit(g[1+c(i,(1:nx)[-i])],mode=mode,diag=diag,tol=tol,rety=TRUE, tstat=tstat)
      xres<-g[[1+i]]
      xres[xsel]<-qr.resid(xfit[[1]],xfit[[2]])  #Get residuals of x_i
      if(mode=="graph")
        xres[upper.tri(xres)]<-t(xres)[upper.tri(xres)]
      #Draw replicate coefs using permuted x residuals
      for(j in 1:reps)
        repdist[j,i]<-gfit(c(g[-(1+i)],list(rmperm(xres))),mode=mode,diag=diag, tol=tol,rety=FALSE,tstat=tstat)[nx]
    }
    #Prepare output
    fit$dist<-repdist
    fit$pleeq<-apply(sweep(fit$dist,2,fit$tstat,"<="),2,mean)
    fit$pgreq<-apply(sweep(fit$dist,2,fit$tstat,">="),2,mean)
    fit$pgreqabs<-apply(sweep(abs(fit$dist),2,abs(fit$tstat),">="),2, mean)
  }
  #Finalize the results
  fit$nullhyp<-nullhyp
  fit$names<-paste("x",1:(nx-intercept),sep="")
  if(intercept)
    fit$names<-c("(intercept)",fit$names)
  fit$intercept<-intercept
  class(fit)<-"netlm"
  #Return the result
  fit
}


#netlogit - God help me, it's a network regression routine using a 
#binomial/logit GLM.  It's also frighteningly slow, since it's essentially a 
#front end to the builtin GLM routine with a bunch of network hypothesis testing
#stuff thrown in for good measure.
netlogit<-function(y,x,intercept=TRUE,mode="digraph",diag=FALSE,nullhyp=c("qap", "qapspp","qapy","qapx","qapallx","cugtie","cugden","cuguman","classical"), tol=1e-7,reps=1000){
  #Define an internal routine to quickly fit logit models to graphs
  gfit<-function(glist,mode,diag){
    y<-gvectorize(glist[[1]],mode=mode,diag=diag,censor.as.na=TRUE)
    x<-vector()
    for(i in 2:length(glist))
      x<-cbind(x,gvectorize(glist[[i]],mode=mode,diag=diag,censor.as.na=TRUE))
    if(!is.matrix(x))
      x<-matrix(x,ncol=1)
    mis<-is.na(y)|apply(is.na(x),1,any)
    glm.fit(x[!mis,],y[!mis],family=binomial(),intercept=FALSE)
  }
  #Repeat the above, for strictly linear models
  gfitlm<-function(glist,mode,diag,tol){
    y<-gvectorize(glist[[1]],mode=mode,diag=diag,censor.as.na=TRUE)
    x<-vector()
    for(i in 2:length(glist))
      x<-cbind(x,gvectorize(glist[[i]],mode=mode,diag=diag,censor.as.na=TRUE))
    if(!is.matrix(x))
      x<-matrix(x,ncol=1)
    mis<-is.na(y)|apply(is.na(x),1,any)
    list(qr(x[!mis,],tol=tol),y[!mis])
  }
  #Get the data in order
  y<-as.sociomatrix.sna(y)
  x<-as.sociomatrix.sna(x)
  if(is.list(y)||((length(dim(y))>2)&&(dim(y)[1]>1))) 
    stop("y must be a single graph in netlogit.")
  if(length(dim(y))>2)
    y<-y[1,,]
  if(is.list(x)||(dim(x)[2]!=dim(y)[2]))
    stop("Homogeneous graph orders required in netlogit.")
  nx<-stackcount(x)+intercept    #Get number of predictors
  n<-dim(y)[2]                   #Get graph order
  g<-list(y)                     #Put graphs into a list
  if(intercept)
    g[[2]]<-matrix(1,n,n)
  if(nx-intercept==1)
    g[[2+intercept]]<-x
  else
    for(i in 1:(nx-intercept))
      g[[i+1+intercept]]<-x[i,,]
  if(any(sapply(lapply(g,is.na),any)))
    warning("Missing data supplied to netlogit; this may pose problems for certain null hypotheses.  Hope you know what you're doing....")
  #Fit the initial baseline model
  fit.base<-gfit(g,mode=mode,diag=diag)
  fit<-list()  #Initialize output
  fit$coefficients<-fit.base$coefficients
  fit$fitted.values<-fit.base$fitted.values
  fit$residuals<-fit.base$residuals
  fit$linear.predictors<-fit.base$linear.predictors
  fit$n<-length(fit.base$y)
  fit$df.model<-fit.base$rank
  fit$df.residual<-fit.base$df.residual
  fit$deviance<-fit.base$deviance
  fit$null.deviance<-fit.base$null.deviance
  fit$df.null<-fit.base$df.null
  fit$aic<-fit.base$aic
  fit$bic<-fit$deviance+fit$df.model*log(fit$n)
  fit$qr<-fit.base$qr
  fit$ctable<-table(as.numeric(fit$fitted.values>=0.5),fit.base$y, dnn=c("Predicted","Actual"))  #Get the contingency table 
  if(NROW(fit$ctable)==1){
    if(rownames(fit$ctable)=="0")
      fit$ctable<-rbind(fit$ctable,c(0,0))
    else
      fit$ctable<-rbind(c(0,0),fit$ctable)
    rownames(fit$ctable)<-c("0","1")
  }

  #Proceed based on selected null hypothesis
  nullhyp<-match.arg(nullhyp)
  if((nullhyp%in%c("qap","qapspp"))&&(nx==1))  #No partialling w/one predictor
    nullhyp<-"qapy"
  if(nullhyp=="classical"){
    cvm<-chol2inv(fit$qr$qr)
    se<-sqrt(diag(cvm))
    tval<-fit$coefficients/se
    #Prepare output
    fit$dist<-NULL
    fit$pleeq<-pt(tval,fit$df.residual)
    fit$pgreq<-pt(tval,fit$df.residual,lower.tail=FALSE)
    fit$pgreqabs<-2*pt(abs(tval),fit$df.residual,lower.tail=FALSE)
  }else if(nullhyp%in%c("cugtie","cugden","cuguman")){
    #Generate replicates for each predictor
    repdist<-matrix(0,reps,nx)
    for(i in 1:nx){
      gr<-g
      for(j in 1:reps){
        #Modify the focal x
        gr[[i+1]]<-switch(nullhyp,
          cugtie<-rgraph(n,mode=mode,diag=diag,replace=FALSE,tielist=g[[i+1]]),
          cugden<-rgraph(n,tprob=gden(g[[i+1]],mode=mode,diag=diag),mode=mode, diag=diag),
          cuguman<-(function(dc,n){rguman(1,n,mut=x[1],asym=x[2],null=x[3], method="exact")})(dyad.census(g[[i+1]]),n)
        )
        #Fit model with modified x
        repdist[j,i]<-gfit(gr,mode=mode,diag=diag)$coef[i]
      }
    }
    #Prepare output
    fit$dist<-repdist
    fit$pleeq<-apply(sweep(fit$dist,2,fit$coefficients,"<="),2,mean)
    fit$pgreq<-apply(sweep(fit$dist,2,fit$coefficients,">="),2,mean)
    fit$pgreqabs<-apply(sweep(abs(fit$dist),2,abs(fit$coefficients),">="),2, mean)
  }else if(nullhyp=="qapy"){
    #Generate replicates for each predictor
    repdist<-matrix(0,reps,nx)
    gr<-g
    for(i in 1:reps){
      gr[[1]]<-rmperm(g[[1]])  #Permute y
      #Fit the model under replication
      repdist[i,]<-gfit(gr,mode=mode,diag=diag)$coef
    }
    #Prepare output
    fit$dist<-repdist
    fit$pleeq<-apply(sweep(fit$dist,2,fit$coefficients,"<="),2,mean)
    fit$pgreq<-apply(sweep(fit$dist,2,fit$coefficients,">="),2,mean)
    fit$pgreqabs<-apply(sweep(abs(fit$dist),2,abs(fit$coefficients),">="),2, mean)
  }else if(nullhyp=="qapx"){
    #Generate replicates for each predictor
    repdist<-matrix(0,reps,nx)
    for(i in 1:nx){
      gr<-g
      for(j in 1:reps){
        gr[[i+1]]<-rmperm(gr[[i+1]]) #Modify the focal x
        #Fit model with modified x
        repdist[j,i]<-gfit(gr,mode=mode,diag=diag)$coef[i]
      }
    }
    #Prepare output
    fit$dist<-repdist
    fit$pleeq<-apply(sweep(fit$dist,2,fit$coefficients,"<="),2,mean)
    fit$pgreq<-apply(sweep(fit$dist,2,fit$coefficients,">="),2,mean)
    fit$pgreqabs<-apply(sweep(abs(fit$dist),2,abs(fit$coefficients),">="),2, mean)
  }else if(nullhyp=="qapallx"){
    #Generate replicates for each predictor
    repdist<-matrix(0,reps,nx)
    gr<-g
    for(i in 1:reps){
      for(j in 1:nx)
        gr[[1+j]]<-rmperm(g[[1+j]])  #Permute each x
      #Fit the model under replication
      repdist[i,]<-gfit(gr,mode=mode,diag=diag)$coef
    }
    #Prepare output
    fit$dist<-repdist
    fit$pleeq<-apply(sweep(fit$dist,2,fit$coefficients,"<="),2,mean)
    fit$pgreq<-apply(sweep(fit$dist,2,fit$coefficients,">="),2,mean)
    fit$pgreqabs<-apply(sweep(abs(fit$dist),2,abs(fit$coefficients),">="),2, mean)
  }else if((nullhyp=="qap")||(nullhyp=="qapspp")){
    xsel<-matrix(TRUE,n,n)
    if(!diag)
      diag(xsel)<-FALSE
    if(mode=="graph")
      xsel[upper.tri(xsel)]<-FALSE
    #Generate replicates for each predictor
    repdist<-matrix(0,reps,nx)
    for(i in 1:nx){
      #Regress x_i on other x's
      xfit<-gfitlm(g[1+c(i,(1:nx)[-i])],mode=mode,diag=diag,tol=tol)
      xres<-g[[1+i]]
      xres[xsel]<-qr.resid(xfit[[1]],xfit[[2]])  #Get residuals of x_i
      if(mode=="graph")
        xres[upper.tri(xres)]<-t(xres)[upper.tri(xres)]
      #Draw replicate coefs using permuted x residuals
      for(j in 1:reps)
        repdist[j,i]<-gfit(c(g[-(1+i)],list(rmperm(xres))),mode=mode, diag=diag)$coef[nx]
    }
    #Prepare output
    fit$dist<-repdist
    fit$pleeq<-apply(sweep(fit$dist,2,fit$coefficients,"<="),2,mean)
    fit$pgreq<-apply(sweep(fit$dist,2,fit$coefficients,">="),2,mean)
    fit$pgreqabs<-apply(sweep(abs(fit$dist),2,abs(fit$coefficients),">="),2, mean)
  }
  #Finalize the results
  fit$nullhyp<-nullhyp
  fit$names<-paste("x",1:(nx-intercept),sep="")
  if(intercept)
    fit$names<-c("(intercept)",fit$names)
  fit$intercept<-intercept
  class(fit)<-"netlogit"
  #Return the result
  fit
}


#npostpred - Take posterior predictive draws for functions of networks.
npostpred<-function(b,FUN,...){
   #Find the desired function
   fun<-match.fun(FUN)
   #Take the draws
   out<-apply(b$net,1,fun,...)
   out
}


#plot.bbnam - Plot method for bbnam
plot.bbnam<-function(x,mode="density",intlines=TRUE,...){
   UseMethod("plot",x)
}


#plot.bbnam.actor - Plot method for bbnam.actor
plot.bbnam.actor<-function(x,mode="density",intlines=TRUE,...){
   #Get the initial graphical settings, so we can restore them later
   oldpar<-par(no.readonly=TRUE)
   #Change plotting params
   par(ask=dev.interactive())
   #Initial plot: global error distribution
   par(mfrow=c(2,1))
   if(mode=="density"){   #Approximate the pdf using kernel density estimation
      #Plot marginal population (i.e. across actors) density of p(false negative)
      plot(density(x$em),main=substitute(paste("Estimated Marginal Population Density of ",{e^{"-"}},", ",draws," Draws"),list(draws=x$draws)),xlab=expression({e^{"-"}}),xlim=c(0,1),...)
      #Plot interval lines if required.
      if(intlines)
         abline(v=quantile(x$em,c(0.05,0.5,0.95)),lty=c(3,2,3))
      #Plot marginal population (i.e. across actors) density of p(false positive)
      plot(density(x$ep),main=substitute(paste("Estimated Marginal Population Density of ",{e^{"+"}},", ",draws," Draws"),list(draws=x$draws)),xlab=expression({e^{"+"}}),xlim=c(0,1),...)
      #Plot interval lines if required.
      if(intlines)
         abline(v=quantile(x$ep,c(0.05,0.5,0.95)),lty=c(3,2,3))
   }else{     #Use histograms to plot the estimated density
      #Plot marginal population (i.e. across actors) density of p(false negative)
      hist(x$em,main=substitute(paste("Histogram of ",{e^{"-"}},", ",draws," Draws"),list(draws=x$draws)),xlab=expression({e^{"-"}}),xlim=c(0,1),...)
      #Plot interval lines if required.
      if(intlines)
         abline(v=quantile(x$em,c(0.05,0.5,0.95)),lty=c(3,2,3))
      #Plot marginal population (i.e. across actors) density of p(false positive)
      hist(x$ep,main=substitute(paste("Histogram of ",{e^{"+"}},", ",draws," Draws"),list(draws=x$draws)),xlab=expression({e^{"+"}}),xlim=c(0,1),...)
      #Plot interval lines if required.
      if(intlines)
         abline(v=quantile(x$ep,c(0.05,0.5,0.95)),lty=c(3,2,3))
   }
   #Plot e- next
   par(mfrow=c(floor(sqrt(x$nobservers)),ceiling(sqrt(x$nobservers))))
   for(i in 1:x$nobservers){
      if(mode=="density"){
         plot(density(x$em[,i]),main=substitute({e^{"-"}}[it],list(it=i)), xlab=substitute({e^{"-"}}[it],list(it=i)),xlim=c(0,1),...)
         #Plot interval lines if required.
         if(intlines)
            abline(v=quantile(x$em[,i],c(0.05,0.5,0.95)),lty=c(3,2,3))
      }else{
         hist(x$em[,i],main=substitute({e^{"-"}}[it],list(it=i)), xlab=substitute({e^{"-"}}[it],list(it=i)),xlim=c(0,1),...)
         #Plot interval lines if required.
         if(intlines)
            abline(v=quantile(x$em[,i],c(0.05,0.5,0.95)),lty=c(3,2,3))
      }
   }
   #Now plot e+
   par(mfrow=c(floor(sqrt(x$nobservers)),ceiling(sqrt(x$nobservers))))
   for(i in 1:x$nobservers){
      if(mode=="density"){
         plot(density(x$ep[,i]),main=substitute({e^{"+"}}[it],list(it=i)), xlab=substitute({e^{"+"}}[it],list(it=i)),xlim=c(0,1),...)
         #Plot interval lines if required.
         if(intlines)
            abline(v=quantile(x$ep[,i],c(0.05,0.5,0.95)),lty=c(3,2,3))
      }else{
         hist(x$ep[,i],main=substitute({e^{"+"}}[it],list(it=i)), xlab=substitute({e^{"+"}}[it],list(it=i)),xlim=c(0,1),...)
         #Plot interval lines if required.
         if(intlines)
            abline(v=quantile(x$ep[,i],c(0.05,0.5,0.95)),lty=c(3,2,3))
      }
   }
   #Finally, try to plot histograms of tie probabilities
   par(mfrow=c(1,1))
   plot.sociomatrix(apply(x$net,c(2,3),mean),labels=list(x$anames,x$anames),main="Marginal Posterior Tie Probability Distribution")
   #Clean up
   par(oldpar)
}


#plot.bbnam.fixed - Plot method for bbnam.fixed
plot.bbnam.fixed<-function(x,mode="density",intlines=TRUE,...){
   #Get the initial graphical settings, so we can restore them later
   oldpar<-par()
   #Perform matrix plot of tie probabilities
   par(mfrow=c(1,1))
   plot.sociomatrix(apply(x$net,c(2,3),mean),labels=list(x$anames,x$anames),main="Marginal Posterior Tie Probability Distribution")
   #Clean up
   par(oldpar)
}


#plot.bbnam.pooled - Plot method for bbnam.pooled
plot.bbnam.pooled<-function(x,mode="density",intlines=TRUE,...){
   #Get the initial graphical settings, so we can restore them later
   oldpar<-par()
   #Change plotting params
   par(ask=dev.interactive())
   #Initial plot: pooled error distribution
   par(mfrow=c(2,1))
   if(mode=="density"){   #Approximate the pdf using kernel density estimation
      #Plot marginal population (i.e. across actors) density of p(false negative)
      plot(density(x$em),main=substitute(paste("Estimated Marginal Posterior Density of ",{e^{"-"}},", ",draws," Draws"),list(draws=x$draws)),xlab=expression({e^{"-"}}),xlim=c(0,1),...)
      #Plot interval lines if required.
      if(intlines)
         abline(v=quantile(x$em,c(0.05,0.5,0.95)),lty=c(3,2,3))
      #Plot marginal population (i.e. across actors) density of p(false positive)
      plot(density(x$ep),main=substitute(paste("Estimated Marginal Posterior Density of ",{e^{"+"}},", ",draws," Draws"),list(draws=x$draws)),xlab=expression({e^{"+"}}),xlim=c(0,1),...)
      #Plot interval lines if required.
      if(intlines)
         abline(v=quantile(x$ep,c(0.05,0.5,0.95)),lty=c(3,2,3))
   }else{     #Use histograms to plot the estimated density
      #Plot marginal population (i.e. across actors) density of p(false negative)
      hist(x$em,main=substitute(paste("Histogram of ",{e^{"-"}},", ",draws," Draws"),list(draws=x$draws)),xlab=expression({e^{"-"}}),xlim=c(0,1),...)
      #Plot interval lines if required.
      if(intlines)
         abline(v=quantile(x$em,c(0.05,0.5,0.95)),lty=c(3,2,3))
      #Plot marginal population (i.e. across actors) density of p(false positive)
      hist(x$ep,main=substitute(paste("Histogram of ",{e^{"+"}},", ",draws," Draws"),list(draws=x$draws)),xlab=expression({e^{"+"}}),xlim=c(0,1),...)
      #Plot interval lines if required.
      if(intlines)
         abline(v=quantile(x$ep,c(0.05,0.5,0.95)),lty=c(3,2,3))
   }
   #Finally, try to plot histograms of tie probabilities
   par(mfrow=c(1,1))
   plot.sociomatrix(apply(x$net,c(2,3),mean),labels=list(x$anames,x$anames),main="Marginal Posterior Tie Probability Distribution")
   #Clean up
   par(oldpar)
}


#plot.bn - Plot method for bn
plot.bn<-function(x,...){
  op<-par(no.readonly=TRUE) #Store old plotting params
  on.exit(par(op))          #Reset when finished
  par(mfrow=c(2,2))
  #Dyad plot
  dc<-sum(x$dyads)          #Get # dyads
  dp<-x$dyads.pred          #Get dyad probs
  dpm<-x$dyads.pred*dc      #Get pred marginals
  dpsd<-sqrt(dp*(1-dp)*dc)  #Get pred SD
  dr<-range(c(x$dyads,dpm+1.96*dpsd,dpm-1.96*dpsd))  #Get range
  if(all(x$dyads>0)&&(all(dpm>0)))
    plot(1:3,dpm,axes=FALSE,ylim=dr,main="Predicted Dyad Census", xlab="Dyad Type",ylab="Count",log="y",col=2,xlim=c(0.5,3.5))
  else
    plot(1:3,dpm,axes=FALSE,ylim=dr,main="Predicted Dyad Census", xlab="Dyad Type",ylab="Count",col=2,xlim=c(0.5,3.5))
  segments(1:3,dpm-1.96*dpsd,1:3,dpm+1.96*dpsd,col=2)
  segments(1:3-0.3,dpm-1.96*dpsd,1:3+0.3,dpm-1.96*dpsd,col=2)
  segments(1:3-0.3,dpm+1.96*dpsd,1:3+0.3,dpm+1.96*dpsd,col=2)
  points(1:3,x$dyads,pch=19)
  axis(2)
  axis(1,at=1:3,labels=names(x$dyads),las=3)
  #Triad plot
  tc<-sum(x$triads)          #Get # triads
  tp<-x$triads.pred          #Get triad probs
  tpm<-x$triads.pred*tc      #Get pred marginals
  tpsd<-sqrt(tp*(1-tp)*tc)   #Get pred SD
  tr<-range(c(x$triads,tpm+1.96*tpsd,tpm-1.96*tpsd))  #Get range
  if(all(x$triads>0)&&(all(tpm>0)))
    plot(1:16,tpm,axes=FALSE,ylim=tr,main="Predicted Triad Census", xlab="Triad Type",ylab="Count",log="y",col=2)
  else
    plot(1:16,tpm,axes=FALSE,ylim=tr,main="Predicted Triad Census", xlab="Triad Type",ylab="Count",col=2)
  segments(1:16,tpm-1.96*tpsd,1:16,tpm+1.96*tpsd,col=2)
  segments(1:16-0.3,tpm-1.96*tpsd,1:16+0.3,tpm-1.96*tpsd,col=2)
  segments(1:16-0.3,tpm+1.96*tpsd,1:16+0.3,tpm+1.96*tpsd,col=2)
  points(1:16,x$triads,pch=19)
  axis(2)
  axis(1,at=1:16,labels=names(x$triads),las=3)
  #Structure statistics
  ssr<-range(c(x$ss,x$ss.pred))
  plot(0:(length(x$ss)-1),x$ss,type="b",xlab="Distance",ylab="Proportion Reached", main="Predicted Structure Statistics",ylim=ssr)
  lines(0:(length(x$ss)-1),x$ss.pred,col=2,lty=2)
}


#plot.lnam - Plot method for lnam
plot.lnam<-function(x,...){
   r<-residuals(x)
   f<-fitted(x)
   d<-x$disturbances
   sdr<-sd(r)
   ci<-c(-1.959964,1.959964)
   old.par <- par(no.readonly = TRUE)
   on.exit(par(old.par))
   par(mfrow=c(2,2))
   #Plot residual versus actual values
   plot(x$y,f,ylab=expression(hat(y)),xlab=expression(y),main="Fitted vs. Observed Values")
   abline(ci[1]*sdr,1,lty=3)
   abline(0,1,lty=2)
   abline(ci[2]*sdr,1,lty=3)
   #Plot disturbances versus fitted values
   plot(f,d,ylab=expression(hat(nu)),xlab=expression(hat(y)), ylim=c(min(ci[1]*x$sigma,d),max(ci[2]*x$sigma,d)),main="Fitted Values vs. Estimated Disturbances")
   abline(h=c(ci[1]*x$sigma,0,ci[2]*x$sigma),lty=c(3,2,3))
   #QQ-Plot the residuals
   qqnorm(r,main="Normal Q-Q Residual Plot")
   qqline(r)
   #Plot an influence diagram
   if(!(is.null(x$W1)&&is.null(x$W2))){
      inf<-matrix(0,ncol=x$df.total,nrow=x$df.total)
      if(!is.null(x$W1))
         inf<-inf+qr.solve(diag(x$df.total)-apply(sweep(x$W1,1,x$rho1,"*"), c(2,3),sum))
      if(!is.null(x$W2))
         inf<-inf+qr.solve(diag(x$df.total)-apply(sweep(x$W2,1,x$rho2,"*"), c(2,3),sum))
      syminf<-abs(inf)+abs(t(inf))
      diag(syminf)<-0
      infco<-cmdscale(as.dist(max(syminf)-syminf),k=2)
      diag(inf)<-NA
      stdinf<-inf-mean(inf,na.rm=TRUE)
      infsd<-sd(as.vector(stdinf),na.rm=TRUE)
      stdinf<-stdinf/infsd
      gplot(abs(stdinf),thresh=1.96,coord=infco,main="Net Influence Plot",edge.lty=1,edge.lwd=abs(stdinf)/2,edge.col=2+(inf>0)) 
   }
   #Restore plot settings
   invisible()
}


#potscalered.mcmc - Potential scale reduction (sqrt(Rhat)) for scalar estimands.
#Input must be a matrix whose columns correspond to replicate chains.  This, 
#clearly, doesn't belong here, but lacking a better place to put it I have 
#included it nonetheless.
potscalered.mcmc<-function(psi){
   #Use Gelman et al. notation, for convenience
   J<-dim(psi)[2]
   n<-dim(psi)[1]
   #Find between-group variance estimate
   mpsij<-apply(psi,2,mean)
   mpsitot<-mean(mpsij)
   B<-(n/(J-1))*sum((mpsij-mpsitot)^2)
   #Find within-group variance estimate
   s2j<-apply(psi,2,var)
   W<-mean(s2j)
   #Calculate the (estimated) marginal posterior variance of the estimand
   varppsi<-((n-1)/n)*W+(1/n)*B
   #Return the potential scale reduction estimate
   sqrt(varppsi/W)
}


#print.bayes.factor - A fairly generic routine for printing bayes factors, here used for the bbnam routine.
print.bayes.factor<-function(x,...){
   tab<-x$int.lik
   rownames(tab)<-x$model.names
   colnames(tab)<-x$model.names
   cat("Log Bayes Factors by Model:\n\n(Diagonals indicate raw integrated log likelihood estimates.)\n\n")
   print(tab)
   cat("\n")   
}


#print.bbnam - Print method for bbnam
print.bbnam<-function(x,...){
   UseMethod("print",x)
}


#print.bbnam.actor - Print method for bbnam.actor
print.bbnam.actor<-function(x,...){
   cat("\nButts' Hierarchical Bayes Model for Network Estimation/Informant Accuracy\n\n")
   cat("Multiple Error Probability Model\n\n")
   #Dump marginal posterior network
   cat("Marginal Posterior Network Distribution:\n\n")
   d<-apply(x$net,c(2,3),mean)
   rownames(d)<-as.vector(x$anames)
   colnames(d)<-as.vector(x$anames)
   print.table(d,digits=2)
   cat("\n")
   #Dump summary of error probabilities
   cat("Marginal Posterior Global Error Distribution:\n\n")
   d<-matrix(ncol=2,nrow=6)
   d[1:3,1]<-quantile(x$em,c(0,0.25,0.5),names=FALSE,na.rm=TRUE)
   d[4,1]<-mean(x$em,na.rm=TRUE)
   d[5:6,1]<-quantile(x$em,c(0.75,1.0),names=FALSE,na.rm=TRUE)
   d[1:3,2]<-quantile(x$ep,c(0,0.25,0.5),names=FALSE,na.rm=TRUE)
   d[4,2]<-mean(x$ep,na.rm=TRUE)
   d[5:6,2]<-quantile(x$ep,c(0.75,1.0),names=FALSE,na.rm=TRUE)
   colnames(d)<-c("e^-","e^+")
   rownames(d)<-c("Min","1stQ","Median","Mean","3rdQ","Max")
   print.table(d,digits=4)
   cat("\n")
}


#print.bbnam.fixed - Print method for bbnam.fixed
print.bbnam.fixed<-function(x,...){
   cat("\nButts' Hierarchical Bayes Model for Network Estimation/Informant Accuracy\n\n")
   cat("Fixed Error Probability Model\n\n")
   #Dump marginal posterior network
   cat("Marginal Posterior Network Distribution:\n\n")
   d<-apply(x$net,c(2,3),mean)
   rownames(d)<-as.vector(x$anames)
   colnames(d)<-as.vector(x$anames)
   print.table(d,digits=2)
   cat("\n")
}


#print.bbnam.pooled - Print method for bbnam.pooled
print.bbnam.pooled<-function(x,...){
   cat("\nButts' Hierarchical Bayes Model for Network Estimation/Informant Accuracy\n\n")
   cat("Pooled Error Probability Model\n\n")
   #Dump marginal posterior network
   cat("Marginal Posterior Network Distribution:\n\n")
   d<-apply(x$net,c(2,3),mean)
   rownames(d)<-as.vector(x$anames)
   colnames(d)<-as.vector(x$anames)
   print.table(d,digits=2)
   cat("\n")
   #Dump summary of error probabilities
   cat("Marginal Posterior Global Error Distribution:\n\n")
   d<-matrix(ncol=2,nrow=6)
   d[1:3,1]<-quantile(x$em,c(0,0.25,0.5),names=FALSE,na.rm=TRUE)
   d[4,1]<-mean(x$em,na.rm=TRUE)
   d[5:6,1]<-quantile(x$em,c(0.75,1.0),names=FALSE,na.rm=TRUE)
   d[1:3,2]<-quantile(x$ep,c(0,0.25,0.5),names=FALSE,na.rm=TRUE)
   d[4,2]<-mean(x$ep,na.rm=TRUE)
   d[5:6,2]<-quantile(x$ep,c(0.75,1.0),names=FALSE,na.rm=TRUE)
   colnames(d)<-c("e^-","e^+")
   rownames(d)<-c("Min","1stQ","Median","Mean","3rdQ","Max")
   print.table(d,digits=4)
   cat("\n")
}


#print.bn - Print method for summary.bn
print.bn<-function(x, digits=max(4,getOption("digits")-3), ...){
  cat("\nBiased Net Model\n\n")
  cat("Parameters:\n\n")
  cmat<-matrix(c(x$d,x$pi,x$sigma,x$rho),ncol=1)
  colnames(cmat)<-"Estimate"
  rownames(cmat)<-c("d","pi","sigma","rho")
  printCoefmat(cmat,digits=digits,...)
  cat("\n")
}


#print.lnam - Print method for lnam
print.lnam<-function(x, digits = max(3, getOption("digits") - 3), ...){
   cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
   cat("Coefficients:\n")
   print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
   cat("\n")
}


#print.netcancor - Print method for netcancor
print.netcancor<-function(x,...){
   cat("\nCanonical Network Correlation\n\n")

   cat("Canonical Correlations:\n\n")
   cmat<-matrix(data=x$cor,ncol=length(x$cor),nrow=1)
   rownames(cmat)<-""
   colnames(cmat)<-as.vector(x$cnames)
   print.table(cmat)
   cat("\n")
   cat("Pr(>=cor):\n\n")
   cmat <- matrix(data=format(x$cpgreq),ncol=length(x$cpgreq),nrow=1)
   colnames(cmat) <- as.vector(x$cnames)
   rownames(cmat)<- ""
   print.table(cmat)
   cat("\n")
   cat("Pr(<=cor):\n\n")
   cmat <- matrix(data=format(x$cpleeq),ncol=length(x$cpleeq),nrow=1)
   colnames(cmat) <- as.vector(x$cnames)
   rownames(cmat)<- ""
   print.table(cmat)
   cat("\n")

   cat("X Coefficients:\n\n")
   cmat <- format(x$xcoef)
   colnames(cmat) <- as.vector(x$xnames)
   rownames(cmat)<- as.vector(x$xnames)
   print.table(cmat)
   cat("\n")
   cat("Pr(>=xcoef):\n\n")
   cmat <- format(x$xpgreq)
   colnames(cmat) <- as.vector(x$xnames)
   rownames(cmat)<- as.vector(x$xnames)
   print.table(cmat)
   cat("\n")
   cat("Pr(<=xcoef):\n\n")
   cmat <- format(x$xpleeq)
   colnames(cmat) <- as.vector(x$xnames)
   rownames(cmat)<- as.vector(x$xnames)
   print.table(cmat)
   cat("\n")

   cat("Y Coefficients:\n\n")
   cmat <- format(x$ycoef)
   colnames(cmat) <- as.vector(x$ynames)
   rownames(cmat)<- as.vector(x$ynames)
   print.table(cmat)
   cat("\n")
   cat("Pr(>=ycoef):\n\n")
   cmat <- format(x$ypgreq)
   colnames(cmat) <- as.vector(x$ynames)
   rownames(cmat)<- as.vector(x$ynames)
   print.table(cmat)
   cat("\n")
   cat("Pr(<=ycoef):\n\n")
   cmat <- format(x$ypleeq)
   colnames(cmat) <- as.vector(x$ynames)
   rownames(cmat)<- as.vector(x$ynames)
   print.table(cmat)
   cat("\n")
}


#print.netlm - Print method for netlm
print.netlm<-function(x,...){
   cat("\nOLS Network Model\n\n")
   cat("Coefficients:\n")
   cmat <- as.vector(format(as.numeric(x$coefficients)))
   cmat <- cbind(cmat, as.vector(format(x$pleeq)))
   cmat <- cbind(cmat, as.vector(format(x$pgreq)))
   cmat <- cbind(cmat, as.vector(format(x$pgreqabs)))
   colnames(cmat) <- c("Estimate", "Pr(<=b)", "Pr(>=b)","Pr(>=|b|)")
   rownames(cmat)<- as.vector(x$names)
   print.table(cmat)
   #Goodness of fit measures
   mss<-if(x$intercept)
      sum((fitted(x)-mean(fitted(x)))^2)
   else
      sum(fitted(x)^2)
   rss<-sum(resid(x)^2)
   qn<-NROW(x$qr$qr)
   df.int<-x$intercept
   rdf<-qn-x$rank
   resvar<-rss/rdf
   fstatistic<-c(value=(mss/(x$rank-df.int))/resvar,numdf=x$rank-df.int, dendf=rdf)
   r.squared<-mss/(mss+rss)
   adj.r.squared<-1-(1-r.squared)*((qn-df.int)/rdf)
   sigma<-sqrt(resvar)
   cat("\nResidual standard error:",format(sigma,digits=4),"on",rdf,"degrees of freedom\n")
   cat("F-statistic:",formatC(fstatistic[1],digits=4),"on",fstatistic[2],"and", fstatistic[3],"degrees of freedom, p-value:",formatC(1-pf(fstatistic[1],fstatistic[2],fstatistic[3]),digits=4),"\n")
   cat("Multiple R-squared:",format(r.squared,digits=4),"\t")
   cat("Adjusted R-squared:",format(adj.r.squared,digits=4),"\n")
   cat("\n")
}


#print.netlogit - Print method for netlogit
print.netlogit<-function(x,...){
   cat("\nNetwork Logit Model\n\n")
   cat("Coefficients:\n")
   cmat <- as.vector(format(as.numeric(x$coefficients)))
   cmat <- cbind(cmat, as.vector(format(exp(as.numeric(x$coefficients)))))
   cmat <- cbind(cmat, as.vector(format(x$pleeq)))
   cmat <- cbind(cmat, as.vector(format(x$pgreq)))
   cmat <- cbind(cmat, as.vector(format(x$pgreqabs)))
   colnames(cmat) <- c("Estimate", "Exp(b)", "Pr(<=b)", "Pr(>=b)", "Pr(>=|b|)")
   rownames(cmat)<- as.vector(x$names)
   print.table(cmat)
   cat("\nGoodness of Fit Statistics:\n")
   cat("\nNull deviance:",x$null.deviance,"on",x$df.null,"degrees of freedom\n")
   cat("Residual deviance:",x$deviance,"on",x$df.residual,"degrees of freedom\n")
   cat("Chi-Squared test of fit improvement:\n\t",x$null.deviance-x$deviance,"on",x$df.null-x$df.residual,"degrees of freedom, p-value",1-pchisq(x$null.deviance-x$deviance,df=x$df.null-x$df.residual),"\n") 
   cat("AIC:",x$aic,"\tBIC:",x$bic,"\nPseudo-R^2 Measures:\n\t(Dn-Dr)/(Dn-Dr+dfn):",(x$null.deviance-x$deviance)/(x$null.deviance-x$deviance+x$df.null),"\n\t(Dn-Dr)/Dn:",1-x$deviance/x$null.deviance,"\n")
   cat("\n")
}


#print.summary.bayes.factor - Printing for bayes factor summary objects
print.summary.bayes.factor<-function(x,...){
   cat("Log Bayes Factors by Model:\n\n(Diagonals indicate raw integrated log likelihood estimates.)\n\n")
   print(x$int.lik)
   stdtab<-matrix(x$int.lik.std,nrow=1)
   colnames(stdtab)<-x$model.names
   cat("\n\nLog Inverse Bayes Factors:\n\n(Diagonals indicate log posterior probability of model under within-set choice constraints and uniform model priors.\n\n")
   print(x$inv.bf)
   cat("\nEstimated model probabilities (within-set):\n")
   temp<-exp(diag(x$inv.bf))
   names(temp)<-x$model.names
   print(temp)
   cat("\n\nDiagnostics:\n\nReplications - ",x$reps,"\n\nLog std deviations of integrated likelihood estimates:\n")
   names(x$int.lik.std)<-x$model.names
   print(x$int.lik.std)
   cat("\n\nVector of hyperprior parameters:\n\n")
   priortab<-matrix(x$prior.param,nrow=1,ncol=length(x$prior.param))
   colnames(priortab)<-x$prior.param.names
   print(priortab)
   cat("\n\n")   
}


#print.summary.bbnam - Print method for summary.bbnam
print.summary.bbnam<-function(x,...){
   UseMethod("print",x)
}


#print.summary.bbnam.actor - Print method for summary.bbnam.actor
print.summary.bbnam.actor<-function(x,...){
   cat("\nButts' Hierarchical Bayes Model for Network Estimation/Informant Accuracy\n\n")
   cat("Multiple Error Probability Model\n\n")
   #Dump marginal posterior network
   cat("Marginal Posterior Network Distribution:\n\n")
   d<-apply(x$net,c(2,3),mean)
   rownames(d)<-as.vector(x$anames)
   colnames(d)<-as.vector(x$anames)
   print.table(d,digits=2)
   cat("\n")
   #Dump summary of error probabilities
   cat("Marginal Posterior Global Error Distribution:\n\n")
   d<-matrix(ncol=2,nrow=6)
   d[1:3,1]<-quantile(x$em,c(0,0.25,0.5),names=FALSE,na.rm=TRUE)
   d[4,1]<-mean(x$em,na.rm=TRUE)
   d[5:6,1]<-quantile(x$em,c(0.75,1.0),names=FALSE,na.rm=TRUE)
   d[1:3,2]<-quantile(x$ep,c(0,0.25,0.5),names=FALSE,na.rm=TRUE)
   d[4,2]<-mean(x$ep,na.rm=TRUE)
   d[5:6,2]<-quantile(x$ep,c(0.75,1.0),names=FALSE,na.rm=TRUE)
   colnames(d)<-c("e^-","e^+")
   rownames(d)<-c("Min","1stQ","Median","Mean","3rdQ","Max")
   print.table(d,digits=4)
   cat("\n")
   #Dump error probability estimates per observer
   cat("Marginal Posterior Error Distribution (by observer):\n\n")
   cat("Probability of False Negatives (e^-):\n\n")
   d<-matrix(ncol=6)
   for(i in 1:x$nobservers){
      dv<-matrix(c(quantile(x$em[,i],c(0,0.25,0.5),names=FALSE,na.rm=TRUE),mean(x$em[,i],na.rm=TRUE),quantile(x$em[,i],c(0.75,1.0),names=FALSE,na.rm=TRUE)),nrow=1,ncol=6)
      d<-rbind(d,dv)
   }
   d<-d[2:(x$nobservers+1),]
   rownames(d)<-as.vector(x$onames)
   colnames(d)<-c("Min","1stQ","Median","Mean","3rdQ","Max")
   print.table(d,digits=4)
   cat("\n")
   cat("Probability of False Positives (e^+):\n\n")
   d<-matrix(ncol=6)
   for(i in 1:x$nobservers){
      dv<-matrix(c(quantile(x$ep[,i],c(0,0.25,0.5),names=FALSE,na.rm=TRUE),mean(x$ep[,i],na.rm=TRUE),quantile(x$ep[,i],c(0.75,1.0),names=FALSE,na.rm=TRUE)),nrow=1,ncol=6)
      d<-rbind(d,dv)
   }
   d<-d[2:(x$nobservers+1),]
   rownames(d)<-as.vector(x$onames)
   colnames(d)<-c("Min","1stQ","Median","Mean","3rdQ","Max")
   print.table(d,digits=4)
   cat("\n")
   #Dump MCMC diagnostics
   cat("MCMC Diagnostics:\n\n")
   cat("\tReplicate Chains:",x$reps,"\n")
   cat("\tBurn Time:",x$burntime,"\n")
   cat("\tDraws per Chain:",x$draws/x$reps,"Total Draws:",x$draws,"\n")
   if("sqrtrhat" %in% names(x))
      cat("\tPotential Scale Reduction (G&R's sqrt(Rhat)):\n \t\tMax:",max(x$sqrtrhat[!is.nan(x$sqrtrhat)]),"\n\t\tMed:",median(x$sqrtrhat[!is.nan(x$sqrtrhat)]),"\n\t\tIQR:",IQR(x$sqrtrhat[!is.nan(x$sqrtrhat)]),"\n")
   cat("\n")
}


#print.summary.bbnam.fixed - Print method for summary.bbnam.fixed
print.summary.bbnam.fixed<-function(x,...){
   cat("\nButts' Hierarchical Bayes Model for Network Estimation/Informant Accuracy\n\n")
   cat("Fixed Error Probability Model\n\n")
   #Dump marginal posterior network
   cat("Marginal Posterior Network Distribution:\n\n")
   d<-apply(x$net,c(2,3),mean)
   rownames(d)<-as.vector(x$anames)
   colnames(d)<-as.vector(x$anames)
   print.table(d,digits=2)
   cat("\n")
   #Dump model diagnostics
   cat("Model Diagnostics:\n\n")
   cat("\tTotal Draws:",x$draws,"\n\t(Note: Draws taken directly from network posterior.)")
   cat("\n")
}


#print.summary.bbnam.pooled - Print method for summary.bbnam.pooled
print.summary.bbnam.pooled<-function(x,...){
   cat("\nButts' Hierarchical Bayes Model for Network Estimation/Informant Accuracy\n\n")
   cat("Pooled Error Probability Model\n\n")
   #Dump marginal posterior network
   cat("Marginal Posterior Network Distribution:\n\n")
   d<-apply(x$net,c(2,3),mean)
   rownames(d)<-as.vector(x$anames)
   colnames(d)<-as.vector(x$anames)
   print.table(d,digits=2)
   cat("\n")
   #Dump summary of error probabilities
   cat("Marginal Posterior Error Distribution:\n\n")
   d<-matrix(ncol=2,nrow=6)
   d[1:3,1]<-quantile(x$em,c(0,0.25,0.5),names=FALSE,na.rm=TRUE)
   d[4,1]<-mean(x$em,na.rm=TRUE)
   d[5:6,1]<-quantile(x$em,c(0.75,1.0),names=FALSE,na.rm=TRUE)
   d[1:3,2]<-quantile(x$ep,c(0,0.25,0.5),names=FALSE,na.rm=TRUE)
   d[4,2]<-mean(x$ep,na.rm=TRUE)
   d[5:6,2]<-quantile(x$ep,c(0.75,1.0),names=FALSE,na.rm=TRUE)
   colnames(d)<-c("e^-","e^+")
   rownames(d)<-c("Min","1stQ","Median","Mean","3rdQ","Max")
   print.table(d,digits=4)
   cat("\n")
   #Dump MCMC diagnostics
   cat("MCMC Diagnostics:\n\n")
   cat("\tReplicate Chains:",x$reps,"\n")
   cat("\tBurn Time:",x$burntime,"\n")
   cat("\tDraws per Chain:",x$draws/x$reps,"Total Draws:",x$draws,"\n")
   if("sqrtrhat" %in% names(x))
      cat("\tPotential Scale Reduction (G&R's sqrt(Rhat)):\n \t\tMax:",max(x$sqrtrhat[!is.nan(x$sqrtrhat)]),"\n\t\tMed:",median(x$sqrtrhat[!is.nan(x$sqrtrhat)]),"\n\t\tIQR:",IQR(x$sqrtrhat[!is.nan(x$sqrtrhat)]),"\n")
   cat("\n")
}


#print.summary.bn - Print method for summary.bn
print.summary.bn<-function(x, digits=max(4,getOption("digits")-3), signif.stars=getOption("show.signif.stars"), ...){
  cat("\nBiased Net Model\n\n")
  cat("\nParameters:\n\n")
  cmat<-matrix(c(x$d,x$pi,x$sigma,x$rho),ncol=1)
  colnames(cmat)<-"Estimate"
  rownames(cmat)<-c("d","pi","sigma","rho")
  printCoefmat(cmat,digits=digits,...)
  #Diagnostics
  cat("\nDiagnostics:\n\n")
  cat("\tFit method:",x$method,"\n")
  cat("\tPseudolikelihood G^2:",x$G.square,"\n")
  #Plot edge census
  cat("\n\tEdge census comparison:\n\n")
  ec<-sum(x$edges)
  cmat<-cbind(x$edges,x$edges.pred*ec)
  cmat<-cbind(cmat,(cmat[,1]-cmat[,2])/sqrt(x$edges.pred*(1-x$edges.pred)*ec))
  cmat<-cbind(cmat,2*(1-pnorm(abs(cmat[,3]))))
  colnames(cmat)<-c("Observed","Predicted","Z Value","Pr(>|z|)")
  printCoefmat(cmat,digits=digits,signif.stars=signif.stars,...)
  chisq<-sum((cmat[,1]-cmat[,2])^2/cmat[,2])
  cat("\tChi-Square:",chisq,"on 1 degrees of freedom.  p-value:",1-pchisq(chisq,1),"\n\n")
  #Plot dyad census
  cat("\n\tDyad census comparison:\n\n")
  dc<-sum(x$dyads)
  cmat<-cbind(x$dyads,x$dyads.pred*dc)
  cmat<-cbind(cmat,(cmat[,1]-cmat[,2])/sqrt(x$dyads.pred*(1-x$dyads.pred)*dc))
  cmat<-cbind(cmat,2*(1-pnorm(abs(cmat[,3]))))
  colnames(cmat)<-c("Observed","Predicted","Z Value","Pr(>|z|)")
  printCoefmat(cmat,digits=digits,signif.stars=signif.stars,...)
  chisq<-sum((cmat[,1]-cmat[,2])^2/cmat[,2])
  cat("\tChi-Square:",chisq,"on 2 degrees of freedom.  p-value:",1-pchisq(chisq,2),"\n\n")
  #Plot triad census
  cat("\n\tTriad census comparison:\n\n")
  tc<-sum(x$triads)
  cmat<-cbind(x$triads,x$triads.pred*tc)
  cmat<-cbind(cmat,(cmat[,1]-cmat[,2])/sqrt(x$triads.pred*(1-x$triads.pred)*tc))
  cmat<-cbind(cmat,2*(1-pnorm(abs(cmat[,3]))))
  colnames(cmat)<-c("Observed","Predicted","Z Value","Pr(>|z|)")
  printCoefmat(cmat,digits=digits,signif.stars=signif.stars,...)
  chisq<-sum((cmat[,1]-cmat[,2])^2/cmat[,2])
  cat("\tChi-Square:",chisq,"on 15 degrees of freedom.  p-value:",1-pchisq(chisq,15),"\n\n")
}


#print.summary.brokerage - print method for summary.brokerage objects
print.summary.brokerage<-function(x,...){
  cat("Gould-Fernandez Brokerage Analysis\n\n")
  cat("Global Brokerage Properties\n")
  cmat<-cbind(x$raw.gli,x$exp.gli,x$sd.gli,x$z.gli,2*(1-pnorm(abs(x$z.gli))))
  rownames(cmat)<-names(x$raw.gli)
  colnames(cmat)<-c("t","E(t)","Sd(t)","z","Pr(>|z|)")
  printCoefmat(cmat)
  cat("\nIndividual Properties (by Group)\n")
  for(i in x$clid){
    cat("\n\tGroup ID:",i,"\n")
    temp<-x$cl==i
    cmat<-cbind(x$raw.nli,x$z.nli)[temp,,drop=FALSE]
    print(cmat)
  }
  cat("\n")
}


#print.summary.lnam - Print method for summary.lnam
print.summary.lnam<-function(x, digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"), ...){
   cat("\nCall:\n")
   cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
   cat("Residuals:\n")
   nam <- c("Min", "1Q", "Median", "3Q", "Max")
   resid<-x$residuals 
   rq <- if (length(dim(resid)) == 2) 
      structure(apply(t(resid), 1, quantile), dimnames = list(nam, dimnames(resid)[[2]]))
   else structure(quantile(resid), names = nam)
   print(rq, digits = digits, ...)
   cat("\nCoefficients:\n")
   cmat<-cbind(coef(x),se.lnam(x))
   cmat<-cbind(cmat,cmat[,1]/cmat[,2],(1-pnorm(abs(cmat[,1]),0,cmat[,2]))*2)
   colnames(cmat)<-c("Estimate","Std. Error","Z value","Pr(>|z|)")
   #print(format(cmat,digits=digits),quote=FALSE)
   printCoefmat(cmat,digits=digits,signif.stars=signif.stars,...)
   cat("\n")
   cmat<-cbind(x$sigma,x$sigma.se)
   colnames(cmat)<-c("Estimate","Std. Error")
   rownames(cmat)<-"Sigma"
   printCoefmat(cmat,digits=digits,signif.stars=signif.stars,...)
   cat("\nGoodness-of-Fit:\n")
   rss<-sum(x$residuals^2)
   mss<-sum((x$fitted-mean(x$fitted))^2)
   rdfns<-x$df.residual+1
   cat("\tResidual standard error: ",format(sqrt(rss/rdfns),digits=digits)," on ",rdfns," degrees of freedom (w/o Sigma)\n",sep="")
   cat("\tMultiple R-Squared: ",format(mss/(mss+rss),digits=digits),", Adjusted R-Squared: ",format(1-(1-mss/(mss+rss))*x$df.total/rdfns,digits=digits),"\n",sep="")
   cat("\tModel log likelihood:", format(x$lnlik.model,digits=digits), "on", x$df.resid, "degrees of freedom (w/Sigma)\n\tAIC:",format(-2*x$lnlik.model+2*x$df.model,digits=digits),"BIC:",format(-2*x$lnlik.model+log(x$df.total)*x$df.model,digits=digits),"\n")
   cat("\n\tNull model:",x$null.model,"\n")
   cat("\tNull log likelihood:", format(x$lnlik.null,digits=digits), "on", x$df.null.resid, "degrees of freedom\n\tAIC:",format(-2*x$lnlik.null+2*x$df.null,digits=digits),"BIC:",format(-2*x$lnlik.null+log(x$df.total)*x$df.null,digits=digits),"\n")
   cat("\tAIC difference (model versus null):",format(-2*x$lnlik.null+2*x$df.null+2*x$lnlik.model-2*x$df.model,digits=digits),"\n")
   cat("\tHeuristic Log Bayes Factor (model versus null): ",format(-2*x$lnlik.null+log(x$df.total)*x$df.null+2*x$lnlik.model-log(x$df.total)*x$df.model,digits=digits),"\n")
   cat("\n")
}


#print.summary.netcancor - Print method for summary.netcancor
print.summary.netcancor<-function(x,...){
   cat("\nCanonical Network Correlation\n\n")

   cat("Canonical Correlations:\n\n")
   cmat<-as.vector(x$cor)
   cmat<-rbind(cmat,as.vector((x$cor)^2))
   rownames(cmat)<-c("Correlation","Coef. of Det.")
   colnames(cmat)<-as.vector(x$cnames)
   print.table(cmat)
   cat("\n")
   cat("Pr(>=cor):\n\n")
   cmat <- matrix(data=format(x$cpgreq),ncol=length(x$cpgreq),nrow=1)
   colnames(cmat) <- as.vector(x$cnames)
   rownames(cmat)<- ""
   print.table(cmat)
   cat("\n")
   cat("Pr(<=cor):\n\n")
   cmat <- matrix(data=format(x$cpleeq),ncol=length(x$cpleeq),nrow=1)
   colnames(cmat) <- as.vector(x$cnames)
   rownames(cmat)<- ""
   print.table(cmat)
   cat("\n")

   cat("X Coefficients:\n\n")
   cmat <- format(x$xcoef)
   colnames(cmat) <- as.vector(x$xnames)
   rownames(cmat)<- as.vector(x$xnames)
   print.table(cmat)
   cat("\n")
   cat("Pr(>=xcoef):\n\n")
   cmat <- format(x$xpgreq)
   colnames(cmat) <- as.vector(x$xnames)
   rownames(cmat)<- as.vector(x$xnames)
   print.table(cmat)
   cat("\n")
   cat("Pr(<=xcoef):\n\n")
   cmat <- format(x$xpleeq)
   colnames(cmat) <- as.vector(x$xnames)
   rownames(cmat)<- as.vector(x$xnames)
   print.table(cmat)
   cat("\n")

   cat("Y Coefficients:\n\n")
   cmat <- format(x$ycoef)
   colnames(cmat) <- as.vector(x$ynames)
   rownames(cmat)<- as.vector(x$ynames)
   print.table(cmat)
   cat("\n")
   cat("Pr(>=ycoef):\n\n")
   cmat <- format(x$ypgreq)
   colnames(cmat) <- as.vector(x$ynames)
   rownames(cmat)<- as.vector(x$ynames)
   print.table(cmat)
   cat("\n")
   cat("Pr(<=ycoef):\n\n")
   cmat <- format(x$ypleeq)
   colnames(cmat) <- as.vector(x$ynames)
   rownames(cmat)<- as.vector(x$ynames)
   print.table(cmat)
   cat("\n")

   cat("Test Diagnostics:\n\n")
   cat("\tNull Hypothesis:")
   if(x$nullhyp=="qap")
      cat(" QAP\n")
   else
      cat(" CUG\n")
   cat("\tReplications:",dim(x$cdist)[1],"\n")
   cat("\tDistribution Summary for Correlations:\n\n")
   dmat<-apply(x$cdist,2,min,na.rm=TRUE)
   dmat<-rbind(dmat,apply(x$cdist,2,quantile,probs=0.25,names=FALSE,na.rm=TRUE))
   dmat<-rbind(dmat,apply(x$cdist,2,quantile,probs=0.5,names=FALSE,na.rm=TRUE))
   dmat<-rbind(dmat,apply(x$cdist,2,mean,na.rm=TRUE))
   dmat<-rbind(dmat,apply(x$cdist,2,quantile,probs=0.75,names=FALSE,na.rm=TRUE))
   dmat<-rbind(dmat,apply(x$cdist,2,max,na.rm=TRUE))
   colnames(dmat)<-as.vector(x$cnames)
   rownames(dmat)<-c("Min","1stQ","Median","Mean","3rdQ","Max")
   print.table(dmat,digits=4)
   cat("\n")
}


#print.summary.netlm - Print method for summary.netlm
print.summary.netlm<-function(x,...){
   cat("\nOLS Network Model\n\n")
   #Residuals
   cat("Residuals:\n")
   print.table(format(quantile(x$residuals)))
   #Coefficients
   cat("\nCoefficients:\n")
   cmat <- as.vector(format(as.numeric(x$coefficients)))
   cmat <- cbind(cmat, as.vector(format(x$pleeq)))
   cmat <- cbind(cmat, as.vector(format(x$pgreq)))
   cmat <- cbind(cmat, as.vector(format(x$pgreqabs)))
   colnames(cmat) <- c("Estimate", "Pr(<=b)", "Pr(>=b)", "Pr(>=|b|)")
   rownames(cmat)<- as.vector(x$names)
   print.table(cmat)
   #Goodness of fit measures
   mss<-if(x$intercept)
      sum((fitted(x)-mean(fitted(x)))^2)
   else
      sum(fitted(x)^2)
   rss<-sum(resid(x)^2)
   qn<-NROW(x$qr$qr)
   df.int<-x$intercept
   rdf<-qn-x$rank
   resvar<-rss/rdf
   fstatistic<-c(value=(mss/(x$rank-df.int))/resvar,numdf=x$rank-df.int, dendf=rdf)
   r.squared<-mss/(mss+rss)
   adj.r.squared<-1-(1-r.squared)*((qn-df.int)/rdf)
   sigma<-sqrt(resvar)
   cat("\nResidual standard error:",format(sigma,digits=4),"on",rdf,"degrees of freedom\n")
   cat("Multiple R-squared:",format(r.squared,digits=4),"\t")
   cat("Adjusted R-squared:",format(adj.r.squared,digits=4),"\n")
   cat("F-statistic:",formatC(fstatistic[1],digits=4),"on",fstatistic[2],"and", fstatistic[3],"degrees of freedom, p-value:",formatC(1-pf(fstatistic[1],fstatistic[2],fstatistic[3]),digits=4),"\n")
   #Test diagnostics
   cat("\n\nTest Diagnostics:\n\n")
   cat("\tNull Hypothesis:",x$nullhyp,"\n")
   if(!is.null(x$dist)){
     cat("\tReplications:",dim(x$dist)[1],"\n")
     cat("\tCoefficient Distribution Summary:\n\n")
     dmat<-apply(x$dist,2,min,na.rm=TRUE)
     dmat<-rbind(dmat,apply(x$dist,2,quantile,probs=0.25,names=FALSE, na.rm=TRUE))
     dmat<-rbind(dmat,apply(x$dist,2,quantile,probs=0.5,names=FALSE,na.rm=TRUE))
     dmat<-rbind(dmat,apply(x$dist,2,mean,na.rm=TRUE))
     dmat<-rbind(dmat,apply(x$dist,2,quantile,probs=0.75,names=FALSE, na.rm=TRUE))
     dmat<-rbind(dmat,apply(x$dist,2,max,na.rm=TRUE))
     colnames(dmat)<-as.vector(x$names)
     rownames(dmat)<-c("Min","1stQ","Median","Mean","3rdQ","Max")
     print.table(dmat,digits=4)
     cat("\n")
  }
}


#print.summary.netlogit - Print method for summary.netlogit
print.summary.netlogit<-function(x,...){
   cat("\nNetwork Logit Model\n\n")
   cat("Coefficients:\n")
   cmat <- as.vector(format(as.numeric(x$coefficients)))
   cmat <- cbind(cmat, as.vector(format(exp(as.numeric(x$coefficients)))))
   cmat <- cbind(cmat, as.vector(format(x$pleeq)))
   cmat <- cbind(cmat, as.vector(format(x$pgreq)))
   cmat <- cbind(cmat, as.vector(format(x$pgreqabs)))
   colnames(cmat) <- c("Estimate", "Exp(b)", "Pr(<=b)", "Pr(>=b)", "Pr(>=|b|)")
   rownames(cmat)<- as.vector(x$names)
   print.table(cmat)
   cat("\nGoodness of Fit Statistics:\n")
   cat("\nNull deviance:",x$null.deviance,"on",x$df.null,"degrees of freedom\n")
   cat("Residual deviance:",x$deviance,"on",x$df.residual,"degrees of freedom\n")
   cat("Chi-Squared test of fit improvement:\n\t",x$null.deviance-x$deviance,"on",x$df.null-x$df.residual,"degrees of freedom, p-value",1-pchisq(x$null.deviance-x$deviance,df=x$df.null-x$df.residual),"\n") 
   cat("AIC:",x$aic,"\tBIC:",x$bic,"\nPseudo-R^2 Measures:\n\t(Dn-Dr)/(Dn-Dr+dfn):",(x$null.deviance-x$deviance)/(x$null.deviance-x$deviance+x$df.null),"\n\t(Dn-Dr)/Dn:",1-x$deviance/x$null.deviance,"\n")
   cat("Contingency Table (predicted (rows) x actual (cols)):\n\n")
   print.table(x$ctable,print.gap=3)
   cat("\n\tTotal Fraction Correct:",(x$ctable[1,1]+x$ctable[2,2])/sum(x$ctable),"\n\tFraction Predicted 1s Correct:",x$ctable[2,2]/sum(x$ctable[2,]),"\n\tFraction Predicted 0s Correct:",x$ctable[1,1]/sum(x$ctable[1,]),"\n")
   cat("\tFalse Negative Rate:",x$ctable[1,2]/sum(x$ctable[,2]),"\n")
   cat("\tFalse Positive Rate:",x$ctable[2,1]/sum(x$ctable[,1]),"\n")
   cat("\nTest Diagnostics:\n\n")
   cat("\tNull Hypothesis:",x$nullhyp,"\n")
   if(!is.null(x$dist)){
     cat("\tReplications:",dim(x$dist)[1],"\n")
     cat("\tDistribution Summary:\n\n")
     dmat<-apply(x$dist,2,min,na.rm=TRUE)
     dmat<-rbind(dmat,apply(x$dist,2,quantile,probs=0.25,names=FALSE, na.rm=TRUE))
     dmat<-rbind(dmat,apply(x$dist,2,quantile,probs=0.5,names=FALSE,na.rm=TRUE))
     dmat<-rbind(dmat,apply(x$dist,2,mean,na.rm=TRUE))
     dmat<-rbind(dmat,apply(x$dist,2,quantile,probs=0.75,names=FALSE, na.rm=TRUE))
     dmat<-rbind(dmat,apply(x$dist,2,max,na.rm=TRUE))
     colnames(dmat)<-as.vector(x$names)
     rownames(dmat)<-c("Min","1stQ","Median","Mean","3rdQ","Max")
     print.table(dmat,digits=4)
     cat("\n")
   }
}


#pstar - Perform an approximate p* analysis using the logistic regression 
#approximation.  Note that the result of this is returned as a GLM object, and 
#subsequent printing/summarizing/etc. should be treated accordingly.
pstar<-function(dat,effects=c("choice","mutuality","density","reciprocity","stransitivity","wtransitivity","stranstri","wtranstri","outdegree","indegree","betweenness","closeness","degcentralization","betcentralization","clocentralization","connectedness","hierarchy","lubness","efficiency"),attr=NULL,memb=NULL,diag=FALSE,mode="digraph"){
   #First, take care of various details
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat)||(is.array(dat)&&(length(dim(dat))>2)))
     stop("Single graphs required in pstar.")
   n<-dim(dat)[1]
   m<-dim(dat)[2]
   o<-list()
   #Next, add NAs as needed
   d<-dat
   if(!diag)
      d<-diag.remove(d)
   if(mode=="graph")
      d<-upper.tri.remove(d)
   #Make sure that attr and memb are well-behaved
   if(!is.null(attr)){
      if(is.vector(attr))
         attr<-matrix(attr,ncol=1)
      if(is.null(colnames(attr)))
         colnames(attr)<-paste("Attribute",1:dim(attr)[2],sep=".")
   }
   if(!is.null(memb)){
      if(is.vector(memb))
         memb<-matrix(memb,ncol=1)
      if(is.null(colnames(memb)))
         colnames(memb)<-paste("Membership",1:dim(memb)[2],sep=".")
   }
   #Now, evaluate each specified effect given each possible perturbation
   tiedat<-vector()
   for(i in 1:n)
      for(j in 1:m)
         if(!is.na(d[i,j])){
            #Assess the effects
            td<-vector()
            if(!is.na(pmatch("choice",effects))){  #Compute a choice effect
               td<-c(td,1)  #Always constant
            }
            if(!is.na(pmatch("mutuality",effects))){  #Compute a mutuality effect
               td<-c(td,eval.edgeperturbation(d,i,j,"mutuality"))
            }
            if(!is.na(pmatch("density",effects))){  #Compute a density effect
               td<-c(td,eval.edgeperturbation(d,i,j,"gden",mode=mode,diag=diag))
            }
            if(!is.na(pmatch("reciprocity",effects))){  #Compute a reciprocity effect
               td<-c(td,eval.edgeperturbation(d,i,j,"grecip"))
            }
            if(!is.na(pmatch("stransitivity",effects))){  #Compute a strong transitivity effect
               td<-c(td,eval.edgeperturbation(d,i,j,"gtrans",mode=mode,diag=diag,measure="strong"))
            }
            if(!is.na(pmatch("wtransitivity",effects))){  #Compute a weak transitivity effect
               td<-c(td,eval.edgeperturbation(d,i,j,"gtrans",mode=mode,diag=diag,measure="weak"))
            }
            if(!is.na(pmatch("stranstri",effects))){  #Compute a strong trans census effect
               td<-c(td,eval.edgeperturbation(d,i,j,"gtrans",mode=mode,diag=diag,measure="strongcensus"))
            }
            if(!is.na(pmatch("wtranstri",effects))){  #Compute a weak trans census effect
               td<-c(td,eval.edgeperturbation(d,i,j,"gtrans",mode=mode,diag=diag,measure="weakcensus"))
            }
            if(!is.na(pmatch("outdegree",effects))){  #Compute outdegree effects
               td<-c(td,eval.edgeperturbation(d,i,j,"degree",cmode="outdegree",gmode=mode,diag=diag))
            }
            if(!is.na(pmatch("indegree",effects))){  #Compute indegree effects
               td<-c(td,eval.edgeperturbation(d,i,j,"degree",cmode="indegree",gmode=mode,diag=diag))
            }
            if(!is.na(pmatch("betweenness",effects))){  #Compute betweenness effects
               td<-c(td,eval.edgeperturbation(d,i,j,"betweenness",gmode=mode,diag=diag))
            }
            if(!is.na(pmatch("closeness",effects))){  #Compute closeness effects
               td<-c(td,eval.edgeperturbation(d,i,j,"closeness",gmode=mode,diag=diag))
            }
            if(!is.na(pmatch("degcentralization",effects))){  #Compute degree centralization effects
               td<-c(td,eval.edgeperturbation(d,i,j,"centralization","degree",mode=mode,diag=diag))
            }
            if(!is.na(pmatch("betcentralization",effects))){  #Compute betweenness centralization effects
               td<-c(td,eval.edgeperturbation(d,i,j,"centralization","betweenness",mode=mode,diag=diag))
            }
            if(!is.na(pmatch("clocentralization",effects))){  #Compute closeness centralization effects
               td<-c(td,eval.edgeperturbation(d,i,j,"centralization","closeness",mode=mode,diag=diag))
            }
            if(!is.na(pmatch("connectedness",effects))){  #Compute connectedness effects
               td<-c(td,eval.edgeperturbation(d,i,j,"connectedness"))
            }
            if(!is.na(pmatch("hierarchy",effects))){  #Compute hierarchy effects
               td<-c(td,eval.edgeperturbation(d,i,j,"hierarchy"))
            }
            if(!is.na(pmatch("lubness",effects))){  #Compute lubness effects
               td<-c(td,eval.edgeperturbation(d,i,j,"lubness"))
            }
            if(!is.na(pmatch("efficiency",effects))){  #Compute efficiency effects
               td<-c(td,eval.edgeperturbation(d,i,j,"efficiency",diag=diag))
            }
            #Add attribute differences, if needed
            if(!is.null(attr))
               td<-c(td,abs(attr[i,]-attr[j,]))
            #Add membership similarities, if needed
            if(!is.null(memb))
               td<-c(td,as.numeric(memb[i,]==memb[j,]))
            #Add this data to the aggregated tie data
            tiedat<-rbind(tiedat,c(d[i,j],td))
         }
   #Label the tie data matrix
   tiedat.lab<-"EdgeVal"
   if(!is.na(pmatch("choice",effects)))  #Label the choice effect
      tiedat.lab<-c(tiedat.lab,"Choice")
   if(!is.na(pmatch("mutuality",effects)))  #Label the mutuality effect
      tiedat.lab<-c(tiedat.lab,"Mutuality")
   if(!is.na(pmatch("density",effects)))  #Label the density effect
      tiedat.lab<-c(tiedat.lab,"Density")
   if(!is.na(pmatch("reciprocity",effects)))  #Label the reciprocity effect
      tiedat.lab<-c(tiedat.lab,"Reciprocity")
   if(!is.na(pmatch("stransitivity",effects)))  #Label the strans effect
      tiedat.lab<-c(tiedat.lab,"STransitivity")
   if(!is.na(pmatch("wtransitivity",effects)))  #Label the wtrans effect
      tiedat.lab<-c(tiedat.lab,"WTransitivity")
   if(!is.na(pmatch("stranstri",effects)))  #Label the stranstri effect
      tiedat.lab<-c(tiedat.lab,"STransTriads")
   if(!is.na(pmatch("wtranstri",effects)))  #Label the wtranstri effect
      tiedat.lab<-c(tiedat.lab,"WTransTriads")
   if(!is.na(pmatch("outdegree",effects)))  #Label the outdegree effect
      tiedat.lab<-c(tiedat.lab,paste("Outdegree",1:n,sep="."))
   if(!is.na(pmatch("indegree",effects)))  #Label the indegree effect
      tiedat.lab<-c(tiedat.lab,paste("Indegree",1:n,sep="."))
   if(!is.na(pmatch("betweenness",effects)))  #Label the betweenness effect
      tiedat.lab<-c(tiedat.lab,paste("Betweenness",1:n,sep="."))
   if(!is.na(pmatch("closeness",effects)))  #Label the closeness effect
      tiedat.lab<-c(tiedat.lab,paste("Closeness",1:n,sep="."))
   if(!is.na(pmatch("degcent",effects)))  #Label the degree centralization effect
      tiedat.lab<-c(tiedat.lab,"DegCentralization")
   if(!is.na(pmatch("betcent",effects)))  #Label the betweenness centralization effect
      tiedat.lab<-c(tiedat.lab,"BetCentralization")
   if(!is.na(pmatch("clocent",effects)))  #Label the closeness centralization effect
      tiedat.lab<-c(tiedat.lab,"CloCentralization")
   if(!is.na(pmatch("connectedness",effects)))  #Label the connectedness effect
      tiedat.lab<-c(tiedat.lab,"Connectedness")
   if(!is.na(pmatch("hierarchy",effects)))  #Label the hierarchy effect
      tiedat.lab<-c(tiedat.lab,"Hierarchy")
   if(!is.na(pmatch("lubness",effects)))  #Label the lubness effect
      tiedat.lab<-c(tiedat.lab,"LUBness")
   if(!is.na(pmatch("efficiency",effects)))  #Label the efficiency effect
      tiedat.lab<-c(tiedat.lab,"Efficiency")
   if(!is.null(attr))
      tiedat.lab<-c(tiedat.lab,colnames(attr))
   if(!is.null(memb))
      tiedat.lab<-c(tiedat.lab,colnames(memb))
   colnames(tiedat)<-tiedat.lab
   #Having had our fun, it's time to get serious.  Run a GLM on the resulting data.
   fmla<-as.formula(paste("EdgeVal ~ -1 + ",paste(colnames(tiedat)[2:dim(tiedat)[2]],collapse=" + ")))
   o<-glm(fmla,family="binomial",data=as.data.frame(tiedat))
   o$tiedata<-tiedat
   #Return the result
   o
}


#se.lnam - Standard error method for lnam
se.lnam<-function(object, ...){
   se<-vector()
   sen<-vector()
   if(!is.null(object$beta.se)){
      se<-c(se,object$beta.se)
      sen<-c(sen,names(object$beta.se))
   }
   if(!is.null(object$rho1.se)){
      se<-c(se,object$rho1.se)
      sen<-c(sen,"rho1")
   }
   if(!is.null(object$rho2.se)){
      se<-c(se,object$rho2.se)
      sen<-c(sen,"rho2")
   }
   names(se)<-sen
   se
}


#summary.bayes.factor - A fairly generic summary routine for bayes factors.  
#Clearly, this belongs in some other library than sna, but for the moment this 
#will have to do...
summary.bayes.factor<-function(object, ...){
   o<-object
   rownames(o$int.lik)<-o$model.names
   colnames(o$int.lik)<-o$model.names
   o$inv.bf<--o$int.lik
   for(i in 1:dim(o$int.lik)[1])
      o$inv.bf[i,i]<-o$int.lik[i,i]-logSum(diag(o$int.lik))
   class(o)<-c("summary.bayes.factor","bayes.factor")
   o
}


#summary.bbnam - Summary method for bbnam
summary.bbnam<-function(object, ...){
   out<-object
   class(out)<-c("summary.bbnam",class(out))
   out
}


#summary.bbnam.actor - Summary method for bbnam.actor
summary.bbnam.actor<-function(object, ...){
   out<-object
   class(out)<-c("summary.bbnam.actor",class(out))
   out
}


#summary.bbnam.fixed - Summary method for bbnam.fixed
summary.bbnam.fixed<-function(object, ...){
   out<-object
   class(out)<-c("summary.bbnam.fixed",class(out))
   out
}


#summary.bbnam.pooled - Summary method for bbnam.pooled
summary.bbnam.pooled<-function(object, ...){
   out<-object
   class(out)<-c("summary.bbnam.pooled",class(out))
   out
}


#summary.bn - Summary method for bn
summary.bn<-function(object, ...){
  out<-object
  class(out)<-c("summary.bn",class(out))
  out
}


#summary.brokerage - Summary method for brokerage objects
summary.brokerage<-function(object,...){
  class(object)<-"summary.brokerage"
  object
}


#summary.lnam - Summary method lnam
summary.lnam<-function(object, ...){
   ans<-object
   class(ans)<-c("summary.lnam","lnam")
   ans
}


#summary.netcancor - Summary method for netcancor
summary.netcancor<-function(object, ...){
   out<-object
   class(out)<-c("summary.netcancor",class(out))
   out
}


#summary.netlm - Summary method for netlm
summary.netlm<-function(object, ...){
   out<-object
   class(out)<-c("summary.netlm",class(out))
   out
}


#summary.netlogit - Summary method for netlogit
summary.netlogit<-function(object, ...){
   out<-object
   class(out)<-c("summary.netlogit",class(out))
   out
}
