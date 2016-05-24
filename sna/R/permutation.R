######################################################################
#
# permutation.R
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 4/23/05
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains routines related to permutations on graphs.
#
# Contents:
#   lab.optimize
#   lab.optimize.anneal
#   lab.optimize.exhaustive
#   lab.optimize.gumbel
#   lab.optimize.hillclimb
#   lab.optimize.mc
#   numperm
#   rmperm
#   rperm
#
######################################################################


#lab.optimize - Optimize a function over the accessible permutation groups of two or more graphs.  This routine is a front end for various method-specific functions, and is in turn intended to be called from structural distance/covariance routines and the like.  The methods supported at this time include "exhaustive" (exhaustive search - I hope these are _small_ graphs!), "mc" (simple monte carlo search), " 
lab.optimize<-function(d1,d2,FUN,exchange.list=0,seek="min",opt.method=c("anneal","exhaustive","mc","hillclimb","gumbel"),...){
   meth<-match.arg(opt.method)
   if(meth=="anneal")
      lab.optimize.anneal(d1,d2,FUN,exchange.list,seek,...)
   else if(meth=="exhaustive")
      lab.optimize.exhaustive(d1,d2,FUN,exchange.list,seek,...)
   else if(meth=="mc")
      lab.optimize.mc(d1,d2,FUN,exchange.list,seek,...)
   else if(meth=="hillclimb")
      lab.optimize.hillclimb(d1,d2,FUN,exchange.list,seek,...)
   else if(meth=="gumbel"){
      warning("Warning, gumbel method not yet supported. Try at your own risk.\n")
      lab.optimize.gumbel(d1,d2,FUN,exchange.list,seek,...)
   }
}


#lab.optimize.anneal - Annealing method for lab.optimize
lab.optimize.anneal<-function(d1,d2,FUN,exchange.list=0,seek="min",prob.init=1,prob.decay=0.99,freeze.time=1000,full.neighborhood=TRUE,...){
   #Pre-process the raw input data
   d1<-as.sociomatrix.sna(d1)
   d2<-as.sociomatrix.sna(d2)
   if(is.list(d1)||is.list(d2)||(dim(d1)[2]!=dim(d2)[2]))
     stop("lab.optimize routines require input graphs to be of identical order.")
   #End pre-processing
   #Find the data set size
   n<-dim(d1)[2]
   #If exchange list is a single number or vector, expand it via replication in a reasonable manner
   if(is.null(dim(exchange.list))){       #Exchange list was given as a single number or vector
      if(length(exchange.list)==1){                 #Single number case
         el<-matrix(rep(exchange.list,2*n),nrow=2,ncol=n)
      }else{                                                    #Vector case
         el<-sapply(exchange.list,rep,2)
      }  
   }else                         #Exchange list was given as a matrix; keep it.
      el<-exchange.list
   #Initialize various things
   fun<-match.fun(FUN)       #Find the function to be optimized
   d1<-d1[order(el[1,]),order(el[1,])]  #Reorder d1
   d2<-d2[order(el[2,]),order(el[2,])]  #Reorder d2
   el[1,]<-el[1,order(el[1,])]  #Reorder the exchange lists to match
   el[2,]<-el[2,order(el[2,])]
   if(any(el[1,]!=el[2,]))  #Make sure the exlist is legal
      stop("Illegal exchange list; lists must be comparable!\n")
   best<-fun(d1,d2,...)  #Take the seed value (this has to be legal)
   o<-1:n                              #Set the initial ordering
   global.best<-best              #Set global best values
   global.o<-o
   prob<-prob.init                  #Set acceptance prob
   ftime<-freeze.time            #Set time until freezing occurs
   nc<-choose(n,2)               #How many candidate steps?
   candp<-sapply(o,rep,choose(n,2))   #Build the candidate permutation matrix
   ccount<-1
   for(i in 1:n)
      for(j in i:n)
         if(i!=j){                                          #Perform binary exchanges 
            temp<-candp[ccount,i]
            candp[ccount,i]<-candp[ccount,j]
            candp[ccount,j]<-temp
            ccount<-ccount+1
         }         
   #Run the annealer
   flag<-FALSE
   if(any(duplicated(el[2,])))                #If we're dealing with the labeled case, don't bother.
      while((!flag)|(ftime>0)){           #Until we both freeze _and_ reach an optimum...
         #cat("Best: ",o," Perf: ",best," Global best: ",global.o," Global perf: ",global.best," Temp: ",prob,"\n")
         #cat("Perf: ",best," Global perf: ",global.best," Temp: ",prob,"\n")
         flag<-TRUE
         if(full.neighborhood){ #Full neighborhood search method - much slower, but more likely to find the optimum
            candperf<-vector()
            for(i in 1:nc)                           #Use candidate permutation matrix to produce new candidates
               if(all(el[2,]==el[2,o[candp[i,]]]))   #Is this legal?
                  candperf[i]<-fun(d1,d2[o[candp[i,]],o[candp[i,]]],...)
               else
                  candperf[i]<-NA   #If not, put the results in as missing data
            if(seek=="min"){
               bestcand<-(1:nc)[candperf==min(candperf,na.rm=TRUE)]   #Find the best candidate
               bestcand<-bestcand[!is.na(bestcand)]            
               if(length(bestcand)>1)
                  bestcand<-sample(bestcand,1)  #If we have multiple best candidates, choose one at random
               #cat(min(candperf,na.rm=TRUE),bestcand,candperf[bestcand],"\n")
               if(candperf[bestcand]<best){  #If this is better, move on and keep looking...        
                  o<-o[candp[bestcand,]]
                  best<-candperf[bestcand]
                  flag<-FALSE
                  if(best<global.best){   #Check to see if this is better than the global best
                     global.best<-best
                     global.o<-o
                  }
               }else if((ftime>0)&(runif(1,0,1)<prob)){   #...but if not frozen and no better option, take a chance.
                  #print(candperf)
                  bestcand<-sample(1:nc,1)               #Choose randomly from the available options
                  while(!all(el[2,]==el[2,o[candp[bestcand,]]]))  #Make sure we have a legal one...
                     bestcand<-sample(1:nc,1)
                  #cat("Wildcard - perm ",bestcand," perf ",candperf[bestcand],"\n")
                  o<-o[candp[bestcand,]]                #Accept the new candidate
                  best<-candperf[bestcand]
               }
            }else{
               bestcand<-(1:nc)[candperf==max(candperf,na.rm=TRUE)]   #Find the best candidate
               bestcand<-bestcand[!is.na(bestcand)]            
               if(length(bestcand)>1)
                  bestcand<-sample(bestcand,1)          #If we have multiple best candidates, choose one at random
               if((candperf[bestcand]>best)|(runif(1,0,1)<prob)){  #If this is better, move on and keep looking...        
                  o<-o[candp[bestcand,]]
                  best<-candperf[bestcand]
                  flag<-FALSE
                  if(best>global.best){   #Check to see if this is better than the global best
                     global.best<-best
                     global.o<-o
                  }
               }else if((ftime>0)&(runif(1,0,1)<prob)){   #...but if not frozen and no better option, take a chance.
                  bestcand<-sample(1:nc,1)               #Choose randomly from the available options
                  while(!all(el[2,]==el[2,o[candp[bestcand,]]]))  #Make sure we have a legal one...
                     bestcand<-sample(1:nc,1)
                  o<-o[candp[bestcand,]]                #Accept the new candidate
                  best<-candperf[bestcand]
               }
            }
         }else{      #Single candidate method.  Much faster, but less likely to find the optimum.
            #Use candidate permutation matrix to produce new candidates
            i<-sample(1:nc,1)       
            while(!all(el[2,]==el[2,o[candp[i,]]]))   #Is this legal?
               i<-sample(1:nc,1)       #Keep trying till we get it right.
            #Assess candidate performance
            candperf<-fun(d1,d2[o[candp[i,]],o[candp[i,]]],...)
            #Make a decision
            if(seek=="min"){
               if(candperf<best){  #If this is better, move on and keep looking...        
                  o<-o[candp[i,]]
                  best<-candperf
                  flag<-FALSE
                  if(best<global.best){   #Check to see if this is better than the global best
                     global.best<-best
                     global.o<-o
                  }
               }else if((ftime>0)&(runif(1,0,1)<prob)){   #...but if not frozen and no better option, take a chance.
                  i<-sample(1:nc,1)               #Choose randomly from the available options
                  while(!all(el[2,]==el[2,o[candp[i,]]]))  #Make sure we have a legal one...
                     i<-sample(1:nc,1)
                  o<-o[candp[i,]]                #Accept the new candidate
                  best<-candperf
                  if(best<global.best){   #Check to see if this is better than the global best
                     global.best<-best
                     global.o<-o
                  }
               }
            }else{
               if(candperf>best){  #If this is better, move on and keep looking...        
                  o<-o[candp[i,]]
                  best<-candperf
                  flag<-FALSE
                  if(best>global.best){   #Check to see if this is better than the global best
                     global.best<-best
                     global.o<-o
                  }
               }else if((ftime>0)&(runif(1,0,1)<prob)){   #...but if not frozen and no better option, take a chance.
                  i<-sample(1:nc,1)               #Choose randomly from the available options
                  while(!all(el[2,]==el[2,o[candp[i,]]]))  #Make sure we have a legal one...
                     i<-sample(1:nc,1)
                  o<-o[candp[i,]]                #Accept the new candidate
                  best<-candperf
                  if(best>global.best){   #Check to see if this is better than the global best
                     global.best<-best
                     global.o<-o
                  }
               }
            }
         }
         #Set things up for the next iteration (if there is one)
         ftime<-ftime-1               #Continue the countdown to the freezing point
         prob<-prob*prob.decay   #Cool things off a bit
      }
   #Report the results
   global.best
}


#lab.optimize.exhaustive - Exhaustive search method for lab.optimize
lab.optimize.exhaustive<-function(d1,d2,FUN,exchange.list=0,seek="min",...){
   #Pre-process the raw input data
   d1<-as.sociomatrix.sna(d1)
   d2<-as.sociomatrix.sna(d2)
   if(is.list(d1)||is.list(d2)||(dim(d1)[2]!=dim(d2)[2]))
     stop("lab.optimize routines require input graphs to be of identical order.")
   #End pre-processing
   #Find the data set size
   n<-dim(d1)[2]
   #If exchange list is a single number or vector, expand it via replication in a reasonable manner
   if(is.null(dim(exchange.list))){       #Exchange list was given as a single number or vector
      if(length(exchange.list)==1){                 #Single number case
         el<-matrix(rep(exchange.list,2*n),nrow=2,ncol=n)
      }else{                                                    #Vector case
         el<-sapply(exchange.list,rep,2)
      }  
   }else                         #Exchange list was given as a matrix; keep it.
      el<-exchange.list
   #Initialize various things
   fun<-match.fun(FUN)       #Find the function to be optimized
   d1<-d1[order(el[1,]),order(el[1,])]  #Reorder d1
   d2<-d2[order(el[2,]),order(el[2,])]  #Reorder d2
   el[1,]<-el[1,order(el[1,])]  #Reorder the exchange lists to match
   el[2,]<-el[2,order(el[2,])]
   if(any(el[1,]!=el[2,]))  #Make sure the exlist is legal
      stop("Illegal exchange list; lists must be comparable!\n")
   best<-fun(d1,d2,...)  #Take the seed value (this has to be legal)
   #Search exhaustively - I hope you're not in a hurry!
   if(any(duplicated(el[1,])))  #If we're dealing with the labeled case, don't bother.
      for(k in 0:(gamma(n+1)-1)){
         o<-numperm(n,k)
         if(all(el[1,]==el[2,o])){
            if(seek=="min")
               best<-min(best,fun(d1,d2[o,o],...))
            else
               best<-max(best,fun(d1,d2[o,o],...))
         }
      }
   #Report the results
   best
}


#lab.optimize.gumbel - Extreme value method for lab.optimize
lab.optimize.gumbel<-function(d1,d2,FUN,exchange.list=0,seek="min",draws=500,tol=1e-5,estimator="median",...){
   #Pre-process the raw input data
   d1<-as.sociomatrix.sna(d1)
   d2<-as.sociomatrix.sna(d2)
   if(is.list(d1)||is.list(d2)||(dim(d1)[2]!=dim(d2)[2]))
     stop("lab.optimize routines require input graphs to be of identical order.")
   #End pre-processing
   #Find the data set size
   n<-dim(d1)[2]
   #If exchange list is a single number or vector, expand it via replication in a reasonable manner
   if(is.null(dim(exchange.list))){       #Exchange list was given as a single number or vector
      if(length(exchange.list)==1){                 #Single number case
         el<-matrix(rep(exchange.list,2*n),nrow=2,ncol=n)
      }else{                                                    #Vector case
         el<-sapply(exchange.list,rep,2)
      }  
   }else                         #Exchange list was given as a matrix; keep it.
      el<-exchange.list
   #Initialize various things
   fun<-match.fun(FUN)       #Find the function to be optimized
   fg<-vector()                    #Set up the function
   d1<-d1[order(el[1,]),order(el[1,])]  #Reorder d1
   d2<-d2[order(el[2,]),order(el[2,])]  #Reorder d2
   el[1,]<-el[1,order(el[1,])]  #Reorder the exchange lists to match
   el[2,]<-el[2,order(el[2,])]
   if(any(el[1,]!=el[2,]))  #Make sure the exlist is legal
      stop("Illegal exchange list; lists must be comparable!\n")
   #Approximate the distribution using Monte Carlo
   for(i in 1:draws){
      o<-rperm(el[2,])
      fg[i]<-fun(d1,d2[o,o],...)
   }
   #Use the approximated distribution to fit a Gumbel model for the extreme values;
   #this is only approximate, since the extreme value model assumes an unbounded, continuous underlying
   #distribution.  Also, these results are "unproven," in the sense that no actual permutation has been
   #found by the algorithm which results in the predicted value (unlike the other methods); OTOH, in 
   #a world of approximations, this one may not be any worse than the others....
   b<-1
   b.old<-1
   bdiff<-Inf
   mfg<-mean(fg)
   print(quantile(fg))
   while(bdiff>tol){         #Solve iteratively for bhat
      cat("bold=",b.old,"b=",b,"bdiff=",bdiff,"\n")
      b.old<-b
      b<-mfg-sum(fg*exp(-fg/b))/sum(exp(-fg/b))
      bdiff<-abs(b.old-b)
   }
   a<--b*log(sum(exp(-fg/b))/draws)     #Given this, ahat is a function of bhat and the data
   #Report the results
   cat("a=",a,"b=",b,"\n")
   switch(estimator,
      mean=a-b*digamma(1),
      mode=a,
      median=a-b*log(log(2))
   )
}


#lab.optimize.hillclimb - Hill-climbing method for lab.optimize
lab.optimize.hillclimb<-function(d1,d2,FUN,exchange.list=0,seek="min",...){
   #Pre-process the raw input data
   d1<-as.sociomatrix.sna(d1)
   d2<-as.sociomatrix.sna(d2)
   if(is.list(d1)||is.list(d2)||(dim(d1)[2]!=dim(d2)[2]))
     stop("lab.optimize routines require input graphs to be of identical order.")
   #End pre-processing
   #Find the data set size
   n<-dim(d1)[2]
   #If exchange list is a single number or vector, expand it via replication in a reasonable manner
   if(is.null(dim(exchange.list))){       #Exchange list was given as a single number or vector
      if(length(exchange.list)==1){                 #Single number case
         el<-matrix(rep(exchange.list,2*n),nrow=2,ncol=n)
      }else{                                                    #Vector case
         el<-sapply(exchange.list,rep,2)
      }  
   }else                         #Exchange list was given as a matrix; keep it.
      el<-exchange.list
   #Initialize various things
   fun<-match.fun(FUN)       #Find the function to be optimized
   d1<-d1[order(el[1,]),order(el[1,])]  #Reorder d1
   d2<-d2[order(el[2,]),order(el[2,])]  #Reorder d2
   el[1,]<-el[1,order(el[1,])]  #Reorder the exchange lists to match
   el[2,]<-el[2,order(el[2,])]
   if(any(el[1,]!=el[2,]))  #Make sure the exlist is legal
      stop("Illegal exchange list; lists must be comparable!\n")
   best<-fun(d1,d2,...)  #Take the seed value (this has to be legal)
   o<-1:n                              #Set the initial ordering
   nc<-choose(n,2)               #How many candidate steps?
   candp<-sapply(o,rep,choose(n,2))   #Build the candidate permutation matrix
   ccount<-1
   for(i in 1:n)
      for(j in i:n)
         if(i!=j){                                          #Perform binary exchanges 
            temp<-candp[ccount,i]
            candp[ccount,i]<-candp[ccount,j]
            candp[ccount,j]<-temp
            ccount<-ccount+1
         }         
   #Run the hill climber
   flag<-FALSE
   while(!flag){           #Until we reach an optimum...
      #cat("Best: ",o," Perf: ",best,"\n")
      flag<-TRUE
      candperf<-vector()
      for(i in 1:nc)                           #Use candidate permutation matrix to produce new candidates
         if(all(el[2,]==el[2,o[candp[i,]]]))   #Is this legal?
            candperf[i]<-fun(d1,d2[o[candp[i,]],o[candp[i,]]],...)
         else
            candperf[i]<-NA   #If not, put the results in as missing data
      if(seek=="min"){
         bestcand<-(1:nc)[candperf==min(candperf,na.rm=TRUE)]   #Find the best candidate
         if(length(bestcand)>1)
            bestcand<-sample(bestcand,1)          #If we have multiple best candidates, choose one at random
         if(candperf[bestcand]<best){          #If this is better, move on and keep looking...        
            o<-o[candp[bestcand,]]
            best<-candperf[bestcand]
            flag<-FALSE
         }
      }else{
         bestcand<-(1:nc)[candperf==max(candperf,na.rm=TRUE)]   #Find the best candidate
         if(length(bestcand)>1)
            bestcand<-sample(bestcand,1)          #If we have multiple best candidates, choose one at random
         if(candperf[bestcand]>best){          #If this is better, move on and keep looking...        
            o<-o[candp[bestcand,]]
            best<-candperf[bestcand]
            flag<-FALSE
         }
      }
   }
   #Report the results
   best
}


#lab.optimize.mc - Monte Carlo method for lab.optimize
lab.optimize.mc<-function(d1,d2,FUN,exchange.list=0,seek="min",draws=1000,...){
   #Pre-process the raw input data
   d1<-as.sociomatrix.sna(d1)
   d2<-as.sociomatrix.sna(d2)
   if(is.list(d1)||is.list(d2)||(dim(d1)[2]!=dim(d2)[2]))
     stop("lab.optimize routines require input graphs to be of identical order.")
   #End pre-processing
   #Find the data set size
   n<-dim(d1)[2]
   #If exchange list is a single number or vector, expand it via replication in a reasonable manner
   if(is.null(dim(exchange.list))){       #Exchange list was given as a single number or vector
      if(length(exchange.list)==1){                 #Single number case
         el<-matrix(rep(exchange.list,2*n),nrow=2,ncol=n)
      }else{                                                    #Vector case
         el<-sapply(exchange.list,rep,2)
      }  
   }else                         #Exchange list was given as a matrix; keep it.
      el<-exchange.list
   #Initialize various things
   fun<-match.fun(FUN)       #Find the function to be optimized
   d1<-d1[order(el[1,]),order(el[1,])]  #Reorder d1
   d2<-d2[order(el[2,]),order(el[2,])]  #Reorder d2
   el[1,]<-el[1,order(el[1,])]  #Reorder the exchange lists to match
   el[2,]<-el[2,order(el[2,])]
   if(any(el[1,]!=el[2,]))  #Make sure the exlist is legal
      stop("Illegal exchange list; lists must be comparable!\n")
   best<-fun(d1,d2,...)  #Take the seed value (this has to be legal)
   #Search via blind monte carlo - slow, yet ineffectual
   if(any(duplicated(el[1,])))  #If we're dealing with the labeled case, don't bother.
      for(i in 1:draws){
         o<-rperm(el[2,])
         if(seek=="min")
            best<-min(best,fun(d1,d2[o,o],...))
         else
            best<-max(best,fun(d1,d2[o,o],...))
      }
   #Report the results
   best
}


#numperm - Get the nth permutation vector by periodic placement
numperm<-function(olength,permnum){
   if((permnum>gamma(olength+1)-1)|(permnum<0)){
      cat("permnum must be an integer in [0,olength!-1]\n")
   }
   o<-vector(length=olength)
   o[]<--1
   pnum<-permnum
   for(i in 1:olength){
      relpos<-pnum%%(olength-i+1)
      flag<-FALSE
      p<-1
      while(!flag)
         if(o[p]==-1){
            if(relpos==0){
               o[p]<-i
               flag<-TRUE
            }else{
               p<-p+1
               relpos<-relpos-1
            }
         }else
            p<-p+1      
      pnum<-pnum%/%(olength-i+1)
   }
   o
}


#rmperm - Randomly permutes the rows and columns of an input matrix.
rmperm<-function(m){
   #Pre-process the raw input
   m<-as.sociomatrix.sna(m)
   if(is.list(m))
     return(lapply(m,rmperm))
   #End pre-processing
   if(length(dim(m))==2){
      #Only a single matrix is included
      o<-sample(1:dim(m)[1])
      p<-matrix(data=m[o,o],nrow=dim(m)[1],ncol=dim(m)[2])
   }else{
      #Here, we assume a stack of matrices
      p<-array(dim=c(dim(m)[1],dim(m)[2],dim(m)[3]))
      for(i in 1:dim(m)[1]){
         o<-sample(1:dim(m)[2])
         p[i,,]<-array(m[i,o,o])
      }
   }
   p
}


#rperm - Draw a random permutation vector with exchangability constraints
rperm<-function(exchange.list){
   #Note that exchange.list should be a vector whose entries correspond to the class identity
   #of the respective element.  It doesn't matter what the values are, so long as elements have
   #the same value iff they are exchangeable.
   n<-length(exchange.list)   #Get the length of the output vector
   grp<-unique(exchange.list) #Get the groups
   o<-1:n  #Create the initial ordering
   #Randomly scramble orders within groups
   for(i in grp){
      v<-(1:n)[exchange.list==i]
      if(length(v)>1)      #Need this test, because sample is too smart for its own good...
         o[v]<-sample(v)
   }
   #Return the permutation
   o
}

