# TODO: Add comment
# 
# Author: nmorris
###############################################################################
# BASIC HAMILTONIAN MONTE CARLO UPDATE
# This HMC implementation has been greatly modified from the one created by Radford M. Neal, 2012.
# but may still bear a vague simularity

#The copyright for the original code is described bellow:

#>  This directory contains a preliminary version of GRIMS, released 2012-06-07.
#>  
#>  The contents of this directory are Copyright (c) 2011-2012 by Radford M. Neal.
#> 
#>  Permission is granted for anyone to copy, use, modify, or distribute these
#>  programs and accompanying documents for any purpose, provided this copyright
#>  notice is retained and prominently displayed, along with a note saying 
#>  that the original programs are available from Radford Neal's web page, and 
#>    note is made of any changes made to these programs.  These programs and 
#>    documents are distributed without any warranty, express or implied.  As the
#>    programs were written for research purposes only, they have not been tested 
#>    to the degree that would be advisable in any important application.  All use
#>    of these programs is entirely at the user's own risk.
#> 
#>  This software can be obtained from http://www.cs.utoronto.ca/~radford/GRIMS.html

setClass("lpgr",
         representation(
           gr="list",
           lp = "numeric"),
         prototype = list(
           gr=NULL,
           lp=NULL)
) 




.basicHMCOneStep <- function (lpgrf, initial, lpgr.initial, nsteps, step,
                      initial.p = lapply(initial,function(x) gmatrix(grnorm(nrow(x)*ncol(x)), nrow(x),ncol(x), dup=FALSE)),
                      T=1, errorChecks=TRUE)
{
  
	if(errorChecks) {
		#Check and process the arguments.
		if(!("list" %in% class(initial) && "list" %in% class(initial.p)))
		stop("initial and initial.p must be lists")
		nm1=names(initial)
		nm2=names(initial.p)
		if(length(nm1)!=length(nm2))
			stop("initial and initial.p do not match.")
		if(!all(nm1==nm2))
			stop("initial and initial.p do not match.")
		if(!all(sapply(initial, function(x) class(x) %in% c("matrix","gmatrix"))))
		  stop("All elements of 'initial' must be either of class matrix or 'gmatrix.'")
		if(!all(sapply(initial.p, function(x) class(x) %in% c("matrix","gmatrix"))))
		  stop("All elements of 'initial.p' must be either of class matrix or 'gmatrix.'")
		
		nparrallel=ncol(initial[[1]])
		if(!all(sapply(initial, ncol)==nparrallel))
		  stop("'initial' has non matching number of cols")
		if(!all(sapply(initial, ncol)==nparrallel))
		  stop("'initial.p' has non matching number of cols.")
		
		if(!all(sapply(initial, nrow)==sapply(initial.p, nrow)))
		  stop("'initial.p' vs 'initial' has non matching number of rows")
        
		nsteps=as.integer(nsteps)[1]
		if(length(step)==1){
		  step=as.numeric(step)
		  step=sapply(nm1, function(i) step)
		}
		nm3=names(step)
		if(length(nm1)!=length(nm3))
		  stop("initial and step do not match.")
		if(!all(nm1==nm3))
		  stop("initial and step do not match.")
	}
	if(missing(lpgr.initial))
		lpgr.initial <- lpgrf(initial)

	if(errorChecks) {
		if(!("lpgr" %in% class(lpgr.initial))) {
		  stop("lpgr.initial not of class 'lpgr.' Please make sure the lpgr function returns an object of the correct class.")
		}
		if(length(lpgr.initial @lp)!=nparrallel)
		  stop("The lpgr is not returning the the log probability (lp) for the correct numer of parrallel simulations.")
		nm4=names(lpgr.initial@gr)
		if(length(nm1)!=length(nm4))
		  stop("initial names and lpgr function output (for the 'gr' slot) do not match.")
		if(!all(nm1==nm4))
		  stop("initial names and lpgr function output (for the 'gr' slot) do not match.")
		if(!all(sapply(lpgr.initial@gr, ncol)==nparrallel))
		  stop("Number of columns returned by lpgr function in an element of the 'gr' slot is not the parrallel number of simulations.")
	}
	
	# Compute the kinetic energy at the start of the trajectory.
	kinetic.initial <- rowSums(sapply(initial.p,function(x) { x=as.matrix(x); colSums(x*x) / 2 }))
  
	# Compute the trajectory by the leapfrog method.
	q <- initial
	p <- initial.p

	# Make a half step for momentum at the beginning.
	p <- mapply(function(p,gr,step)  return(p + (step/2) * gr), p , lpgr.initial@gr, step, SIMPLIFY=FALSE)

	# Alternate full steps for position and momentum.
	for (i in 1:nsteps)
	{ 
		# Make a full step for the position, and evaluate the gradient at the new 
		# position.
		q <- mapply(function(p,q,step) return( q + step * p), p , q, step, SIMPLIFY = FALSE)	 #q <- q + step * p   
		lpgr.current <- lpgrf(q) 
		if(!("lpgr" %in% class(lpgr.initial))) {
		  stop("lpgr not of class 'lpgr.'")
		}
		# Make a full step for the momentum, except when we're coming to the end of 
		# the trajectory.  		
		if (i!=nsteps)
		{ #p <- p + step * gr  
			p=mapply(function(p,gr,step)  return(p + (step) * gr), p , lpgr.current@gr, step, SIMPLIFY = FALSE)
		}
	}
	
	# Make a half step for momentum at the end.	
	#p <- p + (step/2) * gr
	p <- mapply(function(p,gr,step)  return(p + (step/2) * gr), p , lpgr.current@gr, step, SIMPLIFY = FALSE)
	
	
	# Look at log probability and kinetic energy at the end of the trajectory.
	#lpgr.prop <- lpgr.current$ll
	#kinetic.prop <- sum(p^2) / 2
	kinetic.prop <-  rowSums(sapply(p,function(x) { x=as.matrix(x); colSums(x*x) / 2 }))
	
	# Accept or reject the state at the end of the trajectory.
	H.initial <- -as.numeric(lpgr.initial@lp) + kinetic.initial
	H.prop <- -as.numeric( lpgr.current@lp )+ kinetic.prop
	delta <- H.prop - H.initial
	apr <- pmin(1,exp(-delta/T))
	acc <- runif(nparrallel) < apr
	notacc=which(!acc)
    
	moveBack=function(propl,oldl) {
		mapply(function(prop,old) {
			prop[,notacc]=old[, notacc]
			return(prop)},
		propl,oldl,SIMPLIFY = FALSE)
	}
	if(length(notacc) > 0) {
		q=moveBack(q, initial)
		lpgr.current@gr=moveBack(lpgr.current@gr, lpgr.initial@gr)
		lpgr.current@lp[notacc] = lpgr.initial@lp[notacc]
	}
	return(list (acc=acc, q.current=q, lpgr.current=lpgr.current))
}

.simp=function(lst) {
	return(lapply(lst, function(x) 
		if("gmatrix" %in% class(x))
			return(gmatrix(grnorm(nrow(x)*ncol(x)), nrow(x),ncol(x), dup=FALSE))
		else
			return(matrix(rnorm(nrow(x)*ncol(x)), nrow(x),ncol(x)))
	))
}

keep=function(q) lapply(q, function(x) if(any(class(x) %in% c("gmatrix","gvector"))) h(x) else x)

gBasicHMC = function (lpgrf, initial, nsims, nsteps, step, 
                     burnin=1, nstepsburnin=nsteps,stepburnin=step, Tstart=1, r=1 ,
                     keep=keep,
					 thin=1, report=100){

  sims=list()
  lp=list()
  st=list(q.current=initial, lpgr.current=lpgrf(initial))
  psims=length(st$lpgr.current@lp)
  accBurninSum=numeric(psims)
  accSimsSum=numeric(psims)
  accBurninDenom=0
  accSimsDenom=0
  
  curnsteps=nstepsburnin
  curstep=stepburnin
  T=Tstart
  for(i in 1:nsims) {
    initial.p = .simp(st$q.current)
    st=.basicHMCOneStep(lpgrf, 
                        initial=st$q.current,
                        initial.p=initial.p,
                        lpgr.initial=st$lpgr.current,
                        nsteps=curnsteps, step=curstep, T=T)
	lp[[i]]=st$lpgr.current@lp
    if(i<burnin) {
      T = max(1, r * T)
	  
      accBurninSum=accBurninSum+st$acc
      accBurninDenom=accBurninDenom+1

      
    } else {
      T=1
      curnsteps=nsteps
      curstep=step
      accSimsSum=accSimsSum+st$acc
      accSimsDenom=accSimsDenom+1
      if( (i %% thin)==0) {
        sims[[paste0("Sim",i)]]=keep(st$q.current)
      }
    }
    
    if( (i %% report)==0) {
		cat("i =", i, "\n")
		if(i<burnin) { 
			cat("   Burnin Acceptance Temperature:", T, "\n")
			cat("   Burnin Acceptance rates:", accBurninSum/accBurninDenom,"\n")
		} else {
			cat("   Acceptance rates:", accSimsSum/accSimsDenom,"\n")
		}
		cat("   Current log proability density (not necisarily normalized):", st$lpgr.current@lp,"\n")
    }
    
  }
  return(list(sims=sims,
			  lp=do.call(rbind, lp),
              AcceptanceRate=accSimsSum/accSimsDenom, 
              BurninAcceptanceRate=accBurninSum/accBurninDenom))
}
