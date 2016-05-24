# Author: Maarten Speekenbrink
# With suggestions from Robert McGehee for speeding up the code for large objects
# First version: 23-3-2008
# latest change:

rdirichlet <- function(n, alpha) {
  # taken from gtools
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    x/as.vector(sm)
}

which.is.max <- function(x) {
    # taken from MASS
    y <- seq_along(x)[x == max(x)]
    if(length(y) > 1L) sample(y, 1L) else y
}


ind.max <- function(x) {
    out <- rep(0,length(x))
    out[which.is.max(x)] <- 1
    out
}


emviterbi <- function(A,B,init,ntimes,nstates,homogeneous,na.allow=TRUE) {
    # used for EM with hard classification, so that we don't need to change the object...
    # returns the most likely state sequence
    nt <- sum(ntimes)
    lt <- length(ntimes)
    et <- cumsum(ntimes)
    bt <- c(1,et[-lt]+1)
    
    ns <- nstates
    
    delta <- psi <- matrix(nrow=nt,ncol=ns)
    state <- vector(length=nt)
    pstate <- vector(length=nt)
    
    prior <- init
    
    #if(max(ntimes>1)) A <- object@trDens
    #B <- object@dens
    if(na.allow) B <- replace(B,is.na(B),1)
    B <- apply(B,c(1,3),prod)
    
    for(case in 1:lt) {
      # initialization
      delta[bt[case],] <- prior[case,]*B[bt[case],]
      delta[bt[case],] <- delta[bt[case],]/(sum(delta[bt[case],]))
      psi[bt[case],] <- 0
      # recursion
      if(ntimes[case]>1) {
        for(tt in ((bt[case]+1):et[case])) {
          for(j in 1:ns) {
            if(!homogeneous) {
              delta[tt,j] <- max(delta[tt-1,]*(A[tt,j,]))*B[tt,j]
              k <- which.max(delta[tt-1,]*A[tt,j,])
            } else {
              delta[tt,j] <- max(delta[tt-1,]*(A[1,j,]))*B[tt,j]
              k <- which.max(delta[tt-1,]*A[1,j,])
            }
            if(length(k) == 0) k <- 0 # what's this doing here??? can this ever occur? FIX ME
            psi[tt,j] <- k
          }
          delta[tt,] <- delta[tt,]/(sum(delta[tt,]))
          
        }
      }
      # trace maximum likely state
      state[et[case]] <- which.max(delta[et[case],])
      # this doesn't need a for loop does it???? FIX ME
      if(ntimes[case]>1) {
        for(i in (et[case]-1):bt[case]) {
          state[i] <- psi[i+1,state[i+1]]
        }
      }
    }
    # compute the unconditional probability of the state sequence
    pstate[bt] <- prior[cbind(1:lt,state[bt])]
    btt <- bt[ntimes>1]
    ett <- et[ntimes>1]
    idx <- unlist(mapply(seq,btt+1,ett))
    if(!homogeneous) {
      pstate[idx] <- A[cbind(idx,state[idx],state[idx-1])]
    } else {
      pstate[idx] <- A[cbind(1,state[idx],state[idx-1])]
    }
    delta <- data.frame(state,pstate,delta)
    return(delta)
}


em <- function(object,...) {
	if(!is(object,"mix")) stop("object is not of class '(dep)mix'")
	call <- match.call()
	if(is(object,"depmix")) {
		call[[1]] <- as.name("em.depmix")
	} else {
		call[[1]] <- as.name("em.mix")
	}
	object <- eval(call, parent.frame())
	object
}

# em for lca and mixture models
em.mix <- function(object,maxit=100,tol=1e-8,crit=c("relative","absolute"),random.start=TRUE,verbose=FALSE,classification=c("soft","hard"),na.allow=TRUE,...) {
	
	clsf <- match.arg(classification)
	crit <- match.arg(crit)
	
	if(!is(object,"mix")) stop("object is not of class 'mix'")
		
	ns <- nstates(object)
	ntimes <- ntimes(object)
	lt <- length(ntimes)
	et <- cumsum(ntimes)
	bt <- c(1,et[-lt]+1)
  
	prior <- object@prior
	response <- object@response
	dens <- object@dens
	init <- dens(object@prior)

	converge <- FALSE
	
	if(random.start) {
		nr <- sum(ntimes(object))
		gamma <- rdirichlet(nr,alpha=rep(.01,ns))
		if(clsf == "hard") {
		    gamma <- t(apply(gamma,1,ind.max))
		}
		LL <- -1e10
		
		for(i in 1:ns) {
			for(k in 1:nresp(object)) {
			  response[[i]][[k]] <- fit(response[[i]][[k]],w=gamma[,i])
				# update dens slot of the model
				dens[,k,i] <- dens(response[[i]][[k]])
			}
		}
		# initial expectation
		if(clsf == "hard") {
		  fbo <- list()
		  vstate <- apply(gamma,1,which.max)
      fbo$gamma <- t(apply(gamma,1,ind.max))
		  B <- dens
		  if(na.allow) B[is.na(B)] <- 1
		  #fbo$gamma <- t(apply(gamma,1,ind.max))
		  fbo$logLike <- sum(log((apply(B,c(1,3),prod))[cbind(1:sum(ntimes),vstate)])) + sum(log(init[cbind(1:lt,vstate)]))
		} else {
		  fbo <- fb(init=init,matrix(0,1,1),B=dens,ntimes=ntimes(object))
		}
		LL <- fbo$logLike
		
		if(is.nan(LL)) stop("Cannot find suitable starting values; please provide them.")
		
	} else {
		# initial expectation
		fbo <- fb(init=init,A=matrix(0,1,1),B=dens,ntimes=ntimes(object))
		if(clsf == "hard") {
      fbo$gamma <- t(apply(fbo$gamma,1,ind.max))
		  vstate <- apply(fbo$gamma,1,which.max)
		  B <- object@dens
		  if(na.allow) B[is.na(B)] <- 1
		  fbo$logLike <- sum(log((apply(B,c(1,3),prod))[cbind(1:sum(ntimes),vstate)])) + sum(log(init[cbind(1:lt,vstate)]))
          
		}
		LL <- fbo$logLike
		if(is.nan(LL)) stop("Starting values not feasible; please provide them.")
	}
	
	LL.old <- LL + 1
	
	for(j in 0:maxit) {
		
		# maximization		
		prior@y <- fbo$gamma[bt,,drop=FALSE]
		prior <- fit(prior, w=NULL,ntimes=NULL)
		init <- dens(prior)
		
		for(i in 1:ns) {
			for(k in 1:nresp(object)) {
				if(sum(fbo$gamma[,i]) > 0) {
          response[[i]][[k]] <- fit(response[[i]][[k]],w=fbo$gamma[,i])
				  # update dens slot of the model
				  dens[,k,i] <- dens(response[[i]][[k]])
				}
			}
		}
		
		# expectation
		fbo <- fb(init=init,A=matrix(0,1,1),B=dens,ntimes=ntimes(object))
		if(clsf == "hard") {
		  fbo$gamma <- t(apply(fbo$gamma,1,ind.max))
		  vstate <- apply(fbo$gamma,1,which.max)
		  B <- dens
		  if(na.allow) B[is.na(B)] <- 1
		  fbo$logLike <- sum(log((apply(B,c(1,3),prod))[cbind(1:sum(ntimes),vstate)])) + sum(log(init[cbind(1:lt,vstate)]))
		}
		LL <- fbo$logLike
		# print stuff
		if(verbose&((j%%5)==0)) {
			cat("iteration",j,"logLik:",LL,"\n")
		}
		
		if(LL >= LL.old) {
		  converge <- (crit == "absolute" &&  LL - LL.old < tol) || (crit == "relative" && (LL - LL.old)/abs(LL.old)  < tol)
      if(converge) {
			  cat("converged at iteration",j,"with logLik:",LL,"\n")
			  break
			}
		} else {
			# this should not really happen...
			if(j > 0 && (LL.old - LL) > tol) stop("likelihood decreased on iteration ",j)
		}

		LL.old <- LL
	}
  
	object@prior <- prior
	object@init <- init
	object@response <- response
	object@dens <- dens

	if(clsf == "hard") {
	    object <- as(object,"mix.fitted.classLik") # class(object) <- "mix.fitted.classLik"
	    object@posterior <- data.frame(state=viterbi(object)[,1])
	} else {
	    object <- as(object,"mix.fitted") # class(object) <- "mix.fitted"
	    object@posterior <- viterbi(object)
	}

	if(converge) {
	  if(clsf == "hard") {
		  object@message <- switch(crit,
			  relative = "Log classification likelihood converged to within tol. (relative change)",
			  absolute = "Log classification likelihood converged to within tol. (absolute change)"
		  )	  
	  } else {
		  object@message <- switch(crit,
			  relative = "Log likelihood converged to within tol. (relative change)",
			  absolute = "Log likelihood converged to within tol. (absolute change)"
		  )
		}
	} else object@message <- "'maxit' iterations reached in EM without convergence."

	# no constraints in EM, except for the standard constraints ...
	# which are produced by the following (only necessary for getting df right in logLik and such)
	constraints <- getConstraints(object)
	object@conMat <- constraints$lincon
	object@lin.lower <- constraints$lin.l
	object@lin.upper <- constraints$lin.u
	
	object
	
}

# em for hidden markov models
em.depmix <- function(object,maxit=100,tol=1e-8,crit=c("relative","absolute"),random.start=TRUE,verbose=FALSE,classification=c("soft","hard"),na.allow=TRUE,...) {
	
	if(!is(object,"depmix")) stop("object is not of class 'depmix'")
	
	clsf <- match.arg(classification)
	crit <- match.arg(crit)
  
	ns <- nstates(object)
	
	ntimes <- ntimes(object)
	lt <- length(ntimes)
	et <- cumsum(ntimes)
	bt <- c(1,et[-lt]+1)
  
	prior <- object@prior
	transition <- object@transition
	trDens <- object@trDens
	response <- object@response
	dens <- object@dens
	init <- dens(object@prior)
	
	if(random.start) {
				
		nr <- sum(ntimes(object))
		gamma <- rdirichlet(nr,alpha=rep(.01,ns))
		if(clsf == "hard") {
		    gamma <- t(apply(gamma,1,ind.max))
		}
		LL <- -1e10
		
		for(i in 1:ns) {
			for(k in 1:nresp(object)) {
				response[[i]][[k]] <- fit(response[[i]][[k]],w=gamma[,i])
				# update dens slot of the model
				dens[,k,i] <- dens(response[[i]][[k]])
			}
		}
	}

	# initial expectation
  if(clsf == "hard") {
    fbo <- list()
    vit <- emviterbi(A=trDens,B=dens,init=init,ntimes=object@ntimes,nstates=ns,homogeneous=object@homogeneous,na.allow=na.allow)
	  vstate <- vit[,1]
    pstate <- vit[,2]
	  fbo$gamma <- as.matrix(model.matrix(~ factor(vstate,levels=1:ns) - 1))
	  fbo$xi <- array(0,dim=c(sum(ntimes),ns,ns))
	  fbo$xi[cbind(1:(sum(ntimes)- 1),vstate[-1],vstate[-length(vstate)])] <- 1
	  B <- dens
	  if(na.allow) B[is.na(B)] <- 1
	  fbo$logLike <- sum(log((apply(B,c(1,3),prod))[cbind(1:sum(ntimes),vstate)])) + sum(log(pstate))
	} else {
	  fbo <- fb(init=init,A=trDens,B=dens,ntimes=ntimes(object),homogeneous=object@homogeneous)
  }
	LL <- fbo$logLike
	if(is.nan(LL)) stop("Starting values not feasible; please provide them.")
	LL.old <- LL + 1 # force the "old" likelihood to be larger...
	
	converge <- FALSE
	for(j in 0:maxit) {
		
		# maximization
		prior@y <- fbo$gamma[bt,,drop=FALSE]
		prior <- fit(prior, w=NULL, ntimes=NULL)
		init <- dens(prior)
				
		trm <- matrix(0,ns,ns)
		for(i in 1:ns) {
			if(!object@homogeneous) {
        # TODO: check whether fbo$gamma > 0, otherwise set to previous value....
				transition[[i]]@y <- fbo$xi[,,i]/fbo$gamma[,i]
				transition[[i]] <- fit(transition[[i]],w=as.matrix(fbo$gamma[,i]),ntimes=ntimes(object)) # check this
			} else {
			  if(sum(fbo$gamma[-c(et),i]) == 0) {
          # set unidentified transition probs to previous value
          trm[i,] <- trDens[1,,i]
			  } else {
  				for(k in 1:ns) {
  					trm[i,k] <- sum(fbo$xi[-c(et),k,i])/sum(fbo$gamma[-c(et),i])
  				}
			  }
				# FIX THIS; it will only work with specific trinModels
				# should become object@transition = fit(object@transition, xi, gamma)
				transition[[i]]@parameters$coefficients <- switch(transition[[i]]@family$link,
					identity = transition[[i]]@family$linkfun(trm[i,]),
					mlogit = transition[[i]]@family$linkfun(trm[i,],base=transition[[i]]@family$base),
					transition[[i]]@family$linkfun(trm[i,])
				)
			}
			# update trDens slot of the model
			trDens[,,i] <- dens(transition[[i]])
		}
		
		for(i in 1:ns) {
			for(k in 1:nresp(object)) {
				if(sum(fbo$gamma[,i])>0) {
          response[[i]][[k]] <- fit(response[[i]][[k]],w=fbo$gamma[,i])
  				# update dens slot of the model
  				dens[,k,i] <- dens(response[[i]][[k]])
				}
			}
		}
		
		if(clsf == "hard") {
      fbo <- list()
      vit <- emviterbi(A=trDens,B=dens,init=init,ntimes=object@ntimes,nstates=ns,homogeneous=object@homogeneous,na.allow=na.allow)
      vstate <- vit[,1]
      pstate <- vit[,2]
		  #vstate <- viterbi(object)[,1]
		  fbo$gamma <- as.matrix(model.matrix(~ factor(vstate,levels=1:ns) - 1))
		  fbo$xi <- array(0,dim=c(sum(ntimes),ns,ns))
		  fbo$xi[cbind(1:(sum(ntimes)- 1),vstate[-1],vstate[-length(vstate)])] <- 1
		  B <- dens
		  if(na.allow) B[is.na(B)] <- 1
		  fbo$logLike <- sum(log((apply(B,c(1,3),prod))[cbind(1:sum(ntimes),vstate)])) + sum(log(pstate))
		} else {
		  # expectation
		  fbo <- fb(init=init,A=trDens,B=dens,ntimes=ntimes(object),homogeneous=object@homogeneous)	  
	  }
	  
	  LL <- fbo$logLike	
		

		if( (LL >= LL.old)) {
		  converge <- (crit == "absolute" &&  LL - LL.old < tol) || (crit == "relative" && (LL - LL.old)/abs(LL.old)  < tol) 
      if(converge) {
			  cat("converged at iteration",j,"with logLik:",LL,"\n")
        break
			}
		} else {
		  # this should not really happen...
		  if(j > 0 && (LL.old - LL) > tol) stop("likelihood decreased on iteration ",j)
		}
		
		LL.old <- LL
		#j <- j+1
		
		if(verbose&((j%%5)==0)) cat("iteration",j,"logLik:",LL,"\n")
	}
                      
  object@prior <- prior
  object@init <- init
  object@transition <- transition
  object@trDens <- trDens
  object@response <- response
  object@dens <- dens
		
	if(clsf == "hard") {
	    object <- as(object,"depmix.fitted.classLik") # class(object) <- "depmix.fitted.classLik"
	    object@posterior <- data.frame(state=viterbi(object)[,1])
	} else {
	    object <- as(object,"depmix.fitted") #  class(object) <- "depmix.fitted"
	    object@posterior <- viterbi(object)
	}
	
	if(converge) {
	    if(clsf == "hard") {
		    object@message <- switch(crit,
			    relative = "Log classification likelihood converged to within tol. (relative change)",
			    absolute = "Log classification likelihood converged to within tol. (absolute change)"
		    )	    
	    } else {
		    object@message <- switch(crit,
			    relative = "Log likelihood converged to within tol. (relative change)",
			    absolute = "Log likelihood converged to within tol. (absolute change)"
		    )
		}
	} else object@message <- "'maxit' iterations reached in EM without convergence."
	
	# no constraints in EM, except for the standard constraints ...
	# which are produced by the following (only necessary for getting df right in logLik and such)
	constraints <- getConstraints(object)
	object@conMat <- constraints$lincon
	object@lin.lower <- constraints$lin.l
	object@lin.upper <- constraints$lin.u
	
	object
}
