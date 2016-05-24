# hmm.R
# 
# S4 classes and methods for HMMs
#
# Author: Peter Humburg
###############################################################################

###############################################################################
## HMM generics	                                                             ##
###############################################################################

## A virtual HMM class
setClass("hmm",representation(transition.matrix="matrix",emission="list",init="numeric"),
	prototype(transition.matrix=matrix(),emission=list(),init=NULL))

## Define length of HMM as number of states
setMethod("length","hmm",function(x){dim(x@transition.matrix)[1]})

## Define new generics
## convenience functions to retrieve information about the model
setGeneric("states",def=function(hmm,...){standardGeneric("states")})

## Some default functions for common HMM methods
## Viterbi algorithm
setGeneric("viterbi",
	def=function(hmm,obs,...){standardGeneric("viterbi")},valueClass=c("list","character"))

## Forward algorithm
setGeneric("forward",
	def=function(hmm,obs,...){standardGeneric("forward")},valueClass=c("list","matrix","numeric"))

## Backward algorithm
setGeneric("backward",
	def=function(hmm,obs,...){standardGeneric("backward")},valueClass=c("list","matrix","numeric"))

## Baum-Welch algorithm
setGeneric("baumWelch",
	def=function(hmm,obs,...){standardGeneric("baumWelch")},valueClass=c("hmm"))

## One iteration of the Baum-Welch algorithm
setGeneric(".baumWelchStep", 
	def=function(hmm,obs,...){standardGeneric(".baumWelchStep")},valueClass=c("list"))

## parameter estimates for emission distribution. Select function by hmm and distribution class.
setGeneric(".baumWelchEmission",
	def=function(hmm,dist,obs,...) standardGeneric(".baumWelchEmission"),valueClass=c("list"))

## Viterbi training
setGeneric("viterbiTraining",
	def=function(hmm,obs,...){standardGeneric("viterbiTraining")}, valueClass=c("hmm"))

## emission probability estimates for Viterbi training
setGeneric(".viterbiTrainingEmission",
	def=function(hmm,obs.list,stateSeq,...){standardGeneric(".viterbiTrainingEmission")},valueClass=c("matrix","list"))

## Some other useful functions for HMMs

## Generate a sample from the HMM
setGeneric("sampleSeq",
	def=function(hmm,size,...){standardGeneric("sampleSeq")},valueClass=c("numeric","character","list"))

## Retreiving state names and emission alphabet
setMethod("states","hmm",
		function(hmm){
			rownames(hmm@transition.matrix)
		}		
)

#######################################################################
## HMM implementation for differnt types of emission distributions   ##
#######################################################################

## Class representing HMM with continuous observations
setClass("contHMM",representation(),prototype(),contains="hmm")
## Initialising contHMM
## transition and emission are lists of contDist objects
setMethod("initialize","contHMM",
			function(.Object,transition=list(),emission=list(),init=NULL){
				## Check list entries
				## ensure same number of states in transition and emission
				if(length(transition) != length(emission)) stop("Dimensions of 'transition' and 'emission' do not match!")
				## ensure all transition entries are probability distributions
				chk <- as.logical(lapply(transition,function(x) is(x,"discDist")))
				if(sum(chk) < length(chk)){
					chkLst <- transition[!chk]
					stop(paste("Expected object of class 'discDist' in argument 'transition', found ",class(chkLst[[1]])))
				}
				## ensure alphabet (set of states) of transition entries is making sense
				if(length(transition)) states <- transition[[1]]@alpha
				else states <- ""
				chk <- as.logical(unlist(lapply(transition,function(x) x@alpha == states)))
				if(sum(chk) < length(chk)) stop("Found non matching set of states!")
				
				## ensure all emission entries are probability distributions
				chk <- as.logical(unlist(lapply(emission,function(x) is(x,"contDist"))))
				if(sum(chk) < length(chk)){
					chkLst <- emission[!chk]
					stop(paste("Expected object of class 'contDist' in argument 'emission', found ",class(chkLst[[1]])))
				}
								
				## ensure initial distribution has the right state set
				if(length(init) == 0 || init@alpha != states){
					if(length(init) == 0 && length(states) != 0){
						prob <- numeric(length(states))+1/length(states)
						init <- new("discDist",alpha=states,prob=prob) 
					}
					else{if(length(states) > 0) stop("Illegal initial distribution!")}
				}
				
				## generate matrix of transition probabilities
				if(length(transition)){
					transition.matrix <- t(apply(cbind(lapply(transition,as.matrix)),1,unlist))
					if(dim(transition.matrix)[1] != dim(transition.matrix)[2])
						stop("Matrix of transition probabilities has to be square!")
					colnames(transition.matrix) <- states
					rownames(transition.matrix) <- states
				}
				.Object@transition.matrix <- transition.matrix
				.Object@emission <- emission
				.Object@init <- init@prob
				.Object
			}
)

######################################################################
## Access functions for HMMs                                        ##
######################################################################
"[.contHMM" <- function(x,i,j,transition=TRUE,log=FALSE,sum=TRUE,...){
	if(!missing(i) && !missing(j)){
		if(transition) value <- x@transition.matrix[i,j]
		else{
			value <- x@emission[[i]][,j,log=log]
		}
		if(sum && !transition) value <- sum(value)
		if(log && transition) value <- log(value)
	}
	if(missing(i) && !missing(j)){
		if(transition) value = x@transition.matrix[,j]
		else{
			value <- apply(as.matrix(1:length(x)),1,function(i,j,log,sum){x[i,j,FALSE,log,sum]},j,log,sum)
		}
		if(transition && log) value <- log(value)
	}
	if(!missing(i) && missing(j)){
		if(transition) value = x@transition.matrix[i,]
		else stop("No emission matrix available!")
		if(log) value <- log(value)
	}
	if(missing(i) && missing(j)){
		if(transition) value = x@transition.matrix
		else stop("No emission matrix available!")
		if(log) value <- log(value)
	}
	value
}

#############################################
## Plotting and printing functions         ##
#############################################
## Plotting emission distributions for contHMM objects
plot.contHMM <- function(x, ...){
	N <- length(x)
	names <- states(x)
	rows <- ceiling(sqrt(N))
	columns <- round(sqrt(N))
	xmin.index <- which.min(sapply(x@emission,function(em) min(em@components[,2])))
	xmax.index <- which.max(sapply(x@emission,function(em) max(em@components[,2])))
	xmin.index2 <- which.min(x@emission[[xmin.index]]@components[,2])
	xmax.index2 <- which.max(x@emission[[xmax.index]]@components[,2])
	xmin <- x@emission[[xmin.index]]@components[xmin.index2,2] - 4 * sqrt(x@emission[[xmin.index]]@components[xmin.index2,3])
	xmax <- x@emission[[xmax.index]]@components[xmax.index2,2] + 4 * sqrt(x@emission[[xmax.index]]@components[xmax.index2,3])
	
	par(mfrow=c(rows,columns))
	for(n in 1:N){
		plot(x@emission[[n]],main=names[n],xlim=c(xmin,xmax),...)
	}
}

## printing summary of conHMM objects
setMethod("show", "hmm",
	function(object){
		cat("An object of class \"", class(object), "\"\n", sep='')
		cat("with states:", names(object@emission), "\n")
		cat("\nInitial state distribution:\n")
		show(object@init)	
		cat("\nTransition matrix:\n")
		show(object@transition.matrix)
		cat("\nEmission distributions:\n")
		mapply(
			function(name, x) {
				cat("\"", name, "\":\n", sep='')
				show(x)
			},
			names(object@emission), object@emission
		)
		invisible(NULL)
	}
)


## Sample from continuous observation HMM
setMethod("sampleSeq",c("contHMM","numeric"),
		function(hmm,size,return.states=FALSE){
			if(size <= 0) return(character())
			stateSeq <- sample(states(hmm),1,prob=hmm@init)
			for(i in 2:size){
				stateSeq <- c(stateSeq,sample(states(hmm),1,prob=hmm[stateSeq[i-1],]))
			}
			obs <- numeric()
			for(s in stateSeq){
				obs <- c(obs,sampleObs(hmm@emission[[s]],1))
			}
			if(return.states){
				ret <- list(states=stateSeq,observation=obs)
			} else ret <- obs
			ret
		}
)

##################################################
## Implementation of Algorithms                 ##
##################################################

## Viterbi algorithm for HMMs with scalar observations
## hmm is an object of class hmm and obs is a vector of scalar observations
setMethod("viterbi",c("hmm","ANY"), function(hmm,obs=character(0),names=TRUE){
	if(!is(hmm@emission[[1]],"tDist")){
		N <- length(hmm)
		T <- length(obs)
		M <- matrix(0,nrow=N,ncol=T)
		
		A <- hmm[,,log=TRUE]
		## calculate probability of most likely state sequence, given obs
		M[,1] <- log(hmm@init) + hmm[,obs[1],F,log=TRUE]
		for(t in 2:T){
			for(j in 1:N){
				M[j,t] <- max(M[,t-1] + A[,j]) + hmm[j,obs[t],F,log=TRUE]
			}
		}
		prob <- max(M[,T])
		
		## calculate most likely state sequence given obs
		state.index <- numeric(T)
		state.index[T] <- which.max(M[,T])
		for(t in (T-1):1){
			state.index[t] <- which.max(M[,t] + A[,state.index[t+1]])
		}
		states <- state.index
	
		ret <- list()
		ret[["stateSeq"]] <- states
		ret[["logProb"]] <- prob
		ret[["matrix"]] <- M
	} else{
		ret <- .Call("_viterbi",hmm,obs)		
		if(names){ 
			ret$stateSeq <- states(hmm)[ret$stateSeq]
		}
	}
		ret
	}
)
## Forward algorithm for HMMs with scalar observations
setMethod("forward",c("hmm","ANY"),
		function(hmm,obs=character(0)){
			N <- length(hmm)
			T <- length(obs)
			
			if(is(hmm@emission[[1]],"tDist")){
			 	alpha <- .Call("_forward",hmm,obs,N,T)
			 } else{
				alpha <- matrix(-Inf,nrow=N,ncol=T) # scaled forward variables
				
				## initialisation
				alpha[,1] <- log(hmm@init) + hmm[,obs[1],FALSE,log=TRUE]
				A <- hmm[,,log=TRUE]
				## recursion
				for(t in 2:T){
					alpha[,t]<- t(t(apply(alpha[,t-1] + A,2,logSum)) + hmm[,obs[t],FALSE,log=TRUE])
				}
			}
			
			## log probability of obs given the parameters
			logProb <- logSum(alpha[,T])
			
			ret <- list()
			ret[["logProb"]] <- logProb
			ret[["alpha.scaled"]] <- alpha
			ret			
		}
)

## Backward algorithm for HMMs with scalar observations
setMethod("backward",c("hmm","ANY"),
		function(hmm,obs=character(0)){
			N <- length(hmm)
			T <- length(obs)
			
			if(is(hmm@emission[[1]],"tDist")){
				beta <- .Call("_backward",hmm,obs,N,T)
			} else{
				beta <- matrix(0,nrow=N,ncol=T) # scaled backward variables
				
				## initialisation
				A <- hmm[,,log=TRUE]
				## recursion
				for(t in (T-1):1){
					beta[,t] <- apply(t(t(A) + hmm[,obs[t+1],F,log=TRUE] + beta[,t+1]),1,logSum)
				}
			}	
			beta	
		}
)
## Functions for some parts of the Baum-Welch algorithm that don't change with type of HMM
## calculating forward and backward variables
.baumWelchInit <- function(hmm,obs){
	D <- length(obs)

	## calculating forward and backward variables
	f <- lapply(obs,function(obs,hmm){forward(hmm,obs)},hmm)
	alpha <- lapply(f,function(f){f$alpha.scaled})
	log.prob <- sum(sapply(f,function(f){f$logProb}))
	beta <- lapply(obs,function(obs,hmm){backward(hmm,obs)},hmm)
		
	ret <- list()
	ret[["alpha"]] <- alpha
	ret[["logProb"]] <- log.prob
	ret[["beta"]] <- beta
	ret		
}

## compute new transition probabilities
.baumWelchTransition <- function(hmm, obs, alpha, beta, trans.prior, init.prior){
	.Call("_baumWelch_trans", hmm, obs, alpha, beta, trans.prior, init.prior)
}

## calculate conditional expectation of U
.calcU <- function(density,obs){
	## get dimension of observation
	p <- 1
	if(is(obs,"matrix")) p <- dim(obs)[1]
	## calculate delta
	s <- density[,"variance"]
	delta <- ((obs - density[,"mean"])^2)/s
	df <- density[,"df"]
	
	(df + p)/(df + delta)
}
## function to find estimates for df
## the estimate is a root of this function
.df.fun <- function(df,p,n.tau.u){
	df1 <- 0.5 * df
	df.p <- (df + p) * 0.5
	- digamma(df1) + log(df1) + 1 + n.tau.u + digamma(df.p) - log(df.p) 
}

## Use root finding to determine degrees of freedom
.find.df <- function(tau,u,tau.sum=NULL){
	if(is.null(tau.sum)){
		tau.sum <- apply(sapply(tau,function(tau) apply(tau,1,sum)),1,sum)
	}
	tau.u <- mapply(function(tau,u) tau * (log(u) - u),tau,u,SIMPLIFY=FALSE)
	sum.tau.u <- sapply(tau.u,function(tau.u) apply(tau.u,1,sum))
	n.tau.u <- apply(sum.tau.u,1,sum) / tau.sum
	
	## find interval for root finding
	## we start with (0,20] and increase if necessary
	## (starting to close to 0 results in harmless but irritating warnings 
	## because .df.fun returns Inf)
	lower <- 10*.Machine$double.eps
	upper <- 20
	f.lower <- .df.fun(lower,1,n.tau.u)
	if(sum(is.na(f.lower))) stop("Cannot estimate df: NAs are not allowed")
	while(sum((f.lower * .df.fun(upper,1,n.tau.u)) > 0)){
		upper <- upper + 20
	}
	df <- mapply(function(n.tau.u) uniroot(.df.fun,lower=lower,upper=upper,p=1,n.tau.u=n.tau.u)$root, 
				n.tau.u)
	df		
}
## Method to estimate parameters of t distributions
## This implementation follows the EM-Algorithem for t mixtures as described 
## by Peel and McLachlan (2000)
## Degrees of freedom can be fixed by assigning a vector to df (one entry for each state)
setMethod(".baumWelchEmission",c("contHMM","tDist","list"),
		function(hmm,dist,obs,gamma,alpha,df=NULL){
			## transform back to linear space to handle negative observations
			tau <- lapply(gamma, exp)
			u <- lapply(obs,function(obs,hmm) t(sapply(hmm@emission,.calcU,obs)),hmm)
			tau.u <- mapply(function(tau,u) tau * u,tau,u, SIMPLIFY=FALSE)
			
			## get estimates for mu
			mean.num <- mapply(function(tau.u,obs)apply(t(t(tau.u) * obs),1,sum),tau.u,obs)
			mean.denom <- sapply(tau.u,function(tau.u)apply(tau.u,1,sum))
			## combine estimates from different observation sequences
			## use equal weights
			## TODO remove weights once we are sure they won't be used 
			weight <- 1 
	
			means <- numeric()
			for(i in 1:dim(mean.num)[1]){
				means <- c(means,sum(mean.num[i,]*weight)/sum(mean.denom[i,]*weight))
			}
			## get estimates for variance
			obs.diff <- lapply(obs,
				function(obs,means) t(apply(as.matrix(means),1,function(mu,obs)(obs - mu)^2, obs)),means)
			var.num <- mapply(function(tau.u,obs.diff)apply(tau.u * obs.diff,1,sum),tau.u,obs.diff)
			var.denom <- sapply(tau,function(tau) apply(tau,1,sum))
			## combine estimates from different observation sequences		
			var.denom.all <- apply(var.denom,1,function(vd,w)sum(vd*w),weight)
			var.num.all <- apply(var.num,1,function(vn,w)sum(vn*w),weight)
			vars <- var.num.all/var.denom.all

			## get estimates for df if required
			if(is.null(df)){
				## obtain single overall estimate for each state
				df <- .find.df(tau,u,var.denom.all)
			}
			emission <- mapply(function(mean,var,df) new("tDist",mean=mean,var=var,df=df),means,vars,df)
			
			emission			
		}
)
## Method for one iteration of Baum-Welch algorithm for HMMs with continuous observations
## returns a list of updated parameters
setMethod(".baumWelchStep",c("contHMM","list"),
		function(hmm,obs, trans.prior, init.prior, df=NULL){
			N <- length(hmm)
			D <- length(obs)
			T <- sapply(obs,length)
			
			init <- .baumWelchInit(hmm,obs)
			alpha <- init[["alpha"]]
			beta <- init[["beta"]]
			log.prob <- init[["logProb"]]
			## compute xi and gamma and estimate new transition probabilities
			trans <- .baumWelchTransition(hmm, obs, alpha, beta, trans.prior, init.prior)
			pi <- trans[["pi"]]
			transition <- trans[["transition"]]
			gamma <- trans[["gamma"]]
			rownames(transition) <- states(hmm)
			colnames(transition) <- states(hmm)			
			## estimate new emission distribution parameters
			emission <- .baumWelchEmission(hmm,hmm@emission[[1]],obs,gamma,alpha,df)
			names(emission) <- states(hmm)
			
			ret <- list()
			ret[["pi"]] <- pi
			ret[["transition"]] <- transition
			ret[["emission"]] <- emission
			ret[["prob"]] <- log.prob
			
			ret
	}
)
## Baum-Welch algorithm for HMMs with scalar emissions
setMethod("baumWelch",c("hmm","list"),
		function(hmm,obs,max.iter=FALSE,eps=0.01,df=NULL,trans.prior=NULL, init.prior=NULL,verbose=1){
			if(is.null(trans.prior)){
				## no prior
				trans.prior <- matrix(0,ncol=length(hmm),nrow=length(hmm))
			} else if(is.logical(trans.prior) && trans.prior){
				trans.prior <- hmm@transition.matrix
			}
			if(!is.matrix(trans.prior) || dim(trans.prior) != dim(hmm@transition.matrix)){
				stop("Illegal prior transition distribution")
			}
			## check prior initial state distribution
			if(is.null(init.prior)){
				## no prior
				init.prior <- numeric(length(hmm))
			} else if(is.logical(init.prior) && init.prior){
				init.prior <- hmm@init
			} 
			if(!is.vector(init.prior) || length(init.prior) != length(hmm)){
				stop("Illegal prior initial state distribution")
			}
			
			iter <- 1
			log.prob <- 1
			if(!is.null(df) && is(hmm@emission[[1]],"tDist")){
				for(i in 1:length(hmm)){
					hmm@emission[[i]][,"df"] <- df[((i-1)%%length(df))+1]
				}
			}
			repeat{
				params <- .baumWelchStep(hmm, obs, trans.prior, init.prior, df)
				if(log.prob == 1) {
					diff <- +Inf
				}else diff <- params$prob - log.prob
								
				if( log.prob != 1 & diff < 0){warning(paste("diff = ",as.character(params$prob),"-",as.character(log.prob),
						"=",as.character(params$prob - log.prob),"\niter = ",as.character(iter)))}
				hmm@init <- params$pi
				hmm@transition.matrix <- params$transition
				
				if(class(params$emission) == "matrix"){
					hmm@emission.matrix <- params$emission
				} else{
					hmm@emission <- params$emission
				}
				
				if(abs(diff) < eps || max.iter && iter >= max.iter) {
					if(verbose >= 1){ 
						message("Number of iterations: ",as.character(iter))
						message("Log likelihood of final model: ",params$prob)
						message("Last change in log likelihood: ",diff,"\n")
					}
					break
				}
				
				log.prob <- params$prob 
				if(verbose >= 2) message("iter: ",iter,"  llk: ",log.prob)
				iter <- iter + 1
			}
			
			hmm
		}
)
## Train an HMM using the Viterbi algorithm and successive maximum likelihood estimates
setMethod("viterbiTraining",c("hmm","list"),
		function(hmm,obs,max.iter=10,eps=0.01,df=NULL,trans.prior=NULL, init.prior=NULL,keep.models=NULL,verbose=1){
			## check transition prior
			if(is.null(trans.prior)){
				## no prior
				trans.prior <- matrix(0,ncol=length(hmm),nrow=length(hmm))
			} else if(is.logical(trans.prior) && trans.prior){
				trans.prior <- hmm@transition.matrix
			} 
			if(!is.matrix(trans.prior) || dim(trans.prior) != dim(hmm@transition.matrix)){
				stop("Illegal prior transition distribution")
			}
			## check prior initial state distribution
			if(is.null(init.prior)){
				## no prior
				init.prior <- numeric(length(hmm))
			} else if(is.logical(init.prior) && init.prior){
				init.prior <- hmm@init
			} 
			if(!is.vector(init.prior) || length(init.prior) != length(hmm)){
				stop("Illegal prior initial state distribution")
			}
			
			iter <- 1
			diff <- +Inf

			if(is.character(keep.models)){
				iter.models <- list()
				iter.llk <- list()
			}
			if(!is.null(df) && is(hmm@emission[[1]],"tDist")){
				for(i in 1:length(hmm)){
					hmm@emission[[i]]@components[,"df"] <- df[(i-1)%%length(df) +1]
				}
			}
			prob.new <- sum(sapply(obs,function(obs,hmm){f <- forward(hmm,obs);f$logProb},hmm))
			hmm.best <- hmm
			llk.best <- prob.new
			observed.best <- logical(length(hmm))
			while(iter <= max.iter && abs(diff) > eps){
				if(verbose >= 2) message(paste("iter:",iter,"  llk:",prob.new))
				## get state sequence
				seq <- lapply(obs,function(obs,hmm){vit <- viterbi(hmm,obs,FALSE);vit$stateSeq},hmm)
				## estimate for initial state distribution
				start <- sapply(seq,function(seq){seq[1]})
				start.count <- init.prior
				for(i in start){
					start.count[i] <- start.count[i] + 1
				}
				pi <- start.count/sum(start.count)
				## estimates for transition probabilities
				transitions <- lapply(seq,
						function(seq,len){
							counter <- matrix(0,nrow=len,ncol=len)
							for(t in 1:(length(seq)-1)){
								counter[seq[t],seq[t+1]] <- counter[seq[t],seq[t+1]]+1
							}
							counter 			
						},length(hmm))
				transitions.all <- trans.prior
				for(i in 1:length(transitions)){
					transitions.all <- transitions.all + transitions[[i]]
				}
				transitions.sum <- apply(transitions.all,1,sum)
				if(sum(transitions.sum == 0)){
					transitions.all[which(transitions.sum == 0), ] <- hmm@transition.matrix[which(transitions.sum == 0), ]
					if(sum(transitions.sum == 0) == 1){
						transitions.sum[which(transitions.sum == 0)] <- sum(trans.prior[which(transitions.sum == 0), ]) + 1
					} else{
						transitions.sum[which(transitions.sum == 0)] <- apply(as.matrix(trans.prior[which(transitions.sum == 0), ]), 1, sum) + 1
					}
				}
				A <- transitions.all/transitions.sum
				rownames(A) <- states(hmm)
				colnames(A) <- states(hmm)
				
				## estimates for emission probabilities
				B <- .viterbiTrainingEmission(hmm,obs,seq,df=df)
				## update HMM
				hmm@init <- pi
				hmm@transition.matrix <- A
				
				if(is(hmm,"discHMM")){
					hmm@emission.matrix <- B
				}
				if(is(hmm,"contHMM")){
					for(i in 1:length(hmm)){
						hmm@emission[[i]]@components <- B[[i]]
					}
				}
				prob.old <- prob.new
				prob.new <- sum(sapply(obs,function(obs,hmm){forward(hmm,obs)$logProb},hmm))
				diff <- prob.new - prob.old
				iter <- iter + 1
				if(is.character(keep.models)){
					iter.models[[iter]] <- hmm
					iter.llk[[iter]] <- prob.new
					save(iter.models,iter.llk,file=keep.models)
				}
				if(prob.new > llk.best){
					hmm.best <- hmm
					llk.best <- prob.new
					observed.best <- B$observed
				}
			}
			if(sum(!observed.best) && sum(observed.best)){
				warning("Data contains insufficient observations for the following states:\n\t", 
						states(hmm)[!observed.best],
						"\nConsider reducing the number of states in the model or increasing prior probabilities ",
						"for transitions to these states.", call.=FALSE)
			}
			if(verbose >= 1){ 
						message("Number of iterations: ",as.character(iter-1))
						message("Log likelihood of best model: ",llk.best)
						message("Last change in log likelihood: ",diff,"\n")
					}
			hmm.best
		}
)

## TODO: Remove code not relevant for t distributions
## emission probability estimates for Gaussian mixture distributions
setMethod(".viterbiTrainingEmission",c("contHMM","list","list"),
		function(hmm,obs.list,stateSeq,multi.cycle=FALSE,df=NULL,...){
			obs <- c(obs.list,recursive=TRUE)
			stateSeq <- c(stateSeq,recursive=TRUE)
			index <- matrix(FALSE,nrow=length(hmm),ncol=length(obs))
			for(i in 1:length(hmm)){
				index[i,] <- stateSeq == i
			}
			## estimates for mixture components
			components <- list()
			observed <- logical(length(hmm))
			if(is(hmm@emission[[1]],"tDist")){
				for(i in 1:length(hmm)){
					if(sum(index[i,]) >= 2){
						k <- dim(hmm@emission[[i]]@components)[1]
						components[[i]] <- matrix(0,ncol=4,nrow=k)
						## estimate location parameter
						mu <- mean(obs[index[i,]])
						variance <- var(obs[index[i,]])
						df.est <- hmm@emission[[i]]@components[1,4]
						
						components[[i]][1,] <- c(1,mu,variance,df.est)
						colnames(components[[i]]) <- c("weight","mean","variance","df")
						observed[i] <- TRUE
					} else{
						## insufficien observations for state i, re-use previous estimates
						components[[i]] <- hmm@emission[[i]]@components
					}
				}
				## Update state sequence estimate for multi cycle estimation
				if(multi.cycle){
					for(i in 1:length(hmm)){
						hmm@emission[[i]]@components <- components[[i]]
					}
					stateSeq <- c(lapply(obs.list,function(obs,hmm)viterbi(hmm,obs)$stateSeq,hmm),recursive=TRUE)
				}
				if(is.null(df)){
					for(i in 1:length(hmm)){
						if(observed[i]){
							## estimate degrees of freedom
							U <- .calcU(hmm@emission[[i]],obs)
							idx.u <- log(U[index[i,]]) - U[index[i,]]
							idx.u.sum <- sum(idx.u)
							n.idx.u <- idx.u.sum / sum(index[i,])
							## find interval for root finding
							## we start with (0,20] and increase if necessary
							## (starting to close to 0 results in harmless but irritating warnings 
							## because .df.fun returns Inf)
							lower <- 10*.Machine$double.eps
							upper <- 20
							f.lower <- .df.fun(lower,1,n.idx.u)
							
							while(f.lower * .df.fun(upper,1,n.idx.u) > 0) upper <- upper + 20
							df <- uniroot(.df.fun,lower=lower,upper=upper,p=1,n.tau.u=n.idx.u)$root
							components[[i]][1,4] <- df
						}
					}
				}
			} else{		
				for(i in 1:length(hmm)){
					k <- dim(hmm@emission[[i]]@components)[1]
					## simple mixture
					if(k == 1){
						components[[i]] <- matrix(0,ncol=3,nrow=k)
						components[[i]][1,] <- c(1,mean(obs[index[i,]]),var(obs[index[i,]]))
						colnames(components[[i]]) <- c("weight","mean","var")
					} else{ ## several mixture components
						components[[i]] <- matrix(0,ncol=3,nrow=k)
						## assign weighted observations to mixture components
						## the probability of observations under the current components is used to derive weights
						prob <- apply(as.matrix(obs[index[i,]]),1,function(obs,i) hmm[i,obs,transition=FALSE,sum=FALSE],i)
						obs.weight <- apply(prob,2,function(x)x/sum(x))
						weights <- apply(obs.weight,1,sum)
						weights <- weights/sum(weights)
						means <- apply(as.matrix(1:dim(obs.weight)[1]),1,
								function(i,x,w) x[i,]/(length(x[i,])*w[i]),t(t(obs.weight) * obs[index[i,]]),weights)
						means <- apply(means,2,sum)
						vars <- apply(as.matrix(means),1,function(m,o) (o - m)*(o - m) ,obs[index[i,]])
						vars <- vars * t(obs.weight)
						vars <- apply(vars,2,sum)/(weights*length(obs[index[i,]])-1)
						
						components[[i]] <- cbind(weights,means,vars)
						colnames(components[[i]]) <- c("weight","mean","var")
					}
									
				}
			}
			components[["observed"]] <- observed
			components
		}
)
                             
