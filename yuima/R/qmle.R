##::quasi-likelihood function

##::extract drift term from yuima
##::para: parameter of drift term (theta2)

### TO BE FIXED: all caculations should be made on a private environment to
### avoid problems.
### I have rewritten drift.term and diff.term instead of calc.drift and
### calc.diffusion to make them independent of the specification of the
### parameters.  S.M.I. 22/06/2010

drift.term <- function(yuima, theta, env){
	r.size <- yuima@model@noise.number
	d.size <- yuima@model@equation.number
	modelstate <- yuima@model@state.variable
	DRIFT <- yuima@model@drift
#	n <- length(yuima)[1]
	n <- dim(env$X)[1]

	drift <- matrix(0,n,d.size)
	tmp.env <- new.env()
	assign(yuima@model@time.variable, env$time, envir=tmp.env)


	for(i in 1:length(theta)){
		assign(names(theta)[i],theta[[i]], envir=tmp.env)
	}

	for(d in 1:d.size){
		assign(modelstate[d], env$X[,d], envir=tmp.env)
	}
	for(d in 1:d.size){
		drift[,d] <- eval(DRIFT[d], envir=tmp.env)
	}

	return(drift)
}


diffusion.term <- function(yuima, theta, env){
	r.size <- yuima@model@noise.number
	d.size <- yuima@model@equation.number
	modelstate <- yuima@model@state.variable
	DIFFUSION <- yuima@model@diffusion
#	n <- length(yuima)[1]
	n <- dim(env$X)[1]
    tmp.env <- new.env()
	assign(yuima@model@time.variable, env$time, envir=tmp.env)
	diff <- array(0, dim=c(d.size, r.size, n))
	for(i in 1:length(theta)){
		assign(names(theta)[i],theta[[i]],envir=tmp.env)
	}

	for(d in 1:d.size){
		assign(modelstate[d], env$X[,d], envir=tmp.env)
	}

	for(r in 1:r.size){
		for(d in 1:d.size){
			diff[d, r, ] <- eval(DIFFUSION[[d]][r], envir=tmp.env)
		}
	}
	return(diff)
}


## Koike's code
##::extract jump term from yuima
##::gamma: parameter of diffusion term (theta3)
measure.term <- function(yuima, theta, env){
    r.size <- yuima@model@noise.number
    d.size <- yuima@model@equation.number
    modelstate <- yuima@model@state.variable
    n <- dim(env$X)[1]

    tmp.env <- new.env()
    assign(yuima@model@time.variable, env$time, envir =tmp.env)
    JUMP <- yuima@model@jump.coeff
    measure <- array(0, dim=c(d.size, r.size, n))
    for(i in 1:length(theta)){
        assign(names(theta)[i],theta[[i]],envir=tmp.env)
    }

    for(d in 1:d.size){
        assign(modelstate[d], env$X[,d],envir=tmp.env)
    }
    for(r in 1:r.size){
        #for(d.tmp in 1:d){
        if(d.size==1){
            measure[1,r,] <- eval(JUMP[[r]],envir=tmp.env)
        }else{
            for(d in 1:d.size){
                measure[d,r,] <- eval(JUMP[[d]][r],envir=tmp.env)
            }
        }
    }
    return(measure)
}


### I have rewritten qmle as a version of ml.ql
### This function has an interface more similar to mle.
### ml.ql is limited in that it uses fixed names for drift and diffusion
### parameters, while yuima model allows for any names.
### also, I am using the same interface of optim to specify upper and lower bounds
### S.M.I. 22/06/2010

is.Poisson <- function(obj){
    if(is(obj,"yuima"))
    return(is(obj@model, "yuima.poisson"))
    if(is(obj,"yuima.model"))
    return(is(obj, "yuima.poisson"))
    return(FALSE)
}

is.CARMA <- function(obj){
 if(is(obj,"yuima"))
    return(is(obj@model, "yuima.carma"))
 if(is(obj,"yuima.model"))
    return(is(obj, "yuima.carma"))
 return(FALSE)
}

qmle <- function(yuima, start, method="BFGS", fixed = list(), print=FALSE,
 lower, upper, joint=FALSE, Est.Incr="Carma.IncPar",aggregation=TRUE, threshold=NULL, ...){
  if(is(yuima@model, "yuima.carma")){
    NoNeg.Noise<-FALSE
    cat("\nStarting qmle for carma ... \n")
  }
  if(is.CARMA(yuima)&& length(yuima@model@info@scale.par)!=0){
    method<-"L-BFGS-B"
  }
	call <- match.call()

	if( missing(yuima))
		yuima.stop("yuima object is missing.")
	if(is.COGARCH(yuima)){
	  if(missing(lower))
	    lower <- list()

	  if(missing(upper))
	    upper <- list()

	 res <- NULL
	 if("grideq" %in% names(as.list(call)[-(1:2)])){
	 res  <- PseudoLogLik.COGARCH(yuima, start, method=method, fixed = list(),
	                       lower, upper, Est.Incr, call, ...)
	 }else{
	   res  <- PseudoLogLik.COGARCH(yuima, start, method=method, fixed = list(),
	                         lower, upper, Est.Incr, call, grideq = FALSE,...)
	 }

	 return(res)
	}

    orig.fixed <- fixed
    orig.fixed.par <- names(orig.fixed)
    if(is.Poisson(yuima))
     threshold <- 0
## param handling

## FIXME: maybe we should choose initial values at random within lower/upper
##        at present, qmle stops
	if( missing(start) )
	 yuima.stop("Starting values for the parameters are missing.")

  #14/12/2013 We modify the QMLE function when the model is a Carma(p,q).
  # In this case we use a two step procedure:
  # First) The Coefficient are obtained by QMLE computed using the Kalman Filter.
  # Second) Using the result in Brockwell, Davis and Yang (2007) we retrieve
  # the underlying Levy. The estimated increments are used to find the L?vy parameters.

#   if(is(yuima@model, "yuima.carma")){
#     yuima.warm("two step procedure for carma(p,q)")
#     return(null)
#   }
#

    yuima.nobs <- as.integer(max(unlist(lapply(get.zoo.data(yuima),length))-1,na.rm=TRUE))

    diff.par <- yuima@model@parameter@diffusion

#	24/12
  if(is.CARMA(yuima) && length(diff.par)==0
	   && length(yuima@model@parameter@jump)!=0){
    diff.par<-yuima@model@parameter@jump
	}

  if(is.CARMA(yuima) && length(yuima@model@parameter@jump)!=0){
    CPlist <- c("dgamma", "dexp")
    codelist <- c("rIG", "rgamma")
    if(yuima@model@measure.type=="CP"){
      tmp <- regexpr("\\(", yuima@model@measure$df$exp)[1]
      measurefunc <- substring(yuima@model@measure$df$exp, 1, tmp-1)
      if(!is.na(match(measurefunc,CPlist))){
        yuima.warn("carma(p,q): the qmle for a carma(p,q) driven by a Compound Poisson with no-negative random size")
        NoNeg.Noise<-TRUE
        # we need to add a new parameter for the mean of the Noise
        if((yuima@model@info@q+1)==(yuima@model@info@q+1)){
          start[["mean.noise"]]<-1
        }
  #      return(NULL)
      }

    }

    if(yuima@model@measure.type=="code"){
      tmp <- regexpr("\\(", yuima@model@measure$df$exp)[1]
      measurefunc <- substring(yuima@model@measure$df$exp, 1, tmp-1)
      if(!is.na(match(measurefunc,codelist))){
        yuima.warn("carma(p,q): the qmle for a carma(p,q) driven by a non-Negative Levy  will be implemented as soon as possible")
        NoNeg.Noise<-TRUE
        if((yuima@model@info@q+1)==(yuima@model@info@q+1)){
          start[["mean.noise"]]<-1
        }
        #return(NULL)
      }
    }


#     yuima.warn("carma(p,q): the qmle for a carma(p,q) driven by a Jump process will be implemented as soon as possible ")
#     return(NULL)
  }

  # 24/12
  if(is.CARMA(yuima) && length(yuima@model@info@lin.par)>0){
    yuima.warn("carma(p,q): the case of lin.par will be implemented as soon as")
    return(NULL)
  }


  drift.par <- yuima@model@parameter@drift
#01/01 we introduce the new variable in order
# to take into account the parameters in the starting conditions

  if(is.CARMA(yuima)){
    #if(length(yuima@model@info@scale.par)!=0){
      xinit.par <- yuima@model@parameter@xinit
    #}
  }

  # SMI-2/9/14: measure.par is used for Compound Poisson
  # and CARMA, jump.par only by CARMA
  jump.par <- NULL
  if(is.CARMA(yuima)){
      jump.par <- yuima@model@parameter@jump
      measure.par <- yuima@model@parameter@measure
  } else {
      if(length(yuima@model@parameter@jump)!=0){
          measure.par <- yuima@model@parameter@jump
      } else {
      measure.par <- yuima@model@parameter@measure
      }
  }
  # jump.par is used for CARMA
  common.par <- yuima@model@parameter@common

  JointOptim <- joint
  if(is.CARMA(yuima) && length(yuima@model@parameter@jump)!=0){
    if(any((match(jump.par, drift.par)))){
      JointOptim <- TRUE
      yuima.warn("Drift and diffusion parameters must be different. Doing
					  joint estimation, asymptotic theory may not hold true.")
    }
  }

	if(length(common.par)>0){
		JointOptim <- TRUE
		yuima.warn("Drift and diffusion parameters must be different. Doing
					  joint estimation, asymptotic theory may not hold true.")
	# 24/12
#     if(is(yuima@model, "yuima.carma")){
#            JointOptim <- TRUE
# 		       yuima.warm("Carma(p.q): The case of common parameters in Drift and Diffusion Term will be implemented as soon as possible,")
# 		       #return(NULL)
# 		     }
	}

# if(!is(yuima@model, "yuima.carma")){
#    	if(length(jump.par)+length(measure.par)>0)
#    		yuima.stop("Cannot estimate the jump models, yet")
#	 }


	if(!is.list(start))
		yuima.stop("Argument 'start' must be of list type.")

	fullcoef <- NULL

	if(length(diff.par)>0)
	 fullcoef <- diff.par

	if(length(drift.par)>0)
	 fullcoef <- c(fullcoef, drift.par)

  if(is.CARMA(yuima) &&
       (length(yuima@model@info@loc.par)!=0)){
    # 01/01 We modify the code for considering
    # the loc.par in yuima.carma model
    fullcoef<-c(fullcoef, yuima@model@info@loc.par)

  }

  if(is.CARMA(yuima) && (NoNeg.Noise==TRUE)){
    if((yuima@model@info@q+1)==yuima@model@info@p){
      mean.noise<-"mean.noise"
      fullcoef<-c(fullcoef, mean.noise)
    }
  }

  #  if(is.CARMA(yuima) && (length(measure.par)>0)){
   fullcoef<-c(fullcoef, measure.par)
  #}

	npar <- length(fullcoef)


	fixed.par <- names(fixed) # We use Fixed.par when we consider a Carma with scale parameter
  if(is.CARMA(yuima) && (length(measure.par)>0)){
    fixed.carma=NULL
    if(!missing(fixed)){
      if(names(fixed) %in% measure.par){
        idx.fixed.carma<-match(names(fixed),measure.par)
        idx.fixed.carma<-idx.fixed.carma[!is.na(idx.fixed.carma)]
        if(length(idx.fixed.carma)!=0){
          fixed.carma<-as.numeric(fixed[measure.par[idx.fixed.carma]])
          names(fixed.carma)<-measure.par[idx.fixed.carma]
        }
      }
    }
    upper.carma=NULL
    if(!missing(upper)){
      if(names(upper) %in% measure.par){
        idx.upper.carma<-match(names(upper),measure.par)
        idx.upper.carma<-idx.upper.carma[!is.na(idx.upper.carma)]
        if(length(idx.upper.carma)!=0){
          upper.carma<-as.numeric(upper[measure.par[idx.upper.carma]])
          names(upper.carma)<-measure.par[idx.upper.carma]
        }
      }
    }
    lower.carma=NULL
    if(!missing(lower)){
      if(names(lower) %in% measure.par){
        idx.lower.carma<-match(names(lower),measure.par)
        idx.lower.carma<-idx.lower.carma[!is.na(idx.lower.carma)]
        if(length(idx.lower.carma)!=0){
          lower.carma<-as.numeric(lower[measure.par[idx.lower.carma]])
          names(lower.carma)<-measure.par[idx.lower.carma]
        }
      }
    }




    for( j in c(1:length(measure.par))){
          if(is.na(match(measure.par[j],names(fixed)))){
          fixed.par <- c(fixed.par,measure.par[j])
          fixed[measure.par[j]]<-start[measure.par[j]]
      }
    }

  }
	if (any(!(fixed.par %in% fullcoef)))
	 yuima.stop("Some named arguments in 'fixed' are not arguments to the supplied yuima model")

    nm <- names(start)

    oo <- match(nm, fullcoef)

    if(any(is.na(oo)))
		yuima.stop("some named arguments in 'start' are not arguments to the supplied yuima model")
    start <- start[order(oo)]
    nm <- names(start)

	idx.diff <- match(diff.par, nm)
	idx.drift <- match(drift.par, nm)
    # SMI-2/9/14: idx.measure for CP
    idx.measure <- match(measure.par, nm)
	#01/01
  if(is.CARMA(yuima)){
   # if(length(yuima@model@info@scale.par)!=0){
      idx.xinit <- as.integer(na.omit(match(xinit.par,nm)))# We need to add idx if NoNeg.Noise is TRUE
    #}
  }

  idx.fixed <- match(fixed.par, nm)
  orig.idx.fixed <- idx.fixed

	tmplower <- as.list( rep( -Inf, length(nm)))
	names(tmplower) <- nm
	if(!missing(lower)){
	   idx <- match(names(lower), names(tmplower))
	   if(any(is.na(idx)))
		yuima.stop("names in 'lower' do not match names fo parameters")
	   tmplower[ idx ] <- lower
	}
	lower <- tmplower

	tmpupper <- as.list( rep( Inf, length(nm)))
	names(tmpupper) <- nm
	if(!missing(upper)){
		idx <- match(names(upper), names(tmpupper))
		if(any(is.na(idx)))
			yuima.stop("names in 'lower' do not match names fo parameters")
		tmpupper[ idx ] <- upper
	}
	upper <- tmpupper





	d.size <- yuima@model@equation.number
	if (is.CARMA(yuima)){
	  # 24/12
    d.size <-1
	}
	n <- length(yuima)[1]

	env <- new.env()

    assign("X",  as.matrix(onezoo(yuima)), envir=env)
  	assign("deltaX",  matrix(0, n-1, d.size), envir=env)
    # SMI-2/9/14: for CP
    assign("Cn.r", numeric(n-1), envir=env)
    if(length(measure.par)==0)
        threshold <- 0  # there are no jumps, we take all observations

  if (is.CARMA(yuima)){
    #24/12 If we consider a carma model,
    # the observations are only the first column of env$X
#     assign("X",  as.matrix(onezoo(yuima)), envir=env)
#     env$X<-as.matrix(env$X[,1])
#     assign("deltaX",  matrix(0, n-1, d.size)[,1], envir=env)
     	  env$X<-as.matrix(env$X[,1])
#     	  env$X<-na.omit(as.matrix(env$X[,1]))
     	  env$deltaX<-as.matrix(env$deltaX[,1])
    assign("time.obs",length(env$X),envir=env)
#   env$time.obs<-length(env$X)
    #p <-yuima@model@info@p
    assign("p", yuima@model@info@p, envir=env)
    assign("q", yuima@model@info@q, envir=env)
    assign("V_inf0", matrix(diag(rep(1,env$p)),env$p,env$p), envir=env)


#     env$X<-as.matrix(env$X[,1])
# 	  env$deltaX<-as.matrix(env$deltaX[,1])
#     assign("time.obs",length(env$X), envir=env)
# 	  p <-yuima@model@info@p
# 	  assign("V_inf0", matrix(diag(rep(1,p)),p,p), envir=env)
	}
  assign("time", as.numeric(index(yuima@data@zoo.data[[1]])), envir=env)

    for(t in 1:(n-1)){
        env$deltaX[t,] <- env$X[t+1,] - env$X[t,]
        if(!is.CARMA(yuima))
            env$Cn.r[t] <- ((sqrt( env$deltaX[t,] %*% env$deltaX[t,])) <= threshold)
    }

    if(length(measure.par)==0)
        env$Cn.r <- rep(1, length(env$Cn.r))  # there are no jumps, we take all observations

	assign("h", deltat(yuima@data@zoo.data[[1]]), envir=env)

#SMI: 2/9/214 jump
if(length(measure.par)>0){


    args <- unlist(strsplit(suppressWarnings(sub("^.+?\\((.+)\\)", "\\1",yuima@model@measure$df$expr,perl=TRUE)), ","))
    idx.intensity <- numeric(0)
    for(i in 1:length(measure.par)){
        if(sum(grepl(measure.par[i],yuima@model@measure$intensity)))
        idx.intensity <- append(idx.intensity,i)
    }

    assign("idx.intensity", idx.intensity, envir=env)
    assign("measure.var", args[1], envir=env)
}


	f <- function(p) {
        mycoef <- as.list(p)
        if(!is.CARMA(yuima)){
            if(length(c(idx.fixed,idx.measure))>0) ## SMI 2/9/14
                names(mycoef) <- nm[-c(idx.fixed,idx.measure)] ## SMI 2/9/14
                else
                    names(mycoef) <- nm
        } else {
            if(length(idx.fixed)>0)
            names(mycoef) <- nm[-idx.fixed]
            else
            names(mycoef) <- nm
        }
        mycoef[fixed.par] <- fixed
	    minusquasilogl(yuima=yuima, param=mycoef, print=print, env)
    }

# SMI-2/9/14:
        fpsi <- function(p){
            mycoef <- as.list(p)

            idx.cont <- c(idx.diff,idx.drift)
            if(length(c(idx.fixed,idx.cont))>0)
            names(mycoef) <- nm[-c(idx.fixed,idx.cont)]
            else
            names(mycoef) <- nm
            mycoef[fixed.par] <- fixed
            #            print(mycoef)
            #print(p)
            minusquasipsi(yuima=yuima, param=mycoef, print=print, env=env)
        }


	 fj <- function(p) {
		 mycoef <- as.list(p)
         #		 names(mycoef) <- nm
         if(!is.CARMA(yuima)){
          idx.fixed <- orig.idx.fixed
          if(length(c(idx.fixed,idx.measure))>0) ## SMI 2/9/14
           names(mycoef) <- nm[-c(idx.fixed,idx.measure)] ## SMI 2/9/14
          else
		  names(mycoef) <- nm
         } else {
             names(mycoef) <- nm
             mycoef[fixed.par] <- fixed
	     }
         mycoef[fixed.par] <- fixed
		 minusquasilogl(yuima=yuima, param=mycoef, print=print, env)
	 }

	 oout <- NULL
     HESS <- matrix(0, length(nm), length(nm))
	 colnames(HESS) <- nm
	 rownames(HESS) <- nm


     HaveDriftHess <- FALSE
	 HaveDiffHess <- FALSE
	 HaveMeasHess <- FALSE


    if(length(start)){
		if(JointOptim){ ### joint optimization
            old.fixed <- fixed
            new.start <- start
            old.start <- start
            if(!is.CARMA(yuima)){
             if(length(c(idx.fixed,idx.measure))>0)
              new.start <- start[-c(idx.fixed,idx.measure)] # considering only initial guess for
            }

            if(length(new.start)>1){ #??multidimensional optim # Adjust lower for no negative Noise
                if(is.CARMA(yuima) && (NoNeg.Noise==TRUE))
                    if(mean.noise %in% names(lower)){lower[mean.noise]<-10^-7}
				oout <- optim(new.start, fj, method = method, hessian = TRUE, lower=lower, upper=upper)

            if(is.CARMA(yuima)){
                HESS <- oout$hessian
            } else {
                HESS[names(new.start),names(new.start)] <- oout$hessian
            }


				if(is.CARMA(yuima) && length(yuima@model@info@scale.par)!=0){
				  b0<-paste0(yuima@model@info@ma.par,"0",collapse="")
				  idx.b0<-match(b0,rownames(HESS))
				  HESS<-HESS[-idx.b0,]
				  HESS<-HESS[,-idx.b0]
				}
				if(is.CARMA(yuima) && length(yuima@model@parameter@measure)!=0){
				  for(i in c(1:length(fixed.par))){
                      indx.fixed<-match(fixed.par[i],rownames(HESS))
                      HESS<-HESS[-indx.fixed,]
                      HESS<-HESS[,-indx.fixed]
				  }
                  if(is.CARMA(yuima) && (NoNeg.Noise==TRUE)){
                      idx.noise<-(match(mean.noise,rownames(HESS)))
                      HESS<-HESS[-idx.noise,]
                      HESS<-HESS[,-idx.noise]
                  }
				}
				HaveDriftHess <- TRUE
				HaveDiffHess <- TRUE
			} else { ### one dimensional optim
				opt1 <- optimize(f, ...) ## an interval should be provided
                oout <- list(par = opt1$minimum, value = opt1$objective)
			} ### endif( length(start)>1 )
          theta1 <- oout$par[diff.par]
          theta2 <- oout$par[drift.par]

		} else {  ### first diffusion, then drift
			theta1 <- NULL

			old.fixed <- fixed
			old.start <- start

			if(length(idx.diff)>0){
## DIFFUSION ESTIMATIOn first
			old.fixed <- fixed
			old.start <- start
            old.fixed.par <- fixed.par
			new.start <- start[idx.diff] # considering only initial guess for diffusion
			new.fixed <- fixed
			if(length(idx.drift)>0)
			 new.fixed[nm[idx.drift]] <- start[idx.drift]
			fixed <- new.fixed
			fixed.par <- names(fixed)
			idx.fixed <- match(fixed.par, nm)
			names(new.start) <- nm[idx.diff]

			mydots <- as.list(call)[-(1:2)]
			mydots$print <- NULL
			mydots$fixed <- NULL
			mydots$fn <- as.name("f")
			mydots$start <- NULL
			mydots$par <- unlist(new.start)
			mydots$hessian <- FALSE
			mydots$upper <- as.numeric(unlist( upper[ nm[idx.diff] ]))
			mydots$lower <- as.numeric(unlist( lower[ nm[idx.diff] ]))
            mydots$threshold <- NULL #SMI 2/9/14

           if((length(mydots$par)>1) | any(is.infinite(c(mydots$upper,mydots$lower)))){
            oout <- do.call(optim, args=mydots)
           } else {
			 mydots$f <- mydots$fn
			 mydots$fn <- NULL
			 mydots$par <- NULL
			 mydots$hessian <- NULL
			 mydots$method <- NULL
			 mydots$interval <- as.numeric(c(unlist(lower[diff.par]),unlist(upper[diff.par])))


             mydots$lower <- NULL
			 mydots$upper <- NULL
			 opt1 <- do.call(optimize, args=mydots)
			 theta1 <- opt1$minimum
			 names(theta1) <- diff.par
			 oout <- list(par = theta1, value = opt1$objective)
			}
			theta1 <- oout$par

            fixed <- old.fixed
			start <- old.start
            fixed.par <- old.fixed.par

			} ## endif(length(idx.diff)>0)

			theta2 <- NULL


			if(length(idx.drift)>0){
## DRIFT estimation with first state diffusion estimates
			fixed <- old.fixed
			start <- old.start
            old.fixed.par <- fixed.par
			new.start <- start[idx.drift] # considering only initial guess for drift
			new.fixed <- fixed
			new.fixed[names(theta1)] <- theta1
			fixed <- new.fixed
			fixed.par <- names(fixed)
			idx.fixed <- match(fixed.par, nm)

            names(new.start) <- nm[idx.drift]

			mydots <- as.list(call)[-(1:2)]
			mydots$print <- NULL
			mydots$fixed <- NULL
			mydots$fn <- as.name("f")
            mydots$threshold <- NULL #SMI 2/9/14

			mydots$start <- NULL
			mydots$par <- unlist(new.start)
			mydots$hessian <- FALSE
			mydots$upper <- unlist( upper[ nm[idx.drift] ])
			mydots$lower <- unlist( lower[ nm[idx.drift] ])





			if(length(mydots$par)>1 | any(is.infinite(c(mydots$upper,mydots$lower)))){
			  if(is.CARMA(yuima)){
			    if(NoNeg.Noise==TRUE){
                    if((yuima@model@info@q+1)==yuima@model@info@p){
                        mydots$lower[names(start["NoNeg.Noise"])]<-10^(-7)
                    }

			    }
                if(length(yuima@model@info@scale.par)!=0){
                    name_b0<-paste0(yuima@model@info@ma.par,"0",collapse="")
                    index_b0<-match(name_b0,nm)
                    mydots$lower[index_b0]<-1
                    mydots$upper[index_b0]<-1+10^(-7)
			    }
                if (length(yuima@model@info@loc.par)!=0){
                    mydots$upper <- unlist( upper[ nm ])
                    mydots$lower <- unlist( lower[ nm ])
                    idx.tot<-unique(c(idx.drift,idx.xinit))
                    new.start <- start[idx.tot]
                    names(new.start) <- nm[idx.tot]
                    mydots$par <- unlist(new.start)
                }
			  }  # END if(is.CARMA)



            oout1 <- do.call(optim, args=mydots)


	#		  oout1 <- optim(mydots$par,f,method = "L-BFGS-B" , lower = mydots$lower, upper = mydots$upper)
			} else {
				mydots$f <- mydots$fn
				mydots$fn <- NULL
				mydots$par <- NULL
				mydots$hessian <- NULL
				mydots$method <- NULL
				mydots$interval <- as.numeric(c(lower[drift.par],upper[drift.par]))
				opt1 <- do.call(optimize, args=mydots)
				theta2 <- opt1$minimum
				names(theta2) <- drift.par
				oout1 <- list(par = theta2, value = as.numeric(opt1$objective))
			}
			theta2 <- oout1$par
            fixed <- old.fixed
			start <- old.start
            old.fixed.par <- fixed.par
			} ## endif(length(idx.drift)>0)


			oout1 <- list(par=  c(theta1, theta2))
      if (! is.CARMA(yuima)){
              if(length(c(diff.par, diff.par))>0)
                names(oout1$par) <- c(diff.par,drift.par)
      }


			oout <- oout1

		} ### endif JointOptim
    } else {
		list(par = numeric(0L), value = f(start))
	}


      fMeas <- function(p) {
            mycoef <- as.list(p)
            #  if(! is.CARMA(yuima)){
            #    #                names(mycoef) <- drift.par
            #    mycoef[measure.par] <- coef[measure.par]
            #}
            minusquasipsi(yuima=yuima, param=mycoef, print=print, env=env)
            #            minusquasilogl(yuima=yuima, param=mycoef, print=print, env)
        }


    fDrift <- function(p) {
		 mycoef <- as.list(p)
     if(! is.CARMA(yuima)){
       names(mycoef) <- drift.par
  		 mycoef[diff.par] <- coef[diff.par]
     }
		 minusquasilogl(yuima=yuima, param=mycoef, print=print, env)
	 }

	 fDiff <- function(p) {
	   mycoef <- as.list(p)
     if(! is.CARMA(yuima)){
  		 names(mycoef) <- diff.par
  		 mycoef[drift.par] <- coef[drift.par]
     }
		 minusquasilogl(yuima=yuima, param=mycoef, print=print, env)
	 }

     # coef <- oout$par
	 #control=list()
	 #par <- coef

#names(par) <- unique(c(diff.par, drift.par))
#     nm <- unique(c(diff.par, drift.par))

# START: ESTIMATION OF CP part
      theta3 <- NULL

       if(length(idx.measure)>0 & !is.CARMA(yuima)){
           idx.cont <- c(idx.drift,idx.diff)

           fixed <- old.fixed
           start <- old.start
           old.fixed.par <- fixed.par
           new.fixed <- fixed

           new.start <- start[idx.measure] # considering only initial guess for measure
           new.fixed <- fixed

           new.fixed[names(theta1)] <- theta1
           new.fixed[names(theta2)] <- theta2

           fixed <- new.fixed
           fixed.par <- names(fixed)
           idx.fixed <- match(fixed.par, nm)
           #            names(new.start) <- nm[idx.drift]
           names(new.start) <- nm[idx.measure]

           mydots <- as.list(call)[-(1:2)]
           #    mydots$print <- NULL
           mydots$threshold <- NULL
           mydots$fixed <- NULL
           mydots$fn <- as.name("fpsi")
           mydots$start <- NULL
           mydots$threshold <- NULL #SMI 2/9/14

           mydots$par <- unlist(new.start)
           mydots$hessian <- TRUE
           mydots$joint <- NULL
           mydots$upper <- unlist( upper[ nm[idx.measure] ])
           mydots$lower <- unlist( lower[ nm[idx.measure] ])
           mydots$method  <- method

           oout3 <- do.call(optim, args=mydots)

           theta3 <- oout3$par
           #print(theta3)
           HESS[measure.par,measure.par] <- oout3$hessian
           HaveMeasHess <- TRUE

           fixed <- old.fixed
           start <- old.start
           fixed.par <- old.fixed.par
       }
# END: ESTIMATION OF CP part



 if(!is.CARMA(yuima)){

  oout4 <- list(par=  c(theta1, theta2, theta3))
       names(oout4$par) <- c(diff.par,drift.par,measure.par)
       oout <- oout4
    }

     coef <- oout$par


       control=list()
       par <- coef
       if(!is.CARMA(yuima)){

       names(par) <- unique(c(diff.par, drift.par,measure.par))
       nm <- unique(c(diff.par, drift.par,measure.par))
       } else {
           names(par) <- unique(c(diff.par, drift.par))
           nm <- unique(c(diff.par, drift.par))
       }
#return(oout)


  if(is.CARMA(yuima) && length(yuima@model@parameter@measure)!=0){
      nm <-c(nm,measure.par)
      if((NoNeg.Noise==TRUE)){nm <-c(nm,mean.noise)}

      nm<-unique(nm)
  }
  if(is.CARMA(yuima) && (length(yuima@model@info@loc.par)!=0)){
    nm <-unique(c(nm,yuima@model@info@loc.par))
  }


	 conDrift <- list(trace = 5, fnscale = 1,
				 parscale = rep.int(5, length(drift.par)),
				 ndeps = rep.int(0.001, length(drift.par)), maxit = 100L,
				 abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1,
				 beta = 0.5, gamma = 2, REPORT = 10, type = 1, lmm = 5,
				 factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)
	 conDiff <- list(trace = 5, fnscale = 1,
				  parscale = rep.int(5, length(diff.par)),
				  ndeps = rep.int(0.001, length(diff.par)), maxit = 100L,
				  abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1,
				  beta = 0.5, gamma = 2, REPORT = 10, type = 1, lmm = 5,
				  factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)
      conMeas <- list(trace = 5, fnscale = 1,
                  parscale = rep.int(5, length(measure.par)),
                  ndeps = rep.int(0.001, length(measure.par)), maxit = 100L,
                  abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1,
                  beta = 0.5, gamma = 2, REPORT = 10, type = 1, lmm = 5,
                  factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)
  if(is.CARMA(yuima) && length(yuima@model@info@loc.par)!=0 ){
    conDrift <- list(trace = 5, fnscale = 1,
                     parscale = rep.int(5, length(c(drift.par,yuima@model@info@loc.par))),
                     ndeps = rep.int(0.001, length(c(drift.par,yuima@model@info@loc.par))),
                     maxit = 100L,
                     abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1,
                     beta = 0.5, gamma = 2, REPORT = 10, type = 1, lmm = 5,
                     factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)
    conDiff <- list(trace = 5, fnscale = 1,
                    parscale = rep.int(5, length(diff.par)),
                    ndeps = rep.int(0.001, length(diff.par)), maxit = 100L,
                    abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1,
                    beta = 0.5, gamma = 2, REPORT = 10, type = 1, lmm = 5,
                    factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)
  }



    if(!HaveDriftHess & (length(drift.par)>0)){
	  #hess2 <- .Internal(optimhess(coef[drift.par], fDrift, NULL, conDrift))
	   if(!is.CARMA(yuima)){
       hess2 <- optimHess(coef[drift.par], fDrift, NULL, control=conDrift)
  	  HESS[drift.par,drift.par] <- hess2
	   } else{
         names(coef) <- c(drift.par,yuima@model@info@loc.par)
	     hess2 <- optimHess(coef, fDrift, NULL, control=conDrift)
	     HESS <- hess2
	   }
     if(is.CARMA(yuima) && length(yuima@model@info@scale.par)!=0){
       b0<-paste0(yuima@model@info@ma.par,"0",collapse="")
       idx.b0<-match(b0,rownames(HESS))
       HESS<-HESS[-idx.b0,]
       HESS<-HESS[,-idx.b0]
     }
    }

	 if(!HaveDiffHess  & (length(diff.par)>0)){
	   hess1 <- optimHess(coef[diff.par], fDiff, NULL, control=conDiff)
		 HESS[diff.par,diff.par] <- hess1
	 }

	 oout$hessian <- HESS


     if(!HaveMeasHess & (length(measure.par)>0) & !is.CARMA(yuima)){
        hess1 <- optimHess(coef[measure.par], fMeas, NULL, control=conMeas)
        oout$hessian[measure.par,measure.par] <- hess1
     }

    vcov <- if (length(coef))
	  solve(oout$hessian)
    else matrix(numeric(0L), 0L, 0L)



    mycoef <- as.list(coef)

    if(!is.CARMA(yuima)){
        names(mycoef) <- nm
    }
    idx.fixed <- orig.idx.fixed



mycoef.cont <- mycoef
if(length(c(idx.fixed,idx.measure)>0))  # SMI 2/9/14
        mycoef.cont <- mycoef[-c(idx.fixed,idx.measure)]  # SMI 2/9/14


    min.diff <- 0
    min.jump <- 0


    if(length(c(diff.par,drift.par))>0 & !is.CARMA(yuima)){ # LM 04/09/14
	    min.diff <- minusquasilogl(yuima=yuima, param=mycoef[c(diff.par,drift.par)], print=print, env)
    }else{
      if(length(c(diff.par,drift.par))>0 & is.CARMA(yuima)){
        min.diff <- minusquasilogl(yuima=yuima, param=mycoef, print=print, env)
      }
    }

    if(length(c(measure.par))>0 & !is.CARMA(yuima))
        min.jump <-   minusquasipsi(yuima=yuima, param=mycoef[measure.par], print=print, env=env)



    min <- min.diff + min.jump
    if(min==0)
     min <- NA


  dummycov<-matrix(0,length(coef),length(coef))
  rownames(dummycov)<-names(coef)
  colnames(dummycov)<-names(coef)
  dummycov[rownames(vcov),colnames(vcov)]<-vcov
  vcov<-dummycov


#     new("mle", call = call, coef = coef, fullcoef = unlist(mycoef),
#        vcov = vcov, min = min, details = oout, minuslogl = minusquasilogl,
#        method = method)
#LM 11/01
  if(!is.CARMA(yuima)){
      if(length(measure.par)>0){
          final_res<-new("yuima.CP.qmle",
          Jump.times=env$time[env$Cn.r==0],
          Jump.values=env$deltaX[env$Cn.r==0,],
          X.values=env$X[env$Cn.r==0,],
          model=yuima@model,
          call = call, coef = coef, fullcoef = unlist(mycoef),
                   vcov = vcov, min = min, details = oout, minuslogl = minusquasilogl,
                   method = method, nobs=yuima.nobs, threshold=threshold)
      } else {
          final_res<-new("yuima.qmle", call = call, coef = coef, fullcoef = unlist(mycoef),
          vcov = vcov, min = min, details = oout, minuslogl = minusquasilogl,
          method = method, nobs=yuima.nobs, model=yuima@model)
      }
  } else {
    if( Est.Incr=="Carma.IncPar" || Est.Incr=="Carma.Inc" ){
    final_res<-new("yuima.carma.qmle", call = call, coef = coef, fullcoef = unlist(mycoef),
                   vcov = vcov, min = min, details = oout, minuslogl = minusquasilogl,
                    method = method, nobs=yuima.nobs, logL.Incr = NULL)
    }else{
      if(Est.Incr=="Carma.Par"){
      final_res<-new("mle", call = call, coef = coef, fullcoef = unlist(mycoef),
                     vcov = vcov, min = min, details = oout, minuslogl = minusquasilogl,
                     method = method, nobs=yuima.nobs)
      }else{
        yuima.warn("The variable Est.Incr is not correct. See qmle documentation for the allowed values ")
        final_res<-new("mle", call = call, coef = coef, fullcoef = unlist(mycoef),
                       vcov = vcov, min = min, details = oout, minuslogl = minusquasilogl,
                       method = method, nobs=yuima.nobs)
        return(final_res)
      }
    }
  }

if(!is.CARMA(yuima)){
    return(final_res)
 }else {

    param<-coef(final_res)

    observ<-yuima@data
    model<-yuima@model
    info<-model@info

    numb.ar<-info@p
    name.ar<-paste(info@ar.par,c(numb.ar:1),sep="")
    ar.par<-param[name.ar]

    numb.ma<-info@q
    name.ma<-paste(info@ma.par,c(0:numb.ma),sep="")
    ma.par<-param[name.ma]

    loc.par=NULL
    if (length(info@loc.par)!=0){
      loc.par<-param[info@loc.par]
    }

    scale.par=NULL
    if (length(info@scale.par)!=0){
      scale.par<-param[info@scale.par]
    }

    lin.par=NULL
    if (length(info@lin.par)!=0){
      lin.par<-param[info@lin.par]
    }
    if(min(yuima.PhamBreton.Alg(ar.par[numb.ar:1]))>=0){
      cat("\n Stationarity condition is satisfied...\n Starting Estimation Increments ...\n")
    }else{
      yuima.warn("Insert constraints in Autoregressive parameters for enforcing stationarity" )
      cat("\n Starting Estimation Increments ...\n")
    }

    ttt<-observ@zoo.data[[1]]
    tt<-index(ttt)
    y<-coredata(ttt)
    if(NoNeg.Noise==TRUE && (info@p==(info@q+1))){final_res@coef[mean.noise]<-mean(y)/tail(ma.par,n=1)*ar.par[1]}

    levy<-yuima.CarmaNoise(y,tt,ar.par,ma.par, loc.par, scale.par, lin.par, NoNeg.Noise)
    inc.levy<-NULL
    if (!is.null(levy)){
      inc.levy<-diff(t(levy))
    }
    # INSERT HERE THE NECESSARY STEPS FOR FINDING THE PARAMETERS OF LEVY
   if(Est.Incr=="Carma.Inc"){
     # inc.levy.fin<-zoo(inc.levy,tt,frequency=1/env$h)
     inc.levy.fin<-zoo(inc.levy,tt[(1+length(tt)-length(inc.levy)):length(tt)])
     carma_final_res<-new("yuima.carma.qmle", call = call, coef = coef, fullcoef = unlist(mycoef),
                          vcov = vcov, min = min, details = oout, minuslogl = minusquasilogl,
                          method = method, Incr.Lev = inc.levy.fin,
                          model = yuima@model, nobs=yuima.nobs, logL.Incr = NULL)
     return(carma_final_res)
   }

   cat("\nStarting Estimation parameter Noise ...\n")

    dummycovCarmapar<-vcov[unique(c(drift.par,diff.par)),unique(c(drift.par,diff.par))]
    if(!is.null(loc.par)){
      dummycovCarmapar<-vcov[unique(c(drift.par,diff.par,info@loc.par)),
                             unique(c(drift.par,diff.par,info@loc.par))]
    }



dummycovCarmaNoise<-vcov[unique(measure.par),unique(c(measure.par))] #we need to adjusted
    dummycoeffCarmapar<-coef[unique(c(drift.par,diff.par))]
    if(!is.null(loc.par)){
      dummycoeffCarmapar<-coef[unique(c(drift.par,diff.par,info@loc.par))]
    }

    dummycoeffCarmaNoise<-coef[unique(c(measure.par))]
     coef<-NULL
    coef<-c(dummycoeffCarmapar,dummycoeffCarmaNoise)
    names.par<-c(unique(c(drift.par,diff.par)),unique(c(measure.par)))
    if(!is.null(loc.par)){
      names.par<-c(unique(c(drift.par,diff.par,info@loc.par)),unique(c(measure.par)))
    }

    names(coef)<-names.par
    cov<-NULL
    cov<-matrix(0,length(names.par),length(names.par))
    rownames(cov)<-names.par
    colnames(cov)<-names.par
    if(is.null(loc.par)){
      cov[unique(c(drift.par,diff.par)),unique(c(drift.par,diff.par))]<-dummycovCarmapar
    }else{
      cov[unique(c(drift.par,diff.par,info@loc.par)),unique(c(drift.par,diff.par,info@loc.par))]<-dummycovCarmapar
    }

    cov[unique(c(measure.par)),unique(c(measure.par))]<-dummycovCarmaNoise

    if(length(model@measure.type)!=0){
      if(model@measure.type=="CP"){
        name.func.dummy <- as.character(model@measure$df$expr[1])
        name.func<- substr(name.func.dummy,1,(nchar(name.func.dummy)-1))
        names.measpar<-as.vector(strsplit(name.func,', '))[[1]][-1]
        valuemeasure<-as.numeric(names.measpar)
        name.int.dummy <- as.character(model@measure$intensity)
        valueintensity<-as.numeric(name.int.dummy)
        NaIdx<-which(!is.na(c(valueintensity,valuemeasure)))

        if(length(NaIdx)!=0){
          yuima.warn("the constrained MLE for levy increment will be implemented as soon as possible")
          carma_final_res<-new("yuima.carma.qmle", call = call, coef = coef, fullcoef = unlist(mycoef),
                               vcov = vcov, min = min, details = oout, minuslogl = minusquasilogl,
                               method = method, Incr.Lev = inc.levy,
                               model = yuima@model, logL.Incr = NULL)
          return(carma_final_res)
        }

        if(aggregation==TRUE){
          if(floor(yuima@sampling@n/yuima@sampling@Terminal)!=yuima@sampling@n/yuima@sampling@Terminal){
            yuima.stop("the n/Terminal in sampling information is not an integer. Set Aggregation=FALSE")
          }
          inc.levy1<-diff(cumsum(c(0,inc.levy))[seq(from=1,
                                               to=yuima@sampling@n[1],
                                               by=(yuima@sampling@n/yuima@sampling@Terminal)[1]
                                               )])
        }else{
          inc.levy1<-inc.levy
        }

        names.measpar<-c(name.int.dummy, names.measpar)

        if(measurefunc=="dnorm"){

#           result.Lev<-yuima.Estimation.CPN(Increment.lev=inc.levy1,param0=coef[ names.measpar],
#                                            fixed.carma=fixed.carma,
#                                            lower.carma=lower.carma,
#                                            upper.carma=upper.carma)

          result.Lev<-yuima.Estimation.Lev(Increment.lev=inc.levy1,
                                           param0=coef[ names.measpar],
                                           fixed.carma=fixed.carma,
                                           lower.carma=lower.carma,
                                           upper.carma=upper.carma,
                                           measure=measurefunc,
                                           measure.type=model@measure.type,
                                           dt=env$h,
                                           aggregation=aggregation)

        }
        if(measurefunc=="dgamma"){
#           result.Lev<-yuima.Estimation.CPGam(Increment.lev=inc.levy1,param0=coef[ names.measpar],
#                                              fixed.carma=fixed.carma,
#                                              lower.carma=lower.carma,
#                                              upper.carma=upper.carma)

          result.Lev<-yuima.Estimation.Lev(Increment.lev=inc.levy1,
                                           param0=coef[ names.measpar],
                                           fixed.carma=fixed.carma,
                                           lower.carma=lower.carma,
                                           upper.carma=upper.carma,
                                           measure=measurefunc,
                                           measure.type=model@measure.type,
                                           dt=env$h,
                                           aggregation=aggregation)
        }
        if(measurefunc=="dexp"){
#           result.Lev<-yuima.Estimation.CPExp(Increment.lev=inc.levy1,param0=coef[ names.measpar],
#                                              fixed.carma=fixed.carma,
#                                              lower.carma=lower.carma,
#                                              upper.carma=upper.carma)

          result.Lev<-yuima.Estimation.Lev(Increment.lev=inc.levy1,
                                           param0=coef[ names.measpar],
                                           fixed.carma=fixed.carma,
                                           lower.carma=lower.carma,
                                           upper.carma=upper.carma,
                                           measure=measurefunc,
                                           measure.type=model@measure.type,
                                           dt=env$h,
                                           aggregation=aggregation)

        }
        Inc.Parm<-result.Lev$estLevpar
        IncVCOV<-result.Lev$covLev

        names(Inc.Parm)[NaIdx]<-measure.par
        rownames(IncVCOV)[NaIdx]<-as.character(measure.par)
        colnames(IncVCOV)[NaIdx]<-as.character(measure.par)

        coef<-NULL
        coef<-c(dummycoeffCarmapar,Inc.Parm)

        names.par<-names(coef)
        cov<-NULL
        cov<-matrix(0,length(names.par),length(names.par))
        rownames(cov)<-names.par
        colnames(cov)<-names.par
        if(is.null(loc.par)){
          cov[unique(c(drift.par,diff.par)),unique(c(drift.par,diff.par))]<-dummycovCarmapar
        }else{
          cov[unique(c(drift.par,diff.par,info@loc.par)),unique(c(drift.par,diff.par,info@loc.par))]<-dummycovCarmapar
        }
        cov[names(Inc.Parm),names(Inc.Parm)]<-IncVCOV


      }
      if(yuima@model@measure.type=="code"){
  #     #  "rIG", "rNIG", "rgamma", "rbgamma", "rngamma"
        name.func.dummy <- as.character(model@measure$df$expr[1])
        name.func<- substr(name.func.dummy,1,(nchar(name.func.dummy)-1))
        names.measpar<-as.vector(strsplit(name.func,', '))[[1]][-1]
        valuemeasure<-as.numeric(names.measpar)
        NaIdx<-which(!is.na(valuemeasure))
        if(length(NaIdx)!=0){
          yuima.warn("the constrained MLE for levy increment will be implemented as soon as possible")
          carma_final_res<-new("yuima.carma.qmle", call = call, coef = coef, fullcoef = unlist(mycoef),
                               vcov = vcov, min = min, details = oout, minuslogl = minusquasilogl,
                               method = method, Incr.Lev = inc.levy,
                               model = yuima@model, logL.Incr = NULL)
          return(carma_final_res)
        }
        if(aggregation==TRUE){
          if(floor(yuima@sampling@n/yuima@sampling@Terminal)!=yuima@sampling@n/yuima@sampling@Terminal){
            yuima.stop("the n/Terminal in sampling information is not an integer. Aggregation=FALSE is recommended")
          }
         inc.levy1<-diff(cumsum(c(0,inc.levy))[seq(from=1,
                                              to=yuima@sampling@n[1],
                                              by=(yuima@sampling@n/yuima@sampling@Terminal)[1]
         )])
        }else{
        inc.levy1<-inc.levy
        }

        if(measurefunc=="rIG"){

#           result.Lev<-list(estLevpar=coef[ names.measpar],
#                            covLev=matrix(NA,
#                                          length(coef[ names.measpar]),
#                                          length(coef[ names.measpar]))
#           )
#           result.Lev<-yuima.Estimation.IG(Increment.lev=inc.levy1,param0=coef[ names.measpar],
#                                           fixed.carma=fixed.carma,
#                                           lower.carma=lower.carma,
#                                           upper.carma=upper.carma)

          result.Lev<-yuima.Estimation.Lev(Increment.lev=inc.levy1,
                                           param0=coef[ names.measpar],
                                           fixed.carma=fixed.carma,
                                           lower.carma=lower.carma,
                                           upper.carma=upper.carma,
                                           measure=measurefunc,
                                           measure.type=model@measure.type,
                                           dt=env$h,
                                           aggregation=aggregation)
          #         result.Levy<-gigFit(inc.levy)
  #         Inc.Parm<-coef(result.Levy)
  #         IncVCOV<--solve(gigHessian(inc.levy, param=Inc.Parm))
        }
        if(measurefunc=="rNIG"){
#            result.Lev<-yuima.Estimation.NIG(Increment.lev=inc.levy1,param0=coef[ names.measpar],
#                                             fixed.carma=fixed.carma,
#                                             lower.carma=lower.carma,
#                                             upper.carma=upper.carma)

          result.Lev<-yuima.Estimation.Lev(Increment.lev=inc.levy1,
                                           param0=coef[ names.measpar],
                                           fixed.carma=fixed.carma,
                                           lower.carma=lower.carma,
                                           upper.carma=upper.carma,
                                           measure=measurefunc,
                                           measure.type=model@measure.type,
                                           dt=env$h,
                                           aggregation=aggregation)
        }
        if(measurefunc=="rbgamma"){
          result.Lev<-list(estLevpar=coef[ names.measpar],
                           covLev=matrix(NA,
                                         length(coef[ names.measpar]),
                                         length(coef[ names.measpar]))
                           )
        }
        if(measurefunc=="rngamma"){
#           result.Lev<-yuima.Estimation.VG(Increment.lev=inc.levy1,param0=coef[ names.measpar],
#                                           fixed.carma=fixed.carma,
#                                           lower.carma=lower.carma,
#                                           upper.carma=upper.carma)

          result.Lev<-yuima.Estimation.Lev(Increment.lev=inc.levy1,
                                           param0=coef[ names.measpar],
                                           fixed.carma=fixed.carma,
                                           lower.carma=lower.carma,
                                           upper.carma=upper.carma,
                                           measure=measurefunc,
                                           measure.type=model@measure.type,
                                           dt=env$h,
                                           aggregation=aggregation)

        }

        Inc.Parm<-result.Lev$estLevpar
        IncVCOV<-result.Lev$covLev

        names(Inc.Parm)[NaIdx]<-measure.par
        rownames(IncVCOV)[NaIdx]<-as.character(measure.par)
        colnames(IncVCOV)[NaIdx]<-as.character(measure.par)

        coef<-NULL
        coef<-c(dummycoeffCarmapar,Inc.Parm)

        names.par<-names(coef)
        cov<-NULL
        cov<-matrix(0,length(names.par),length(names.par))
        rownames(cov)<-names.par
        colnames(cov)<-names.par
        if(is.null(loc.par)){
          cov[unique(c(drift.par,diff.par)),unique(c(drift.par,diff.par))]<-dummycovCarmapar
        }else{
          cov[unique(c(drift.par,diff.par,info@loc.par)),unique(c(drift.par,diff.par,info@loc.par))]<-dummycovCarmapar
        }
        cov[names(Inc.Parm),names(Inc.Parm)]<-IncVCOV

      }
    }
#     dummycovCarmapar<-vcov[unique(c(drift.par,diff.par)),unique(c(drift.par,diff.par))]
#     dummycovCarmaNoise<-vcov[unique(measure.par),unique(c(measure.par))] #we need to adjusted
#     dummycoeffCarmapar<-coef[unique(c(drift.par,diff.par))]
#     dummycoeffCarmaNoise<-coef[unique(c(measure.par))]
#     coef<-NULL
#     coef<-c(dummycoeffCarmapar,dummycoeffCarmaNoise)
#     names.par<-c(unique(c(drift.par,diff.par)),unique(c(measure.par)))
#     names(coef)<-names.par
#     cov<-NULL
#     cov<-matrix(0,length(names.par),length(names.par))
#     rownames(cov)<-names.par
#     colnames(cov)<-names.par
#     cov[unique(c(drift.par,diff.par)),unique(c(drift.par,diff.par))]<-dummycovCarmapar
#     cov[unique(c(measure.par)),unique(c(measure.par))]<-dummycovCarmaNoise

#    carma_final_res<-list(mle=final_res,Incr=inc.levy,model=yuima)
    if(Est.Incr=="Carma.IncPar"){
      #inc.levy.fin<-zoo(inc.levy,tt,frequency=1/env$h)
      inc.levy.fin<-zoo(inc.levy,tt[(1+length(tt)-length(inc.levy)):length(tt)])
      carma_final_res<-new("yuima.carma.qmle", call = call, coef = coef, fullcoef = unlist(coef),
                     vcov = cov, min = min, details = oout, minuslogl = minusquasilogl,
                     method = method, Incr.Lev = inc.levy.fin,
                           model = yuima@model, nobs=yuima.nobs,
                     logL.Incr = tryCatch(-result.Lev$value,error=function(theta){NULL}))
    }else{
      if(Est.Incr=="Carma.Par"){
        carma_final_res<-new("mle", call = call, coef = coef, fullcoef = unlist(coef),
            vcov = cov, min = min, details = oout, minuslogl = minusquasilogl,
            method = method, nobs=yuima.nobs)
      }
    }
    return(carma_final_res)
  }
}

# SMI-2/9/14 CP
minusquasipsi <- function(yuima, param, print=FALSE, env){

    idx.intensity <- env$idx.intensity

    fullcoef <- yuima@model@parameter@all
    measurecoef <- param[yuima@model@parameter@measure]

    npar <- length(fullcoef)
    nm <- names(param)
    oo <- match(nm, fullcoef)

    if(any(is.na(oo)))
        yuima.stop("some named arguments in 'param' are not arguments to the supplied yuima model")
    param <- param[order(oo)]

    h <- env$h
    Dn.r <- !env$Cn.r

    #    if(length(idx.intensity)){
    #    intensity <- unlist(measurecoef[idx.intensity])
    #}else{
    #    intensity <- eval(yuima@model@measure$intensity, envir=env)
    #}

    #	print(intensity)
    #print(str(env$time))

#  tmp.env <- new.env()
#for(i in 1:length(param)){
#    assign(names(param)[i],param[[i]],envir=tmp.env)
#}
#print(ls(env))

    d.size <- yuima@model@equation.number
    n <- length(yuima)[1]
    myidx <- which(Dn.r)[-n]

    measure <- measure.term(yuima, param, env)

    QL <- 0

    dx <- env$deltaX
    measure.var <- env$measure.var

    for(i in 1:length(measurecoef))
    #if(!is.Poisson(yuima)){
        #      if(is.na(match(i,idx.intensity)))
        #   assign(names(measurecoef)[i],measurecoef[i][[1]], envir=env)
        # } else {
        assign(names(measurecoef)[i],measurecoef[i][[1]], envir=env)
        # }

    #    print("### ls(env)")
    #       print(ls(env))
    if(is.null(dim(measure[,,1]))){  # one-dimensional
        for(t in myidx){
            iC <- 1/measure[, , t]
            assign(measure.var,iC%*%dx[t,],envir=env)
            assign(yuima@model@time.variable, env$time[t], envir=env)
            #       print("### t")
            #print(t)
            #print(env$time[t])
            intensity <- eval(yuima@model@measure$intensity, envir=env)
            #print("intensity")
            #print(intensity)
            dF <- intensity*eval(yuima@model@measure$df$expr,envir=env)/iC
            logpsi <- 0
            if(dF>0)
                logpsi <- log(dF)
            QL <- QL + logpsi
        }
    } else {
        for(t in myidx){
            iC <- solve(measure[, , t])
            assign(measure.var,iC%*%dx[t,], envir=env)
            assign(yuima@model@time.variable, env$time[t], envir=env)
            intensity <- eval(yuima@model@measure$intensity, envir=env)
            dF <- intensity*eval(yuima@model@measure$df$expr,envir=env)*det(iC)
            logpsi <- 0
            if(dF>0)
                logpsi <- log(dF)
            QL <- QL + logpsi
        }
    }

    myf <- function(x) {
        f1 <- function(u){
         assign(yuima@model@time.variable, u, envir=env)
         intensity <- eval(yuima@model@measure$intensity, envir=env)
        }
        sapply(x, f1)
    }
    #    print(myf(1))
    #  print(str( try(integrate(f=myf, lower=yuima@sampling@Initial, upper=yuima@sampling@Terminal,subdivisions=100),silent=TRUE )))

    myint <- integrate(f=myf, lower=yuima@sampling@Initial, upper=yuima@sampling@Terminal,subdivisions=100)$value
    #  print(myint)
    #print(-h*intensity*(n-1))
    #    QL <- QL -h*intensity*(n-1)
    QL <- QL -myint


    if(!is.finite(QL)){
        yuima.warn("quasi likelihood is too small to calculate.")
        return(1e10)
    }
    if(print==TRUE){
        yuima.warn(sprintf("NEG-QL: %f, %s", -QL, paste(names(param),param,sep="=",collapse=", ")))
    }
    if(is.infinite(QL)) return(1e10)
    return(as.numeric(-QL))

}


quasilogl <- function(yuima, param, print=FALSE){

	d.size <- yuima@model@equation.number
	if (is(yuima@model, "yuima.carma")){
	  # 24/12
	  d.size <-1
	}

	n <- length(yuima)[1]

	env <- new.env()
    assign("X",  as.matrix(onezoo(yuima)), envir=env)
    assign("deltaX",  matrix(0, n-1, d.size), envir=env)
    assign("Cn.r", rep(1,n-1), envir=env)

    if(is.CARMA(yuima)){
        env$X<-as.matrix(env$X[,1])
        env$deltaX<-as.matrix(env$deltaX[,1])
        env$time.obs<-length(env$X)
        assign("p", yuima@model@info@p, envir=env)
        assign("q", yuima@model@info@q, envir=env)
        assign("V_inf0", matrix(diag(rep(1,env$p)),env$p,env$p), envir=env)
	}


	for(t in 1:(n-1))
        env$deltaX[t,] <- env$X[t+1,] - env$X[t,]

	assign("h", deltat(yuima@data@zoo.data[[1]]), envir=env)
	assign("time", as.numeric(index(yuima@data@zoo.data[[1]])), envir=env)

	-minusquasilogl(yuima=yuima, param=param, print=print, env)
}


minusquasilogl <- function(yuima, param, print=FALSE, env){

	diff.par <- yuima@model@parameter@diffusion

	drift.par <- yuima@model@parameter@drift
		if(is.CARMA(yuima)){
		  if(length(yuima@model@info@scale.par)!=0){
	      xinit.par <- yuima@model@parameter@xinit
		  }
		}


    if(is.CARMA(yuima) && length(yuima@model@info@lin.par)==0
	   && length(yuima@model@parameter@jump)!=0){
	  diff.par<-yuima@model@parameter@jump
	 # measure.par<-yuima@model@parameter@measure
	}

	if(is.CARMA(yuima) && length(yuima@model@info@lin.par)==0
	   && length(yuima@model@parameter@measure)!=0){
	  measure.par<-yuima@model@parameter@measure
	}

	# 24/12
	if(is.CARMA(yuima) && length(yuima@model@info@lin.par)>0  ){
	  yuima.warn("carma(p,q): the case of lin.par will be implemented as soon as")
	  return(NULL)
	}

	if(is.CARMA(yuima)){
	   xinit.par <- yuima@model@parameter@xinit
	}


    drift.par <- yuima@model@parameter@drift

	fullcoef <- NULL

	if(length(diff.par)>0)
	fullcoef <- diff.par

	if(length(drift.par)>0)
        fullcoef <- c(fullcoef, drift.par)

    if(is.CARMA(yuima)){
	    if(length(xinit.par)>0)
	        fullcoef <- c(fullcoef, xinit.par)
    }

  	if(is.CARMA(yuima) && (length(yuima@model@parameter@measure)!=0))
        fullcoef<-c(fullcoef, measure.par)

    if(is.CARMA(yuima)){
        if("mean.noise" %in% names(param)){
            mean.noise<-"mean.noise"
            fullcoef <- c(fullcoef, mean.noise)
            NoNeg.Noise<-TRUE
        }
    }


    npar <- length(fullcoef)

    nm <- names(param)
    oo <- match(nm, fullcoef)

    if(any(is.na(oo)))
        yuima.stop("some named arguments in 'param' are not arguments to the supplied yuima model")
    param <- param[order(oo)]
    nm <- names(param)

	idx.diff <- match(diff.par, nm)
	idx.drift <- match(drift.par, nm)


	if(is.CARMA(yuima)){
	    idx.xinit <-as.integer(na.omit(match(xinit.par, nm)))
	}

	h <- env$h

    Cn.r <- env$Cn.r

    theta1 <- unlist(param[idx.diff])
    theta2 <- unlist(param[idx.drift])


	n.theta1 <- length(theta1)
	n.theta2 <- length(theta2)
	n.theta <- n.theta1+n.theta2


	if(is.CARMA(yuima)){
	    theta3 <- unlist(param[idx.xinit])
	    n.theta3 <- length(theta3)
	    n.theta <- n.theta1+n.theta2+n.theta3
	}


  d.size <- yuima@model@equation.number


	n <- length(yuima)[1]


  if (is.CARMA(yuima)){
	  # 24/12
	  d.size <-1
    # We build the two step procedure as described in
  #  if(length(yuima@model@info@scale.par)!=0){
       prova<-as.numeric(param)
       #names(prova)<-fullcoef[oo]
	     names(prova)<-names(param)
       param<-prova[c(length(prova):1)]
       time.obs<-env$time.obs
       y<-as.numeric(env$X)
       u<-env$h
       p<-env$p
       q<-env$q
#         p<-yuima@model@info@p
	  ar.par <- yuima@model@info@ar.par
	  name.ar<-paste0(ar.par, c(1:p))
# 	  q <- yuima@model@info@q
	  ma.par <- yuima@model@info@ma.par
	  name.ma<-paste0(ma.par, c(0:q))
       if (length(yuima@model@info@loc.par)==0){

         a<-param[name.ar]
  #        a_names<-names(param[c(1:p)])
  #        names(a)<-a_names
         b<-param[name.ma]
  #        b_names<-names(param[c((p+1):(length(param)-p+1))])
  #        names(b)<-b_names
         if(length(yuima@model@info@scale.par)!=0){
          if(length(b)==1){
             b<-1
          } else{
             indx_b0<-paste0(yuima@model@info@ma.par,"0",collapse="")
             b[indx_b0]<-1
          }
       sigma<-tail(param,1)
         }else {sigma<-1}
         NoNeg.Noise<-FALSE
         if(is.CARMA(yuima)){
           if("mean.noise" %in% names(param)){

             NoNeg.Noise<-TRUE
           }
         }
         if(NoNeg.Noise==TRUE){
           if (length(b)==p){
             #mean.noise<-param[mean.noise]
             # Be useful for carma driven by a no negative levy process
             mean.y<-mean(y)
             #mean.y<-mean.noise*tail(b,n=1)/tail(a,n=1)*sigma
             #param[mean.noise]<-mean.y/(tail(b,n=1)/tail(a,n=1)*sigma)
           }else{
             mean.y<-0
           }
           y<-y-mean.y
         }
       # V_inf0<-matrix(diag(rep(1,p)),p,p)
        V_inf0<-env$V_inf0
        p<-env$p
        q<-env$q
        strLog<-yuima.carma.loglik1(y, u, a, b, sigma,time.obs,V_inf0,p,q)
       }else{
         # 01/01
#          ar.par <- yuima@model@info@ar.par
#          name.ar<-paste0(ar.par, c(1:p))
         a<-param[name.ar]
#          ma.par <- yuima@model@info@ma.par
#          q <- yuima@model@info@q
         name.ma<-paste0(ma.par, c(0:q))
         b<-param[name.ma]
         if(length(yuima@model@info@scale.par)!=0){
            if(length(b)==1){
                b<-1
              } else{
               indx_b0<-paste0(yuima@model@info@ma.par,"0",collapse="")
               b[indx_b0]<-1
              }
              scale.par <- yuima@model@info@scale.par
              sigma <- param[scale.par]
         } else{sigma <- 1}
         loc.par <- yuima@model@info@loc.par
         mu <- param[loc.par]

         NoNeg.Noise<-FALSE
         if(is.CARMA(yuima)){
           if("mean.noise" %in% names(param)){

             NoNeg.Noise<-TRUE
           }
         }

# Lines 883:840 work if we have a no negative noise
        if(is.CARMA(yuima)&&(NoNeg.Noise==TRUE)){
           if (length(b)==p){
             mean.noise<-param[mean.noise]
           # Be useful for carma driven by levy process
          #   mean.y<-mean.noise*tail(b,n=1)/tail(a,n=1)*sigma
             mean.y<-mean(y-mu)

           }else{
             mean.y<-0
           }
           y<-y-mean.y
        }


         y.start <- y-mu
         #V_inf0<-matrix(diag(rep(1,p)),p,p)
         V_inf0<-env$V_inf0
         p<-env$p
         q<-env$q
         strLog<-yuima.carma.loglik1(y.start, u, a, b, sigma,time.obs,V_inf0,p,q)
       }

       QL<-strLog$loglikCdiag
#       }else {
#         yuima.warn("carma(p,q): the scale parameter is equal to 1. We will implemented as soon as possible")
#         return(NULL)
#     }
	} else{
  	drift <- drift.term(yuima, param, env)
  	diff <- diffusion.term(yuima, param, env)

  	QL <- 0

  	pn <- 0


  	vec <- env$deltaX-h*drift[-n,]

  	K <- -0.5*d.size * log( (2*pi*h) )

  	dimB <- dim(diff[, , 1])

   if(is.null(dimB)){  # one dimensional X
  	  for(t in 1:(n-1)){
  		yB <- diff[, , t]^2
  		logdet <- log(yB)
  		pn <- Cn.r[t]*(K - 0.5*logdet-0.5*vec[t, ]^2/(h*yB))
  		QL <- QL+pn

  		}
  	} else {  # multidimensional X
  	 for(t in 1:(n-1)){
  		yB <- diff[, , t] %*% t(diff[, , t])
  		logdet <- log(det(yB))
  		if(is.infinite(logdet) ){ # should we return 1e10?
  			pn <- log(1)
  			yuima.warn("singular diffusion matrix")
  			return(1e10)
  		}else{
  			pn <- (K - 0.5*logdet +
  					  ((-1/(2*h))*t(vec[t, ])%*%solve(yB)%*%vec[t, ]))*Cn.r[t]
  			QL <- QL+pn
  		}
  	 }
  	}
  }


	if(!is.finite(QL)){
		yuima.warn("quasi likelihood is too small to calculate.")
		return(1e10)
	}
	if(print==TRUE){
		yuima.warn(sprintf("NEG-QL: %f, %s", -QL, paste(names(param),param,sep="=",collapse=", ")))
	}
	if(is.infinite(QL)) return(1e10)
	return(as.numeric(-QL))

}




MatrixA<-function (a)
{
  #Build Matrix A in the state space representation of Carma(p,q)
  #given the autoregressive coefficient
  pp = length(a)
  af = cbind(rep(0, pp - 1), diag(pp - 1))
  af = rbind(af, -a[pp:1])
  return(af)
}


# yuima.Vinfinity<-function(elForVInf,v){
#   # We find the infinity stationary variance-covariance matrix
#   A<-elForVInf$A
#   sigma<-elForVInf$sigma
# #   #p<-dim(A)[1]
# #   p<-elForVInf$p
#   ATrans<-elForVInf$ATrans
#   matrixV<-elForVInf$matrixV
#   matrixV[upper.tri(matrixV,diag=TRUE)]<-v
#   matrixV<-as.matrix(forceSymmetric(matrixV))
# #matrixV[lower.tri(matrixV)]<-matrixV[upper.tri(matrixV)]
# #  l<-rbind(matrix(rep(0,p-1),p-1,1),1)
# #  matrixV<-matrix(v,p,p)
#
#   lTrans<-elForVInf$lTrans
#   l<-elForVInf$l
#
#
#   RigSid<-l%*%elForVInf$lTrans
#   Matrixobj<-A%*%matrixV+matrixV%*%ATrans+sigma^2*RigSid
#   obj<-sum(Matrixobj^2)
#   obj
# }


#carma.kalman<-function(y, tt, p, q, a,bvector, sigma){
carma.kalman<-function(y, u, p, q, a,bvector, sigma, times.obs, V_inf0){
  #new Code
  A<-MatrixA(a)
  expA<-expm(A*u,method="Pade",order=6, trySym=FALSE, do.sparseMsg = FALSE)

  V_inf<-V0inf(a,p,sigma)

  expAT<-t(expA)

  Qmatr <- V_inf - expA %*% V_inf %*% expAT

  statevar<-numeric(length=p)

  SigMatr <- V_inf+0

  sd_2<-0
  Result<-numeric(length=2)
  Kgain<-numeric(length=p)
  dum_zc<-numeric(length=p)
  Mat22int<-numeric(length=(p*p))

  loglstar<- .Call("Cycle_Carma", y, statevar, expA, as.integer(length(y)),
                    as.integer(p), Qmatr, SigMatr, bvector, Result, Kgain,
                   dum_zc, Mat22int,
                    PACKAGE="yuima")
  return(list(loglstar=loglstar[1]-0.5*log(2*pi)*times.obs,s2hat=loglstar[2]))

#   # Old version
#
#
#   V_inf0<-matrix(diag(rep(1,p)),p,p)
#
#   A<-MatrixA(a)
#   # u<-diff(tt)[1]
#
#
#   #  Amatx<-yuima.carma.eigen(A)
#   #   expA<-Amatx$vectors%*%expm(diag(Amatx$values*u),
#   #                              method="Pade",
#   #                              order=6,
#   #                              trySym=TRUE,
#   #                              do.sparseMsg = TRUE)%*%solve(Amatx$vectors)
#
#   #   if(!is.complex(Amatx$values)){
#   #     expA<-Amatx$vectors%*%diag(exp(Amatx$values*u))%*%solve(Amatx$vectors)
#   #     }else{
#   expA<-expm(A*u,method="Pade",order=6, trySym=FALSE, do.sparseMsg = FALSE)
#   #    }
#   #expA<-yuima.exp(A*u)
#
#   v<-as.numeric(V_inf0[upper.tri(V_inf0,diag=TRUE)])
#
#   ATrans<-t(A)
#   matrixV<-matrix(0,p,p)
#   #l.dummy<-c(rep(0,p-1),1)
#   l<-rbind(matrix(rep(0,p-1),p-1,1),1)
#   #l<-matrix(l.dummy,p,1)
#   #lTrans<-matrix(l.dummy,1,p)
#   lTrans<-t(l)
#   elForVInf<-new.env()
#   elForVInf$A<-A
#   elForVInf$ATrans<-ATrans
#   elForVInf$lTrans<-lTrans
#   elForVInf$l<-l
#   elForVInf$matrixV<-matrixV
#   elForVInf$sigma<-sigma
#   #   elForVInf<-list(A=A,
#   #                   ATrans=ATrans,
#   #                   lTrans=lTrans,
#   #                   l=l,
#   #                   matrixV=matrixV,
#   #                   sigma=sigma)
#   #
#   V_inf_vect<-nlm(yuima.Vinfinity, v,
#                   elForVInf = elForVInf)$estimate
#   #  V_inf_vect<-nlminb(start=v,objective=yuima.Vinfinity, elForVInf = elForVInf)$par
#   #  V_inf_vect<-optim(par=v,fn=yuima.Vinfinity,method="L-BFGS-B", elForVInf = elForVInf)$par
#   V_inf<-matrix(0,p,p)
#
#   V_inf[upper.tri(V_inf,diag=TRUE)]<-V_inf_vect
#   V_inf<-forceSymmetric(V_inf)
#
#   V_inf[abs(V_inf)<= 1.e-06]=0
#
# #      A<-MatrixA(a)
# #      expA<-expm(A*u,method="Pade",order=6, trySym=FALSE, do.sparseMsg = FALSE)
# #
# #      V_inf<-V0inf(a,p,sigma)
#   #
#
#
#   expAT<-t(expA)
#   #SIGMA_err<-V_inf-expA%*%V_inf%*%t(expA)
#   SigMatr<-V_inf-expA%*%V_inf%*%expAT
#   statevar<-matrix(rep(0, p),p,1)
#   Qmatr<-SigMatr
#
#   # set
#   #statevar<-statevar0
#
#   # SigMatr<-expA%*%V_inf%*%t(expA)+Qmatr
#
#   #SigMatr<-Qmatr
#   SigMatr<-V_inf
#
#   zc<-matrix(bvector,1,p)
#   loglstar <- 0
#   loglstar1 <- 0
#
#   #  zcT<-matrix(bvector,p,1)
#   zcT<-t(zc)
#   for(t in 1:times.obs){
#     # prediction
#     statevar<-expA%*%statevar
#     SigMatr<-expA%*%SigMatr%*%expAT+Qmatr
#     # forecast
#     Uobs<-y[t]-zc%*%statevar
#     dum.zc<-zc%*%SigMatr
#     sd_2<-dum.zc%*%zcT
#     # sd_2<-zc%*%SigMatr%*%zcT
#     Inv_sd_2<-1/sd_2
#     #correction
#     Kgain<-SigMatr%*%zcT%*%Inv_sd_2
#     statevar<-statevar+Kgain%*%Uobs
#     #SigMatr<-SigMatr-Kgain%*%zc%*%SigMatr
#     SigMatr<-SigMatr-Kgain%*%dum.zc
#     term_int<--0.5*(log(sd_2)+Uobs%*%Uobs%*%Inv_sd_2)
#     loglstar<-loglstar+term_int
#   }
#   return(list(loglstar=(loglstar-0.5*log(2*pi)*times.obs),s2hat=sd_2))
}

V0inf<-function(a,p,sigma){
  # This code is based on the paper A continuous-time ARMA process Tsai-Chan 2000
  # we need to find the values along the diagonal
  #l<-c(numeric(length=(p-1)),0.5)
  # B_{p*p}V^{*}_{p*1}=-sigma^2*l/2
  B<-matrix(0,nrow=p,ncol=p)
  aa <- -rev(a)
#   B1<-.Call("Coeffdiag_B", as.integer(p), aa, B,
#         PACKAGE="yuima")
#   B<-matrix(0,nrow=p,ncol=p)
  for(i in 1:p){
    # Condition on B
    for(j in 1:p){
      if ((2*j-i) %in% c(1:p)){
        B[i,j]<-(-1)^(j-i)*aa[2*j-i]
      }
      if((2*j-i)==(p+1)){
        B[i,j]<-(-1)^(j-i-1)
      }
    }
  }
  Vdiag <- -solve(B)[,p]*0.5*sigma^2
  V <- diag(Vdiag)
  # we insert the values outside the diagonal
  for(i in 1:p){
    for(j in (i+1):p){
      if((i+j)  %% 2 == 0){ # if even
        V[i,j]=(-1)^((i-j)/2)*V[(i+j)/2,(i+j)/2]
        V[j,i]=V[i,j]
      }
    }
  }
  return(V)
}

# CycleCarma<-function(y, statevar, expA, times.obs=integer(),
#                      p=integer(), Qmatr, SigMatr, zc, loglstar){
# #   expAT=t(expA)
# #   zcT=t(zc)
#  # for(t in 1:times.obs){
# #   t=1
# # #     # prediction
# #     statevar <- expA %*% statevar
# #     SigMatr <- expA %*% SigMatr %*% t(expA) + Qmatr
# #     # forecast
# #     Uobs <- y[t] - zc %*% statevar # 1 - 1Xp px1
# #     dum.zc <- zc %*% SigMatr  # 1xp pxp
# #     sd_2 <- dum.zc %*% t(zc)  # 1xp px1
# #     Inv_sd_2 <- 1/sd_2
# #     #correction
# #     Kgain <- SigMatr %*%  t(zc) %*% Inv_sd_2 # pxp px1*1
# #     statevar <- statevar+Kgain %*% Uobs # px1+px1
# #     SigMatr <- SigMatr - Kgain %*% dum.zc # pxp-px1 1x+
# #     term_int<- -0.5 * (log(sd_2)+ Uobs %*% Uobs %*% Inv_sd_2) # every entries are scalars
# #     loglstar <- loglstar + term_int # every entries are scalars
# #   }
# #   expA=matrix(c(1:16),nrow=4,ncol=4)
# #   SigMatr=matrix(c(1:16),nrow=4,ncol=4)+1
# #   Qmatr=matrix(c(1:16),nrow=4,ncol=4)+2
# #   vvvvv<-expA%*%SigMatr
# #   ppppp<-expA%*%SigMatr%*%t(expA)+Qmatr
#   rY=as.numeric(y)
#   rStateVar=as.numeric(statevar)
#   rExpA=as.numeric(expA)
#   rtime_obs=times.obs
#   p=p
#   rQmatr=as.numeric(Qmatr)
#   rSigMatr=as.numeric(SigMatr)
#   rZc=as.numeric(zc)
#   rLogstar=loglstar
#   In_dum=0
#   sd_2=0
#   rMat21=numeric(length=p)
#   rdum_zc=numeric(length=p)
#   rMat22int=numeric(length=p*p)
#   rMat22est=numeric(length=p*p)
#   rKgain=numeric(length=p)
#   for(t in 1:rtime_obs){
#     # prediction
#     for(i in 1:p){
#       rMat21[(i-1)+1] = 0
#       for(j in 1:p){
#         #     statevar <- expA %*% statevar: px1=pxp px1
#         rMat21[(i-1)+1] = rMat21[(i-1)+1]+rExpA[(i-1)+(j-1)*p+1]*rStateVar[(j-1)+1]
#       }
#         rStateVar[(i-1)+1] = rMat21[(i-1)+1]  #     statevar <- expA %*% statevar
#     }
#
# #   SigMatr <- expA %*% SigMatr %*% expAT + Qmatr: pxp = pxp pxp pxp
#       # First We compute rMat22int <- expA %*% SigMatr : pxp = pxp pxp
#     for(i in 1:p){
#       for(j in 1:p){
#         rMat22int[(i-1)+(j-1)*p+1]=0
#         for(h in 1:p){
#             rMat22int[(i-1)+(j-1)*p+1]=rMat22int[(i-1)+(j-1)*p+1]+rExpA[(i-1)+(h-1)*p+1]*
#             rSigMatr[(h-1)+(j-1)*p+1]
#         }
#       }
#     }
#       # Second We compute rMat22est <- rMat22int %*% t(expA) + Qmatr: pxp = pxp pxp + pxp
#     for(i in 1:p){
#       for(j in 1:p){
#         rMat22est[(i-1)+(j-1)*p+1]=0
#         for(h in 1:p){
#           rMat22est[(i-1)+(j-1)*p+1]=rMat22est[(i-1)+(j-1)*p+1]+rMat22int[(i-1)+(h-1)*p+1]*rExpA[(j-1)+(h-1)*p+1]
#
#         }
#         rSigMatr[(i-1)+(j-1)*p+1]=rMat22est[(i-1)+(j-1)*p+1]+rQmatr[(i-1)+(j-1)*p+1]
#       }
#     }
# #     # forecast
#
# #     Uobs <- y[t] - zc %*% statevar # 1 - 1Xp px1
#       rMat22est[1]=0
#     for(i in c(1:p)){
#       rMat22est[1]=rMat22est[1]+rZc[i]*rStateVar[i]
#     }
#     Uobs=rY[t]-rMat22est[1]
#
#  #   dum.zc <- zc %*% SigMatr  # 1xp pxp
#
#
#     for(i in c(1:p)){
#       rdum_zc[i]=0
#       for(j in c(1:p)){
#         rdum_zc[i]=rdum_zc[i]+rZc[j]*rSigMatr[(i-1)*h+j-1+1]
#       }
#     }
# #     sd_2 <- dum.zc %*% zcT  # 1xp px1
#     sd_2=0
#     for(i in c(1:p)){
#         sd_2=sd_2+rdum_zc[i]*rZc[i]
#     }
# #     #correction
# #   Kgain <- SigMatr %*% zcT %*% 1/sd_2 # pxp px1*1
#     for(i in c(1:p)){
#       rMat21[i]=0
#       for(j in c(1:p)){
#         rMat21[i]=rMat21[i]+rSigMatr[(i-1)+(j-1)*p+1]*rZc[j]
#       }
#       rKgain[i]=rMat21[i]/sd_2
#     }
#
#
# #     statevar <- statevar+Kgain %*% Uobs # px1+px1
#      for(i in c(1:p)){
#        rStateVar[i] = rStateVar[i] + rKgain[i]*Uobs
#      }
# #     SigMatr <- SigMatr - Kgain %*% dum.zc # pxp-px1 1xp
#     for(i in c(1:p)){
#       for(j in c(1:p)){
#         rSigMatr[(i-1)+(j-1)*p+1] =rSigMatr[(i-1)+(j-1)*p+1]-rKgain[i]*rdum_zc[j]
#       }
#     }
#
#      term_int = -0.5 * (log(sd_2)+ Uobs * Uobs * 1/sd_2) # every entries are scalars
#      loglstar = loglstar + term_int # every entries are scalars
#
#
#   }
#   Res<-matrix(c(loglstar,sd_2),nrow=2,ncol=1)
#   return(Res)
# }

yuima.PhamBreton.Alg<-function(a){
  p<-length(a)
  gamma<-a[p:1]
  if(p>2){
    gamma[p]<-a[1]
    alpha<-matrix(NA,p,p)
    for(j in 1:p){
      if(is.integer(as.integer(j)/2)){
        alpha[p,j]<-0
        alpha[p-1,j]<-0
      }else{
        alpha[p,j]<-a[j]
        alpha[p-1,j]<-a[j+1]/gamma[p]
      }
    }
    for(n in (p-1):1){
      gamma[n]<-alpha[n+1,2]-alpha[n,2]
      for(j in 1:n-1){
        alpha[n-1,j]<-(alpha[n+1,j+2]-alpha[n,j+2])/gamma[n]
      }
      alpha[n-1,n-1]<-alpha[n+1,n+1]/gamma[n]
    }
    gamma[1]<-alpha[2,2]
  }
  return(gamma)
}

#yuima.PhamBreton.Inv<-function(gamma){
#   p<-length(gamma)
#   a<-gamma[p:1]
#   if(p>2){
#     x<-polynom()
#     f0<-1*x^0
#     f1<-x
#     f2<-x*f1+gamma[1]*f0
#     for(t in 2:(p-1)){
#       f0<-f1
#       f1<-f2
#       f2<-x*f1+gamma[t]*f0
#     }
#     finpol<-f2+gamma[p]*f1
#     a <- coef(finpol)[p:1]
#   }
#   return(a)
# }



#yuima.carma.loglik1<-function (y, tt, a, b, sigma)
yuima.carma.loglik1<-function (y, u, a, b, sigma,time.obs,V_inf0,p,q)
{
  #This code compute the LogLik using kalman filter

  # if(a_0!=0){we need to correct the Y_t for the mean}
  # if(sigma!=1){we need to write}
  #p <- as.integer(length(a))

#  p <- length(a)

#  bvector <- rep(0, p)
#  q <- length(b)
  bvector <- c(b, rep(0, p - q-1))


  sigma<-sigma
  y<-y

  #xxalt<-carma.kalman(y, tt, p, q, a,bvector,sigma)

  xxalt<-carma.kalman(y, u, p, q, a,bvector,sigma,time.obs,V_inf0)
  list(loglikCdiag = xxalt$loglstar,s2hat=xxalt$s2hat)
}

# returns the vector of log-transitions instead of the final quasilog
quasiloglvec <- function(yuima, param, print=FALSE, env){

	diff.par <- yuima@model@parameter@diffusion
	drift.par <- yuima@model@parameter@drift

	fullcoef <- NULL

	if(length(diff.par)>0)
	fullcoef <- diff.par

	if(length(drift.par)>0)
	fullcoef <- c(fullcoef, drift.par)

	npar <- length(fullcoef)

	nm <- names(param)
    oo <- match(nm, fullcoef)

    if(any(is.na(oo)))
        yuima.stop("some named arguments in 'param' are not arguments to the supplied yuima model")
    param <- param[order(oo)]
    nm <- names(param)

	idx.diff <- match(diff.par, nm)
	idx.drift <- match(drift.par, nm)

	h <- env$h

    theta1 <- unlist(param[idx.diff])
    theta2 <- unlist(param[idx.drift])
	n.theta1 <- length(theta1)
	n.theta2 <- length(theta2)
	n.theta <- n.theta1+n.theta2

	d.size <- yuima@model@equation.number
	n <- length(yuima)[1]


	drift <- drift.term(yuima, param, env)
	diff <- diffusion.term(yuima, param, env)

	QL <- numeric(n-1)  ## here is the difference

	pn <- 0


	vec <- env$deltaX-h*drift[-n,]

	K <- -0.5*d.size * log( (2*pi*h) )

	dimB <- dim(diff[, , 1])

	if(is.null(dimB)){  # one dimensional X
        for(t in 1:(n-1)){
            yB <- diff[, , t]^2
            logdet <- log(yB)
            pn <- K - 0.5*logdet-0.5*vec[t, ]^2/(h*yB)
            QL[t] <- pn

		}
	} else {  # multidimensional X
        for(t in 1:(n-1)){
            yB <- diff[, , t] %*% t(diff[, , t])
            logdet <- log(det(yB))
            if(is.infinite(logdet) ){ # should we return 1e10?
                pn <- log(1)
                yuima.warn("singular diffusion matrix")
                return(1e10)
            }else{
                pn <- K - 0.5*logdet +
                ((-1/(2*h))*t(vec[t, ])%*%solve(yB)%*%vec[t, ])
                QL[t] <- pn
            }
        }
	}
	return(QL)
}




setMethod("summary", "yuima.qmle",
function (object, ...)
{
    cmat <- cbind(Estimate = object@coef, `Std. Error` = sqrt(diag(object@vcov)))
    m2logL <- 2 * object@min

    tmp <- new("summary.yuima.qmle", call = object@call, coef = cmat,
   	m2logL = m2logL,
	model = object@model
    )
    tmp
}
)


setMethod("show", "summary.yuima.qmle",
function (object)
{

    cat("Quasi-Maximum likelihood estimation\n\nCall:\n")
    print(object@call)
    cat("\nCoefficients:\n")
    print(coef(object))
    cat("\n-2 log L:", object@m2logL, "\n")
  }
)

setMethod("plot",signature(x="yuima.CP.qmle"),
function(x, ...){
    t <- x@Jump.times
    X <- x@X.values
    points(x=t,y=X, ...)
}
)


setMethod("summary", "yuima.CP.qmle",
function (object, ...)
{
    cmat <- cbind(Estimate = object@coef, `Std. Error` = sqrt(diag(object@vcov)))
    m2logL <- 2 * object@min
    x <- object@X.values
    j <- object@Jump.values
    t <- object@Jump.times

    tmp <- new("summary.yuima.CP.qmle", call = object@call, coef = cmat,
   	m2logL = m2logL, NJ = length(t),
  	MeanJ = mean(j,na.rm=TRUE),
	SdJ = sd(j,na.rm=TRUE),
	MeanT = mean(diff(t),na.rm=TRUE),
	X.values = x,
	Jump.values = j,
	Jump.times = t,
	model = object@model,
    threshold=object@threshold
    )
    tmp
}
)



setMethod("show", "summary.yuima.CP.qmle",
function (object)
{

    cat("Quasi-Maximum likelihood estimation\n\nCall:\n")
    print(object@call)
    cat("\nCoefficients:\n")
    print(coef(object))
    cat("\n-2 log L:", object@m2logL, "\n")

    cat(sprintf("\n\nNumber of estimated jumps: %d\n",object@NJ))
    cat(sprintf("\nAverage inter-arrival times: %f\n",object@MeanT))
    cat(sprintf("\nAverage jump size: %f\n",object@MeanJ))
    cat(sprintf("\nStandard Dev. of jump size: %f\n",object@SdJ))
    cat(sprintf("\nJump Threshold: %f\n",object@threshold))
    cat("\nSummary statistics for jump times:\n")
    print(summary(object@Jump.times))
    cat("\nSummary statistics for jump size:\n")
    print(summary(object@Jump.values,na.rm=TRUE))
    cat("\n")
}
)
# Utilities for estimation of levy in continuous arma model

setMethod("summary", "yuima.carma.qmle",
          function (object, ...)
          {
            cmat <- cbind(Estimate = object@coef, `Std. Error` = sqrt(diag(object@vcov)))
            m2logL <- 2 * object@min
            data<-Re(coredata(object@Incr.Lev))
            data<- data[!is.na(data)]

            tmp <- new("summary.yuima.carma.qmle", call = object@call, coef = cmat,
                       m2logL = m2logL,
                       MeanI = mean(data),
                       SdI = sd(data),
                       logLI = object@logL.Incr,
                       TypeI = object@model@measure.type,
                       NumbI = length(data),
                       StatI =summary(data)
            )
            tmp
          }
)

setMethod("show", "summary.yuima.carma.qmle",
          function (object)
          {

            cat("Two Stage Quasi-Maximum likelihood estimation\n\nCall:\n")
            print(object@call)
            cat("\nCoefficients:\n")
            print(coef(object))
            cat("\n-2 log L:", object@m2logL, "\n")

            cat(sprintf("\n\nNumber of increments: %d\n",object@NumbI))
            cat(sprintf("\nAverage of increments: %f\n",object@MeanI))
            cat(sprintf("\nStandard Dev. of increments: %f\n",object@SdI))
            if(!is.null(object@logLI)){
              cat(sprintf("\n\n-2 log L of increments: %f\n",-2*object@logLI))
            }
            cat("\nSummary statistics for increments:\n")
            print(object@StatI)
            cat("\n")
          }
)



  # Plot Method for yuima.carma.qmle
setMethod("plot",signature(x="yuima.carma.qmle"),
          function(x, ...){
            Time<-index(x@Incr.Lev)
            Incr.L<-coredata(x@Incr.Lev)
            if(is.complex(Incr.L)){
              yuima.warn("Complex increments. We plot only the real part")
              Incr.L<-Re(Incr.L)
            }
            plot(x=Time,y=Incr.L, ...)
          }
)





#Density code for compound poisson

#CPN

dCPN<-function(x,lambda,mu,sigma){
  a<-min(mu-100*sigma,min(x)-1)
  b<-max(mu+100*sigma,max(x)+1)
  ChFunToDens.CPN <- function(n, a, b, lambda, mu, sigma) {
    i <- 0:(n-1)            # Indices
    dx <- (b-a)/n           # Step size, for the density
    x <- a + i * dx         # Grid, for the density
    dt <- 2*pi / ( n * dx ) # Step size, frequency space
    c <- -n/2 * dt          # Evaluate the characteristic function on [c,d]
    d <-  n/2 * dt          # (center the interval on zero)
    t <- c + i * dt         # Grid, frequency space
    charact.CPN<-function(t,lambda,mu,sigma){
      normal.y<-exp(1i*t*mu-sigma^2*t^2/2)
      y<-exp(lambda*(normal.y-1))
    }
    phi_t <- charact.CPN(t,lambda,mu,sigma)
    X <- exp( -(0+1i) * i * dt * a ) * phi_t
    Y <- fft(X)
    density <- dt / (2*pi) * exp( - (0+1i) * c * x ) * Y
    data.frame(
      i = i,
      t = t,
      characteristic_function = phi_t,
      x = x,
      density = Re(density)
    )
  }
  invFFT<-ChFunToDens.CPN(lambda=lambda,mu=mu,sigma=sigma,n=2^10,a=a,b=b)
  dens<-approx(invFFT$x,invFFT$density,x)
  return(dens$y)
}

# CExp

dCPExp<-function(x,lambda,rate){
  a<-10^-6
  b<-max(1/rate*10 +1/rate^2*10 ,max(x[!is.na(x)])+1)
  ChFunToDens.CPExp <- function(n, a, b, lambda, rate) {
    i <- 0:(n-1)            # Indices
    dx <- (b-a)/n           # Step size, for the density
    x <- a + i * dx         # Grid, for the density
    dt <- 2*pi / ( n * dx ) # Step size, frequency space
    c <- -n/2 * dt          # Evaluate the characteristic function on [c,d]
    d <-  n/2 * dt          # (center the interval on zero)
    t <- c + i * dt         # Grid, frequency space
    charact.CPExp<-function(t,lambda,rate){
      normal.y<-(rate/(1-1i*t))
      # exp(1i*t*mu-sigma^2*t^2/2)
      y<-exp(lambda*(normal.y-1))
    }
    phi_t <- charact.CPExp(t,lambda,rate)
    X <- exp( -(0+1i) * i * dt * a ) * phi_t
    Y <- fft(X)
    density <- dt / (2*pi) * exp( - (0+1i) * c * x ) * Y
    data.frame(
      i = i,
      t = t,
      characteristic_function = phi_t,
      x = x,
      density = Re(density)
    )
  }
  invFFT<-ChFunToDens.CPExp(lambda=lambda,rate=rate,n=2^10,a=a,b=b)
  dens<-approx(invFFT$x[!is.na(invFFT$density)],invFFT$density[!is.na(invFFT$density)],x)
  return(dens$y[!is.na(dens$y)])
}

# CGamma

dCPGam<-function(x,lambda,shape,scale){
  a<-10^-6
  b<-max(shape*scale*10 +shape*scale^2*10 ,max(x[!is.na(x)])+1)
  ChFunToDens.CPGam <- function(n, a, b, lambda, shape,scale) {
    i <- 0:(n-1)            # Indices
    dx <- (b-a)/n           # Step size, for the density
    x <- a + i * dx         # Grid, for the density
    dt <- 2*pi / ( n * dx ) # Step size, frequency space
    c <- -n/2 * dt          # Evaluate the characteristic function on [c,d]
    d <-  n/2 * dt          # (center the interval on zero)
    t <- c + i * dt         # Grid, frequency space
    charact.CPGam<-function(t,lambda,shape,scale){
      normal.y<-(1-1i*t*scale)^(-shape)
      # exp(1i*t*mu-sigma^2*t^2/2)
      y<-exp(lambda*(normal.y-1))
    }
    phi_t <- charact.CPGam(t,lambda,shape,scale)
    X <- exp( -(0+1i) * i * dt * a ) * phi_t
    Y <- fft(X)
    density <- dt / (2*pi) * exp( - (0+1i) * c * x ) * Y
    data.frame(
      i = i,
      t = t,
      characteristic_function = phi_t,
      x = x,
      density = Re(density)
    )
  }
  invFFT<-ChFunToDens.CPGam(lambda=lambda,shape=shape,scale=scale,n=2^10,a=a,b=b)
  dens<-approx(invFFT$x[!is.na(invFFT$density)],invFFT$density[!is.na(invFFT$density)],x)
  return(dens$y[!is.na(dens$y)])
}


minusloglik.Lev<-function(par,env){
  if(env$measure.type=="code"){
    if(env$measure=="rNIG"){
      alpha<-par[1]
      beta<-par[2]
      delta<-par[3]
      mu<-par[4]
      f<-dNIG(env$data,alpha,beta,delta,mu)
      v<-log(as.numeric(na.omit(f)))
      v1<-v[!is.infinite(v)]
      -sum(v1)
    }else{
      if(env$measure=="rngamma"){
        lambda<-par[1]
        alpha<-par[2]
        beta<-par[3]
        mu<-par[4]
        f<-dngamma(env$data,lambda,alpha,beta,mu)
        v<-log(as.numeric(na.omit(f)))
        v1<-v[!is.infinite(v)]
        -sum(v1)
      }else{
        if(env$measure=="rIG"){
          delta<-par[1]
          gamma<-par[2]
          f<-dIG(env$data,delta,gamma)
          v<-log(as.numeric(na.omit(f)))
          v1<-v[!is.infinite(v)]
          -sum(v1)
        }
      }
    }
  }else{
    if(env$measure=="dnorm"){
      lambda<-par[1]
      mu<-par[2]
      sigma<-par[3]
      f<-dCPN(env$data,lambda,mu,sigma)
      v<-log(as.numeric(na.omit(f)))
      v1<-v[!is.infinite(v)]
      -sum(v1)
    }else{
      if(env$measure=="dexp"){
        lambda<-par[1]
        rate<-par[2]
    #    -sum(log(dCPExp(env$data,lambda,rate)))

        f<-dCPExp(env$data,lambda,rate)
        v<-log(as.numeric(na.omit(f)))
        v1<-v[!is.infinite(v)]
        -sum(v1)

      }else{
        if(env$measure=="dgamma"){
          lambda<-par[1]
          shape<-par[2]
          scale<-par[3]
#          -sum(log(dCPGam(env$data,lambda,shape,scale)))

          f<-dCPGam(env$data,lambda,shape,scale)
          v<-log(as.numeric(na.omit(f)))
          v1<-v[!is.infinite(v)]
          -sum(v1)


        }
      }
    }
  }
}



Lev.hessian<-function (params,env){
  logLik.Lev <- function(params){
    if(env$measure.type=="code"){
      if(env$measure=="rNIG"){
        alpha<-params[1]
        beta<-params[2]
        delta<-params[3]
        mu<-params[4]
      #  return(sum(log(dNIG(env$data,alpha,beta,delta,mu))))
        f<-dNIG(env$data,alpha,beta,delta,mu)
        v<-log(as.numeric(na.omit(f)))
        v1<-v[!is.infinite(v)]
        return(sum(v1))
      }else{
        if(env$measure=="rngamma"){
          lambda<-params[1]
          alpha<-params[2]
          beta<-params[3]
          mu<-params[4]
          #return(sum(log(dngamma(env$data,lambda,alpha,beta,mu))))
          f<-dngamma(env$data,lambda,alpha,beta,mu)
          v<-log(as.numeric(na.omit(f)))
          v1<-v[!is.infinite(v)]
          return(sum(v1))
        }else{
          if(env$measure=="rIG"){
            delta<-params[1]
            gamma<-params[2]
            f<-dIG(env$data,delta,gamma)
            v<-log(as.numeric(na.omit(f)))
            v1<-v[!is.infinite(v)]
            return(sum(v1))
          }else{
            if(env$measure=="rgamma"){

               shape<-params[1]
               rate<-params[2]
               f<-dgamma(env$data,shape,rate)
               v<-log(as.numeric(na.omit(f)))
               v1<-v[!is.infinite(v)]
               return(sum(v1))
            }
          }
        }
      }
    }else{
      if(env$measure=="dnorm"){
        lambda<-params[1]
        mu<-params[2]
        sigma<-params[3]
        return(sum(log(dCPN(env$data,lambda,mu,sigma))))
        }else{
          if(env$measure=="dexp"){
            lambda<-params[1]
            rate<-params[2]
            return(sum(log(dCPExp(env$data,lambda,rate))))
            }else{
              if(env$measure=="dgamma"){
                lambda<-params[1]
                shape<-params[2]
                scale<-params[3]
                return(sum(log(dCPGam(env$data,lambda,shape,scale))))
              }
            }
        }
    }
  }
  hessian<-tryCatch(optimHess(par=params, fn=logLik.Lev),
                    error=function(theta){matrix(NA,env$lengpar,env$lengpar)})
  if(env$aggregation==FALSE){
    if(env$measure.type=="CP"){
      Matr.dum<-diag(c(1/env$dt, rep(1, (length(params)-1))))
    }else{
      if(env$measure=="rNIG"){
        Matr.dum<-diag(c(1,1,1/env$dt,1/env$dt))
      }else{
        if(env$measure=="rngamma"){
          Matr.dum<-diag(c(1/env$dt,1,1,1/env$dt))
        }else{
          if(env$measure=="rIG"){
            Matr.dum<-diag(c(1/env$dt,1))
          }else{
            if(env$measure=="rgamma"){
              Matr.dum<-diag(c(1/env$dt,1))
            }
          }
        }
      }
    }
    cov<--Matr.dum%*%solve(hessian)%*%Matr.dum
  }else{
    cov<--solve(hessian)
  }
  return(cov)
}



yuima.Estimation.Lev<-function(Increment.lev,param0,
                               fixed.carma=fixed.carma,
                               lower.carma=lower.carma,
                               upper.carma=upper.carma,
                               measure=measure,
                               measure.type=measure.type,
                               dt=env$h,
                               aggregation=aggregation){


  env<-new.env()
  env$data<-Increment.lev
  env$measure<-measure
  env$measure.type<-measure.type
  # Only one problem
  env$dt<-dt


  if(aggregation==FALSE){
    if(measure.type=="code"){
      if(env$measure=="rNIG"){
        #Matr.dum<-diag(c(1,1,1/env$dt,1/env$dt))
        param0[3]<-param0[3]*dt
        param0[4]<-param0[4]*dt
      }else{
        if(env$measure=="rngamma"){
          #Matr.dum<-diag(c(1/env$dt,1,1,1/env$dt))
          param0[1]<-param0[1]*dt
          param0[4]<-param0[4]*dt
        }else{
          if(env$measure=="rIG"){
            #Matr.dum<-diag(c(1/env$dt,1))
            param0[1]<-param0[1]*dt
          }else{
            if(env$measure=="rgamma"){
              param0[1]<-param0[1]*dt
            }
          }
        }
      }
    }else{
      param0[1]<-param0[1]*dt
    }
  }




  # For NIG
  if(measure.type=="code"){
    if(measure=="rNIG"){
      ui<-rbind(c(1, -1, 0, 0),c(1, 1, 0, 0),c(1, 0, 0, 0),c(0, 0, 1, 0))
      ci<-c(0,0,0,10^(-6))
    }else{
      if(measure=="rngamma"){
        ui<-rbind(c(1,0, 0, 0),c(0, 1, 1, 0),c(0, 1,-1, 0),c(0, 1,0, 0))
        ci<-c(10^-6,10^-6,10^(-6), 0)
      }else{
        if(measure=="rIG"){
          ui<-rbind(c(1,0),c(0, 1))
          ci<-c(10^-6,10^-6)
        }else{
          if(measure=="rgamma"){
            ui<-rbind(c(1,0),c(0, 1))
            ci<-c(10^-6,10^-6)
          }
        }
      }
    }
  }else{
    if(measure=="dnorm"){
      ui<-rbind(c(1,0,0),c(0,0,1))
      ci<-c(10^-6,10^-6)
    }else{
      if(measure=="dexp"){
        ui<-rbind(c(1,0),c(0,1))
        ci<-c(10^-6,10^-6)
      }else{
        if(measure=="dgamma"){
          ui<-rbind(c(1,0,0),c(0,1,0),c(0,0,1))
          ci<-c(10^-6,10^-6,10^-6)
        }
      }
    }
  }



  if(!is.null(lower.carma)){
    lower.con<-matrix(0,length(lower.carma),length(param0))
    rownames(lower.con)<-names(lower.carma)
    colnames(lower.con)<-names(param0)
    numb.lower<-length(lower.carma)
    lower.con[names(lower.carma),names(lower.carma)]<-1*diag(numb.lower)
    dummy.lower.names<-paste0(names(lower.carma),".lower")
    rownames(lower.con)<-dummy.lower.names
    names(lower.carma)<-dummy.lower.names
    ui<-rbind(ui,lower.con)
    ci<-c(ci,lower.carma)
    #idx.lower.carma<-match(names(lower.carma),names(param0))
  }
  if(!is.null(upper.carma)){
    upper.con<-matrix(0,length(upper.carma),length(param0))
    rownames(upper.con)<-names(upper.carma)
    colnames(upper.con)<-names(param0)
    numb.upper<-length(upper.carma)
    upper.con[names(upper.carma),names(upper.carma)]<--1*diag(numb.upper)
    dummy.upper.names<-paste0(names(upper.carma),".upper")
    rownames(upper.con)<-dummy.upper.names
    names(upper.carma)<-dummy.upper.names
    ui<-rbind(ui,upper.con)
    ci<-c(ci,-upper.carma)
  }
  if(!is.null(fixed.carma)){
    names.fixed<-names(fixed.carma)
    numb.fixed<-length(fixed.carma)
    fixed.con<-matrix(0,length(fixed.carma),length(param0))
    rownames(fixed.con)<-names(fixed.carma)
    colnames(fixed.con)<-names(param0)
    fixed.con.bis<-fixed.con
    fixed.con[names(fixed.carma),names(fixed.carma)]<--1*diag(numb.fixed)
    fixed.con.bis[names(fixed.carma),names(fixed.carma)]<-1*diag(numb.fixed)
    dummy.fixed.names<-paste0(names(fixed.carma),".fixed.u")
    dummy.fixed.bis.names<-paste0(names(fixed.carma),".fixed.l")
    rownames(fixed.con)<-dummy.fixed.names
    rownames(fixed.con.bis)<-dummy.fixed.bis.names
    names(fixed.carma)<-dummy.fixed.names
    ui<-rbind(ui,fixed.con,fixed.con.bis)
    ci<-c(ci,-fixed.carma-10^-6,fixed.carma-10^-6)
    #ci<-c(ci,-fixed.carma,fixed.carma)
  }

  lengpar<-length(param0)
  paramLev<-NA*c(1:length(lengpar))

  env$lengpar<-lengpar
  firs.prob<-tryCatch(constrOptim(theta=param0,
                                  f=minusloglik.Lev,grad=NULL,ui=ui,ci=ci,env=env),
                      error=function(theta){NULL})


  if(!is.null(firs.prob)){
    paramLev<-firs.prob$par
    names(paramLev)<-names(param0)
    if(!is.null(fixed.carma)){
      paramLev[names.fixed]<-fixed.carma
      names(paramLev)<-names(param0)
    }
  }else{warning("the start value for levy measure is outside of the admissible region")}

  env$aggregation<-aggregation
  if(is.na(paramLev[1])){
    covLev<-matrix(NA,length(paramLev),length(paramLev))
  }else{
    covLev<-Lev.hessian(params=paramLev,env)
    rownames(covLev)<-names(paramLev)
    if(!is.null(fixed.carma)){
      covLev[names.fixed,]<-matrix(NA,numb.fixed,lengpar)
    }
    colnames(covLev)<-names(paramLev)
    if(!is.null(fixed.carma)){
      covLev[,names.fixed]<-matrix(NA,lengpar,numb.fixed)
    }
  }
  if(aggregation==FALSE){
    if(measure.type=="code"){
      if(env$measure=="rNIG"){
        #Matr.dum<-diag(c(1,1,1/env$dt,1/env$dt))
        paramLev[3]<-paramLev[3]/dt
        paramLev[4]<-paramLev[4]/dt
        }else{
          if(env$measure=="rngamma"){
            #Matr.dum<-diag(c(1/env$dt,1,1,1/env$dt))
            paramLev[1]<-paramLev[1]/dt
            paramLev[4]<-paramLev[4]/dt
            }else{
              if(env$measure=="rIG"){
                #Matr.dum<-diag(c(1/env$dt,1))
                paramLev[1]<-paramLev[1]/dt
                }else{
                  if(env$measure=="rgamma"){
                    paramLev[1]<-paramLev[1]/dt
                  }
                }
            }
        }
      }else{
        paramLev[1]<-paramLev[1]/dt
      }
  }
  results<-list(estLevpar=paramLev,covLev=covLev, value=firs.prob$value)
  return(results)
}





# Normal Inverse Gaussian

# yuima.Estimation.NIG<-function(Increment.lev,param0,
#                                fixed.carma=fixed.carma,
#                                lower.carma=lower.carma,
#                                upper.carma=upper.carma){
#
#   minusloglik.dNIG<-function(par,data){
#     alpha<-par[1]
#     beta<-par[2]
#     delta<-par[3]
#     mu<-par[4]
#     -sum(log(dNIG(data,alpha,beta,delta,mu)))
#   }
#
#   data<-Increment.lev
#
#   # Only one problem
#
#
#   ui<-rbind(c(1, -1, 0, 0),c(1, 1, 0, 0),c(1, 0, 0, 0),c(0, 0, 1, 0))
#   ci<-c(0,0,0,10^(-6))
#
#   if(!is.null(lower.carma)){
#     lower.con<-matrix(0,length(lower.carma),length(param0))
#     rownames(lower.con)<-names(lower.carma)
#     colnames(lower.con)<-names(param0)
#     numb.lower<-length(lower.carma)
#     lower.con[names(lower.carma),names(lower.carma)]<-1*diag(numb.lower)
#     dummy.lower.names<-paste0(names(lower.carma),".lower")
#     rownames(lower.con)<-dummy.lower.names
#     names(lower.carma)<-dummy.lower.names
#     ui<-rbind(ui,lower.con)
#     ci<-c(ci,lower.carma)
#     #idx.lower.carma<-match(names(lower.carma),names(param0))
#   }
#   if(!is.null(upper.carma)){
#     upper.con<-matrix(0,length(upper.carma),length(param0))
#     rownames(upper.con)<-names(upper.carma)
#     colnames(upper.con)<-names(param0)
#     numb.upper<-length(upper.carma)
#     upper.con[names(upper.carma),names(upper.carma)]<--1*diag(numb.upper)
#     dummy.upper.names<-paste0(names(upper.carma),".upper")
#     rownames(upper.con)<-dummy.upper.names
#     names(upper.carma)<-dummy.upper.names
#     ui<-rbind(ui,upper.con)
#     ci<-c(ci,-upper.carma)
#   }
#   if(!is.null(fixed.carma)){
#     names.fixed<-names(fixed.carma)
#     numb.fixed<-length(fixed.carma)
#     fixed.con<-matrix(0,length(fixed.carma),length(param0))
#     rownames(fixed.con)<-names(fixed.carma)
#     colnames(fixed.con)<-names(param0)
#     fixed.con.bis<-fixed.con
#     fixed.con[names(fixed.carma),names(fixed.carma)]<--1*diag(numb.fixed)
#     fixed.con.bis[names(fixed.carma),names(fixed.carma)]<-1*diag(numb.fixed)
#     dummy.fixed.names<-paste0(names(fixed.carma),".fixed.u")
#     dummy.fixed.bis.names<-paste0(names(fixed.carma),".fixed.l")
#     rownames(fixed.con)<-dummy.fixed.names
#     rownames(fixed.con.bis)<-dummy.fixed.bis.names
#     names(fixed.carma)<-dummy.fixed.names
#     ui<-rbind(ui,fixed.con,fixed.con.bis)
#     ci<-c(ci,-fixed.carma-10^-6,fixed.carma-10^-6)
#     #ci<-c(ci,-fixed.carma,fixed.carma)
#   }
#
#
#   firs.prob<-tryCatch(constrOptim(theta=param0,
#                                    f=minusloglik.dNIG,grad=NULL,ui=ui,ci=ci,data=data),
#                        error=function(theta){NULL})
#
#   lengpar<-length(param0)
#   paramLev<-NA*c(1:length(lengpar))
#
#   if(!is.null(firs.prob)){
#     paramLev<-firs.prob$par
#     names(paramLev)<-names(param0)
#     if(!is.null(fixed.carma)){
#       paramLev[names.fixed]<-fixed.carma
#       names(paramLev)<-names(param0)
#     }
#   }else{warning("the start value for levy measure is outside of the admissible region")}
#
#   NIG.hessian<-function (data,params){
#     logLik.NIG <- function(params) {
#
#       alpha<-params[1]
#       beta<-params[2]
#       delta<-params[3]
#       mu<-params[4]
#
#       return(sum(log(dNIG(data,alpha,beta,delta,mu))))
#     }
#     hessian<-optimHess(par=params, fn=logLik.NIG)
#     cov<--solve(hessian)
#     return(cov)
#   }
#
#   if(is.na(paramLev)){
#     covLev<-matrix(NA,length(paramLev),length(paramLev))
#   }else{
#     covLev<-NIG.hessian(data=as.numeric(data),params=paramLev)
#     rownames(covLev)<-names(paramLev)
#     if(!is.null(fixed.carma)){
#       covLev[names.fixed,]<-matrix(NA,numb.fixed,lengpar)
#     }
#     colnames(covLev)<-names(paramLev)
#     if(!is.null(fixed.carma)){
#       covLev[,names.fixed]<-matrix(NA,lengpar,numb.fixed)
#     }
#   }
#   results<-list(estLevpar=paramLev,covLev=covLev)
#   return(results)
# }
#
#
#
# # Variance Gaussian
#
# yuima.Estimation.VG<-function(Increment.lev,param0,
#                               fixed.carma=fixed.carma,
#                               lower.carma=lower.carma,
#                               upper.carma=upper.carma){
#
#   minusloglik.dVG<-function(par,data){
#     lambda<-par[1]
#     alpha<-par[2]
#     beta<-par[3]
#     mu<-par[4]
#     -sum(log(dngamma(data,lambda,alpha,beta,mu)))
#   }
#
#   data<-Increment.lev
#
#   ui<-rbind(c(1,0, 0, 0),c(0, 1, 1, 0),c(0, 1,-1, 0),c(0, 1,0, 0))
#   ci<-c(10^-6,10^-6,10^(-6), 0)
#
#   if(!is.null(lower.carma)){
#     lower.con<-matrix(0,length(lower.carma),length(param0))
#     rownames(lower.con)<-names(lower.carma)
#     colnames(lower.con)<-names(param0)
#     numb.lower<-length(lower.carma)
#     lower.con[names(lower.carma),names(lower.carma)]<-1*diag(numb.lower)
#     dummy.lower.names<-paste0(names(lower.carma),".lower")
#     rownames(lower.con)<-dummy.lower.names
#     names(lower.carma)<-dummy.lower.names
#     ui<-rbind(ui,lower.con)
#     ci<-c(ci,lower.carma)
#     #idx.lower.carma<-match(names(lower.carma),names(param0))
#   }
#   if(!is.null(upper.carma)){
#     upper.con<-matrix(0,length(upper.carma),length(param0))
#     rownames(upper.con)<-names(upper.carma)
#     colnames(upper.con)<-names(param0)
#     numb.upper<-length(upper.carma)
#     upper.con[names(upper.carma),names(upper.carma)]<--1*diag(numb.upper)
#     dummy.upper.names<-paste0(names(upper.carma),".upper")
#     rownames(upper.con)<-dummy.upper.names
#     names(upper.carma)<-dummy.upper.names
#     ui<-rbind(ui,upper.con)
#     ci<-c(ci,-upper.carma)
#   }
#   if(!is.null(fixed.carma)){
#     names.fixed<-names(fixed.carma)
#     numb.fixed<-length(fixed.carma)
#     fixed.con<-matrix(0,length(fixed.carma),length(param0))
#     rownames(fixed.con)<-names(fixed.carma)
#     colnames(fixed.con)<-names(param0)
#     fixed.con.bis<-fixed.con
#     fixed.con[names(fixed.carma),names(fixed.carma)]<--1*diag(numb.fixed)
#     fixed.con.bis[names(fixed.carma),names(fixed.carma)]<-1*diag(numb.fixed)
#     dummy.fixed.names<-paste0(names(fixed.carma),".fixed.u")
#     dummy.fixed.bis.names<-paste0(names(fixed.carma),".fixed.l")
#     rownames(fixed.con)<-dummy.fixed.names
#     rownames(fixed.con.bis)<-dummy.fixed.bis.names
#     names(fixed.carma)<-dummy.fixed.names
#     ui<-rbind(ui,fixed.con,fixed.con.bis)
#     ci<-c(ci,-fixed.carma-10^-6,fixed.carma-10^-6)
#     #ci<-c(ci,-fixed.carma,fixed.carma)
#   }
#
#
#   firs.prob<-tryCatch(constrOptim(theta=param0,
#                                   f=minusloglik.dVG,grad=NULL,ui=ui,ci=ci,data=data),
#                       error=function(theta){NULL})
#
#   lengpar<-length(param0)
#   paramLev<-NA*c(1:length(lengpar))
#
#   if(!is.null(firs.prob)){
#     paramLev<-firs.prob$par
#     names(paramLev)<-names(param0)
#     if(!is.null(fixed.carma)){
#       paramLev[names.fixed]<-fixed.carma
#       names(paramLev)<-names(param0)
#     }
#   }
#
#
#   VG.hessian<-function (data,params){
#     logLik.VG <- function(params) {
#
#       lambda<-params[1]
#       alpha<-params[2]
#       beta<-params[3]
#       mu<-params[4]
#
#       return(sum(log(dngamma(data,lambda,alpha,beta,mu))))
#     }
#     # hessian <- tsHessian(param = params, fun = logLik.VG)
#     #hessian<-optimHess(par, fn, gr = NULL,data=data)
#     hessian<-optimHess(par=params, fn=logLik.VG)
#     cov<--solve(hessian)
#     return(cov)
#   }
#
#   if(is.na(paramLev)){
#     covLev<-matrix(NA,length(paramLev),length(paramLev))
#   }else{
#     covLev<-VG.hessian(data=as.numeric(data),params=paramLev)
#     rownames(covLev)<-names(paramLev)
#     if(!is.null(fixed.carma)){
#       covLev[names.fixed,]<-matrix(NA,numb.fixed,lengpar)
#     }
#     colnames(covLev)<-names(paramLev)
#     if(!is.null(fixed.carma)){
#       covLev[,names.fixed]<-matrix(NA,lengpar,numb.fixed)
#     }
#   }
#   results<-list(estLevpar=paramLev,covLev=covLev)
#   return(results)
# }
#
# # Inverse Gaussian
#
# yuima.Estimation.IG<-function(Increment.lev,param0,
#                               fixed.carma=fixed.carma,
#                               lower.carma=lower.carma,
#                               upper.carma=upper.carma){
#
#   minusloglik.dIG<-function(par,data){
#     delta<-par[1]
#     gamma<-par[2]
#     f<-dIG(data,delta,gamma)
#     v<-log(as.numeric(na.omit(f)))
#     v1<-v[!is.infinite(v)]
#     -sum(v1)
#   }
#
#   data<-Increment.lev
#
#   ui<-rbind(c(1,0),c(0, 1))
#   ci<-c(10^-6,10^-6)
#
#   if(!is.null(lower.carma)){
#     lower.con<-matrix(0,length(lower.carma),length(param0))
#     rownames(lower.con)<-names(lower.carma)
#     colnames(lower.con)<-names(param0)
#     numb.lower<-length(lower.carma)
#     lower.con[names(lower.carma),names(lower.carma)]<-1*diag(numb.lower)
#     dummy.lower.names<-paste0(names(lower.carma),".lower")
#     rownames(lower.con)<-dummy.lower.names
#     names(lower.carma)<-dummy.lower.names
#     ui<-rbind(ui,lower.con)
#     ci<-c(ci,lower.carma)
#     #idx.lower.carma<-match(names(lower.carma),names(param0))
#   }
#   if(!is.null(upper.carma)){
#     upper.con<-matrix(0,length(upper.carma),length(param0))
#     rownames(upper.con)<-names(upper.carma)
#     colnames(upper.con)<-names(param0)
#     numb.upper<-length(upper.carma)
#     upper.con[names(upper.carma),names(upper.carma)]<--1*diag(numb.upper)
#     dummy.upper.names<-paste0(names(upper.carma),".upper")
#     rownames(upper.con)<-dummy.upper.names
#     names(upper.carma)<-dummy.upper.names
#     ui<-rbind(ui,upper.con)
#     ci<-c(ci,-upper.carma)
#   }
#   if(!is.null(fixed.carma)){
#     names.fixed<-names(fixed.carma)
#     numb.fixed<-length(fixed.carma)
#     fixed.con<-matrix(0,length(fixed.carma),length(param0))
#     rownames(fixed.con)<-names(fixed.carma)
#     colnames(fixed.con)<-names(param0)
#     fixed.con.bis<-fixed.con
#     fixed.con[names(fixed.carma),names(fixed.carma)]<--1*diag(numb.fixed)
#     fixed.con.bis[names(fixed.carma),names(fixed.carma)]<-1*diag(numb.fixed)
#     dummy.fixed.names<-paste0(names(fixed.carma),".fixed.u")
#     dummy.fixed.bis.names<-paste0(names(fixed.carma),".fixed.l")
#     rownames(fixed.con)<-dummy.fixed.names
#     rownames(fixed.con.bis)<-dummy.fixed.bis.names
#     names(fixed.carma)<-dummy.fixed.names
#     ui<-rbind(ui,fixed.con,fixed.con.bis)
#     ci<-c(ci,-fixed.carma-10^-6,fixed.carma-10^-6)
#     #ci<-c(ci,-fixed.carma,fixed.carma)
#   }
#
#
#   firs.prob<-tryCatch(constrOptim(theta=param0,
#                                   f=minusloglik.dIG,
#                                   grad=NULL,
#                                   ui=ui,
#                                   ci=ci,
#                                   data=data),
#                       error=function(theta){NULL})
#
#   lengpar<-length(param0)
#   paramLev<-NA*c(1:length(lengpar))
#   if(!is.null(firs.prob)){
#     paramLev<-firs.prob$par
#     names(paramLev)<-names(param0)
#     if(!is.null(fixed.carma)){
#       paramLev[names.fixed]<-fixed.carma
#       names(paramLev)<-names(param0)
#     }
#   }
#
#   IG.hessian<-function (data,params){
#     logLik.IG <- function(params) {
#
#       delta<-params[1]
#       gamma<-params[2]
#       f<-dIG(data,delta,gamma)
#       v<-log(as.numeric(na.omit(f)))
#       v1<-v[!is.infinite(v)]
#       return(sum(v1))
#     }
#     # hessian <- tsHessian(param = params, fun = logLik.VG)
#     #hessian<-optimHess(par, fn, gr = NULL,data=data)
#     hessian<-optimHess(par=params, fn=logLik.IG)
#     cov<--solve(hessian)
#     return(cov)
#   }
#
#   if(is.na(paramLev)){
#     covLev<-matrix(NA,length(paramLev),length(paramLev))
#   }else{
#     covLev<-IG.hessian(data=as.numeric(data),params=paramLev)
#     rownames(covLev)<-names(paramLev)
#     if(!is.null(fixed.carma)){
#       covLev[names.fixed,]<-matrix(NA,numb.fixed,lengpar)
#     }
#     colnames(covLev)<-names(paramLev)
#     if(!is.null(fixed.carma)){
#       covLev[,names.fixed]<-matrix(NA,lengpar,numb.fixed)
#     }
#   }
#   results<-list(estLevpar=paramLev,covLev=covLev)
#   return(results)
# }
#
# # Compound Poisson-Normal
#
# yuima.Estimation.CPN<-function(Increment.lev,param0,
#                                fixed.carma=fixed.carma,
#                                lower.carma=lower.carma,
#                                upper.carma=upper.carma){
#   dCPN<-function(x,lambda,mu,sigma){
#     a<-min(mu-100*sigma,min(x)-1)
#     b<-max(mu+100*sigma,max(x)+1)
#     ChFunToDens.CPN <- function(n, a, b, lambda, mu, sigma) {
#       i <- 0:(n-1)            # Indices
#       dx <- (b-a)/n           # Step size, for the density
#       x <- a + i * dx         # Grid, for the density
#       dt <- 2*pi / ( n * dx ) # Step size, frequency space
#       c <- -n/2 * dt          # Evaluate the characteristic function on [c,d]
#       d <-  n/2 * dt          # (center the interval on zero)
#       t <- c + i * dt         # Grid, frequency space
#       charact.CPN<-function(t,lambda,mu,sigma){
#         normal.y<-exp(1i*t*mu-sigma^2*t^2/2)
#         y<-exp(lambda*(normal.y-1))
#       }
#       phi_t <- charact.CPN(t,lambda,mu,sigma)
#       X <- exp( -(0+1i) * i * dt * a ) * phi_t
#       Y <- fft(X)
#       density <- dt / (2*pi) * exp( - (0+1i) * c * x ) * Y
#       data.frame(
#         i = i,
#         t = t,
#         characteristic_function = phi_t,
#         x = x,
#         density = Re(density)
#       )
#     }
#     invFFT<-ChFunToDens.CPN(lambda=lambda,mu=mu,sigma=sigma,n=2^12,a=a,b=b)
#     dens<-approx(invFFT$x,invFFT$density,x)
#     return(dens$y)
#   }
#
#   minusloglik.dCPN<-function(par,data){
#     lambda<-par[1]
#     mu<-par[2]
#     sigma<-par[3]
#     -sum(log(dCPN(data,lambda,mu,sigma)))
#   }
#
#   data<-Increment.lev
#
#   ui<-rbind(c(1,0,0),c(0,0,1))
#   ci<-c(10^-6,10^-6)
#   if(!is.null(lower.carma)){
#     lower.con<-matrix(0,length(lower.carma),length(param0))
#     rownames(lower.con)<-names(lower.carma)
#     colnames(lower.con)<-names(param0)
#     numb.lower<-length(lower.carma)
#     lower.con[names(lower.carma),names(lower.carma)]<-1*diag(numb.lower)
#     dummy.lower.names<-paste0(names(lower.carma),".lower")
#     rownames(lower.con)<-dummy.lower.names
#     names(lower.carma)<-dummy.lower.names
#     ui<-rbind(ui,lower.con)
#     ci<-c(ci,lower.carma)
#     #idx.lower.carma<-match(names(lower.carma),names(param0))
#   }
#   if(!is.null(upper.carma)){
#     upper.con<-matrix(0,length(upper.carma),length(param0))
#     rownames(upper.con)<-names(upper.carma)
#     colnames(upper.con)<-names(param0)
#     numb.upper<-length(upper.carma)
#     upper.con[names(upper.carma),names(upper.carma)]<--1*diag(numb.upper)
#     dummy.upper.names<-paste0(names(upper.carma),".upper")
#     rownames(upper.con)<-dummy.upper.names
#     names(upper.carma)<-dummy.upper.names
#     ui<-rbind(ui,upper.con)
#     ci<-c(ci,-upper.carma)
#   }
#   if(!is.null(fixed.carma)){
#     names.fixed<-names(fixed.carma)
#     numb.fixed<-length(fixed.carma)
#     fixed.con<-matrix(0,length(fixed.carma),length(param0))
#     rownames(fixed.con)<-names(fixed.carma)
#     colnames(fixed.con)<-names(param0)
#     fixed.con.bis<-fixed.con
#     fixed.con[names(fixed.carma),names(fixed.carma)]<--1*diag(numb.fixed)
#     fixed.con.bis[names(fixed.carma),names(fixed.carma)]<-1*diag(numb.fixed)
#     dummy.fixed.names<-paste0(names(fixed.carma),".fixed.u")
#     dummy.fixed.bis.names<-paste0(names(fixed.carma),".fixed.l")
#     rownames(fixed.con)<-dummy.fixed.names
#     rownames(fixed.con.bis)<-dummy.fixed.bis.names
#     names(fixed.carma)<-dummy.fixed.names
#     ui<-rbind(ui,fixed.con,fixed.con.bis)
#     ci<-c(ci,-fixed.carma-10^-6,fixed.carma-10^-6)
#     #ci<-c(ci,-fixed.carma,fixed.carma)
#   }
#   firs.prob<-tryCatch(constrOptim(theta=param0,
#                                   f=minusloglik.dCPN,
#                                   grad=NULL,
#                                   ui=ui,
#                                   ci=ci,
#                                   data=data),
#                       error=function(theta){NULL})
#
#   lengpar<-length(param0)
#   paramLev<-NA*c(1:lengpar)
#   if(!is.null(firs.prob)){
#     paramLev<-firs.prob$par
#     names(paramLev)<-names(param0)
#     if(!is.null(fixed.carma)){
#       paramLev[names.fixed]<-fixed.carma
#       names(paramLev)<-names(param0)
#     }
#   }
#
#   CPN.hessian<-function (data,params,lengpar){
#     logLik.CPN <- function(params) {
#
#       lambda<-params[1]
#       mu<-params[2]
#       sigma<-params[3]
#       return(sum(log(dCPN(data,lambda,mu,sigma))))
#     }
#     # hessian <- tsHessian(param = params, fun = logLik.VG)
#     #hessian<-optimHess(par, fn, gr = NULL,data=data)
#     hessian<-tryCatch(optimHess(par=params, fn=logLik.CPN),
#                       error=function(theta){matrix(NA,lengpar,lengpar)})
#     cov<--solve(hessian)
#     return(cov)
#   }
#
#   if(is.na(paramLev)){
#     covLev<-matrix(NA, lengpar,lengpar)
#   }else{
#     covLev<-CPN.hessian(data=as.numeric(data),params=paramLev,lengpar=lengpar)
#     rownames(covLev)<-names(paramLev)
#     if(!is.null(fixed.carma)){
#       covLev[names.fixed,]<-matrix(NA,numb.fixed,lengpar)
#     }
#     colnames(covLev)<-names(paramLev)
#     if(!is.null(fixed.carma)){
#       covLev[,names.fixed]<-matrix(NA,lengpar,numb.fixed)
#     }
#   }
#   results<-list(estLevpar=paramLev,covLev=covLev)
#   return(results)
# }
#
# yuima.Estimation.CPExp<-function(Increment.lev,param0,
#                                  fixed.carma=fixed.carma,
#                                  lower.carma=lower.carma,
#                                  upper.carma=upper.carma){
#   dCPExp<-function(x,lambda,rate){
#     a<-10^-6
#     b<-max(1/rate*10 +1/rate^2*10 ,max(x[!is.na(x)])+1)
#     ChFunToDens.CPExp <- function(n, a, b, lambda, rate) {
#       i <- 0:(n-1)            # Indices
#       dx <- (b-a)/n           # Step size, for the density
#       x <- a + i * dx         # Grid, for the density
#       dt <- 2*pi / ( n * dx ) # Step size, frequency space
#       c <- -n/2 * dt          # Evaluate the characteristic function on [c,d]
#       d <-  n/2 * dt          # (center the interval on zero)
#       t <- c + i * dt         # Grid, frequency space
#       charact.CPExp<-function(t,lambda,rate){
#         normal.y<-(rate/(1-1i*t))
#         # exp(1i*t*mu-sigma^2*t^2/2)
#         y<-exp(lambda*(normal.y-1))
#       }
#       phi_t <- charact.CPExp(t,lambda,rate)
#       X <- exp( -(0+1i) * i * dt * a ) * phi_t
#       Y <- fft(X)
#       density <- dt / (2*pi) * exp( - (0+1i) * c * x ) * Y
#       data.frame(
#         i = i,
#         t = t,
#         characteristic_function = phi_t,
#         x = x,
#         density = Re(density)
#       )
#     }
#     invFFT<-ChFunToDens.CPExp(lambda=lambda,rate=rate,n=2^12,a=a,b=b)
#     dens<-approx(invFFT$x[!is.na(invFFT$density)],invFFT$density[!is.na(invFFT$density)],x)
#     return(dens$y[!is.na(dens$y)])
#   }
#
#   minusloglik.dCPExp<-function(par,data){
#     lambda<-par[1]
#     rate<-par[2]
#     -sum(log(dCPExp(data,lambda,rate)))
#   }
#
#   data<-Increment.lev
#
#   ui<-rbind(c(1,0),c(0,1))
#   ci<-c(10^-6,10^-6)
#   if(!is.null(lower.carma)){
#     lower.con<-matrix(0,length(lower.carma),length(param0))
#     rownames(lower.con)<-names(lower.carma)
#     colnames(lower.con)<-names(param0)
#     numb.lower<-length(lower.carma)
#     lower.con[names(lower.carma),names(lower.carma)]<-1*diag(numb.lower)
#     dummy.lower.names<-paste0(names(lower.carma),".lower")
#     rownames(lower.con)<-dummy.lower.names
#     names(lower.carma)<-dummy.lower.names
#     ui<-rbind(ui,lower.con)
#     ci<-c(ci,lower.carma)
#     #idx.lower.carma<-match(names(lower.carma),names(param0))
#   }
#   if(!is.null(upper.carma)){
#     upper.con<-matrix(0,length(upper.carma),length(param0))
#     rownames(upper.con)<-names(upper.carma)
#     colnames(upper.con)<-names(param0)
#     numb.upper<-length(upper.carma)
#     upper.con[names(upper.carma),names(upper.carma)]<--1*diag(numb.upper)
#     dummy.upper.names<-paste0(names(upper.carma),".upper")
#     rownames(upper.con)<-dummy.upper.names
#     names(upper.carma)<-dummy.upper.names
#     ui<-rbind(ui,upper.con)
#     ci<-c(ci,-upper.carma)
#   }
#   if(!is.null(fixed.carma)){
#     names.fixed<-names(fixed.carma)
#     numb.fixed<-length(fixed.carma)
#     fixed.con<-matrix(0,length(fixed.carma),length(param0))
#     rownames(fixed.con)<-names(fixed.carma)
#     colnames(fixed.con)<-names(param0)
#     fixed.con.bis<-fixed.con
#     fixed.con[names(fixed.carma),names(fixed.carma)]<--1*diag(numb.fixed)
#     fixed.con.bis[names(fixed.carma),names(fixed.carma)]<-1*diag(numb.fixed)
#     dummy.fixed.names<-paste0(names(fixed.carma),".fixed.u")
#     dummy.fixed.bis.names<-paste0(names(fixed.carma),".fixed.l")
#     rownames(fixed.con)<-dummy.fixed.names
#     rownames(fixed.con.bis)<-dummy.fixed.bis.names
#     names(fixed.carma)<-dummy.fixed.names
#     ui<-rbind(ui,fixed.con,fixed.con.bis)
#     ci<-c(ci,-fixed.carma-10^-6,fixed.carma-10^-6)
#     #ci<-c(ci,-fixed.carma,fixed.carma)
#   }
#
#   firs.prob<-tryCatch(constrOptim(theta=param0,
#                                   f=minusloglik.dCPExp,
#                                   grad=NULL,
#                                   ui=ui,
#                                   ci=ci,
#                                   data=data),
#                       error=function(theta){NULL})
#
#   lengpar<-length(param0)
#   paramLev<-NA*c(1:length(lengpar))
#   if(!is.null(firs.prob)){
#     paramLev<-firs.prob$par
#     names(paramLev)<-names(param0)
#     if(!is.null(fixed.carma)){
#       paramLev[names.fixed]<-fixed.carma
#       names(paramLev)<-names(param0)
#     }
#   }
#
#
#   CPExp.hessian<-function (data,params){
#     logLik.CPExp <- function(params) {
#
#       lambda<-params[1]
#       rate<-params[2]
#
#       return(sum(log(dCPExp(data,lambda,rate))))
#     }
#     # hessian <- tsHessian(param = params, fun = logLik.VG)
#     #hessian<-optimHess(par, fn, gr = NULL,data=data)
#     hessian<-optimHess(par=params, fn=logLik.CPExp)
#     cov<--solve(hessian)
#     return(cov)
#   }
#
#   if(is.na(paramLev)){
#     covLev<-matrix(NA,length(paramLev),length(paramLev))
#   }else{
#     covLev<-CPExp.hessian(data=as.numeric(data),params=paramLev)
#     rownames(covLev)<-names(paramLev)
#     if(!is.null(fixed.carma)){
#       covLev[names.fixed,]<-matrix(NA,numb.fixed,lengpar)
#     }
#     colnames(covLev)<-names(paramLev)
#     if(!is.null(fixed.carma)){
#       covLev[,names.fixed]<-matrix(NA,lengpar,numb.fixed)
#     }
#   }
#   results<-list(estLevpar=paramLev,covLev=covLev)
#   return(results)
# }
#
# yuima.Estimation.CPGam<-function(Increment.lev,param0,
#                                  fixed.carma=fixed.carma,
#                                  lower.carma=lower.carma,
#                                  upper.carma=upper.carma){
#   dCPGam<-function(x,lambda,shape,scale){
#     a<-10^-6
#     b<-max(shape*scale*10 +shape*scale^2*10 ,max(x[!is.na(x)])+1)
#     ChFunToDens.CPGam <- function(n, a, b, lambda, shape,scale) {
#       i <- 0:(n-1)            # Indices
#       dx <- (b-a)/n           # Step size, for the density
#       x <- a + i * dx         # Grid, for the density
#       dt <- 2*pi / ( n * dx ) # Step size, frequency space
#       c <- -n/2 * dt          # Evaluate the characteristic function on [c,d]
#       d <-  n/2 * dt          # (center the interval on zero)
#       t <- c + i * dt         # Grid, frequency space
#       charact.CPGam<-function(t,lambda,shape,scale){
#         normal.y<-(1-1i*t*scale)^(-shape)
#         # exp(1i*t*mu-sigma^2*t^2/2)
#         y<-exp(lambda*(normal.y-1))
#       }
#       phi_t <- charact.CPGam(t,lambda,shape,scale)
#       X <- exp( -(0+1i) * i * dt * a ) * phi_t
#       Y <- fft(X)
#       density <- dt / (2*pi) * exp( - (0+1i) * c * x ) * Y
#       data.frame(
#         i = i,
#         t = t,
#         characteristic_function = phi_t,
#         x = x,
#         density = Re(density)
#       )
#     }
#     invFFT<-ChFunToDens.CPGam(lambda=lambda,shape=shape,scale=scale,n=2^12,a=a,b=b)
#     dens<-approx(invFFT$x[!is.na(invFFT$density)],invFFT$density[!is.na(invFFT$density)],x)
#     return(dens$y[!is.na(dens$y)])
#   }
#
#   minusloglik.dCPGam<-function(par,data){
#     lambda<-par[1]
#     shape<-par[2]
#     scale<-par[3]
#     -sum(log(dCPGam(data,lambda,shape,scale)))
#   }
#
#   data<-Increment.lev
#
#   ui<-rbind(c(1,0,0),c(0,1,0),c(0,1,0))
#   ci<-c(10^-6,10^-6,10^-6)
#
#   if(!is.null(lower.carma)){
#     lower.con<-matrix(0,length(lower.carma),length(param0))
#     rownames(lower.con)<-names(lower.carma)
#     colnames(lower.con)<-names(param0)
#     numb.lower<-length(lower.carma)
#     lower.con[names(lower.carma),names(lower.carma)]<-1*diag(numb.lower)
#     dummy.lower.names<-paste0(names(lower.carma),".lower")
#     rownames(lower.con)<-dummy.lower.names
#     names(lower.carma)<-dummy.lower.names
#     ui<-rbind(ui,lower.con)
#     ci<-c(ci,lower.carma)
#     #idx.lower.carma<-match(names(lower.carma),names(param0))
#   }
#   if(!is.null(upper.carma)){
#     upper.con<-matrix(0,length(upper.carma),length(param0))
#     rownames(upper.con)<-names(upper.carma)
#     colnames(upper.con)<-names(param0)
#     numb.upper<-length(upper.carma)
#     upper.con[names(upper.carma),names(upper.carma)]<--1*diag(numb.upper)
#     dummy.upper.names<-paste0(names(upper.carma),".upper")
#     rownames(upper.con)<-dummy.upper.names
#     names(upper.carma)<-dummy.upper.names
#     ui<-rbind(ui,upper.con)
#     ci<-c(ci,-upper.carma)
#   }
#   if(!is.null(fixed.carma)){
#     names.fixed<-names(fixed.carma)
#     numb.fixed<-length(fixed.carma)
#     fixed.con<-matrix(0,length(fixed.carma),length(param0))
#     rownames(fixed.con)<-names(fixed.carma)
#     colnames(fixed.con)<-names(param0)
#     fixed.con.bis<-fixed.con
#     fixed.con[names(fixed.carma),names(fixed.carma)]<--1*diag(numb.fixed)
#     fixed.con.bis[names(fixed.carma),names(fixed.carma)]<-1*diag(numb.fixed)
#     dummy.fixed.names<-paste0(names(fixed.carma),".fixed.u")
#     dummy.fixed.bis.names<-paste0(names(fixed.carma),".fixed.l")
#     rownames(fixed.con)<-dummy.fixed.names
#     rownames(fixed.con.bis)<-dummy.fixed.bis.names
#     names(fixed.carma)<-dummy.fixed.names
#     ui<-rbind(ui,fixed.con,fixed.con.bis)
#     ci<-c(ci,-fixed.carma-10^-6,fixed.carma-10^-6)
#     #ci<-c(ci,-fixed.carma,fixed.carma)
#   }
#
#
#   firs.prob<-tryCatch(constrOptim(theta=param0,
#                                   f=minusloglik.dCPGam,
#                                   grad=NULL,
#                                   ui=ui,
#                                   ci=ci,
#                                   data=data),
#                       error=function(theta){NULL})
#
#   lengpar<-length(param0)
#   paramLev<-NA*c(1:length(lengpar))
#   if(!is.null(firs.prob)){
#     paramLev<-firs.prob$par
#     names(paramLev)<-names(param0)
#     if(!is.null(fixed.carma)){
#       paramLev[names.fixed]<-fixed.carma
#       names(paramLev)<-names(param0)
#     }
#   }
#
#
#   CPGam.hessian<-function (data,params){
#     logLik.CPGam <- function(params) {
#
#       lambda<-params[1]
#       shape<-params[2]
#       scale<-params[3]
#
#       return(sum(log(dCPGam(data,lambda,shape,scale))))
#     }
#     # hessian <- tsHessian(param = params, fun = logLik.VG)
#     #hessian<-optimHess(par, fn, gr = NULL,data=data)
#     hessian<-optimHess(par=params, fn=logLik.CPGam)
#     cov<--solve(hessian)
#     return(cov)
#   }
#
#   if(is.na(paramLev)){
#     covLev<-matrix(NA,length(paramLev),length(paramLev))
#   }else{
#     covLev<-CPGam.hessian(data=as.numeric(data),params=paramLev)
#     rownames(covLev)<-names(paramLev)
#     if(!is.null(fixed.carma)){
#       covLev[names.fixed,]<-matrix(NA,numb.fixed,lengpar)
#     }
#     colnames(covLev)<-names(paramLev)
#     if(!is.null(fixed.carma)){
#       covLev[,names.fixed]<-matrix(NA,lengpar,numb.fixed)
#     }
#   }
#   results<-list(estLevpar=paramLev,covLev=covLev)
#   return(results)
# }
