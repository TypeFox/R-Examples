######################################################################
#
# Implementation of GLR and ordinary Poisson/NegBin CUSUM
# -- documentation converted to Rd format.
#
# Author: Michael Hoehle (with contributions by Valentin Wimmer)
# Date:   8 Jan 2008
# History
#  - 2016-01-17 added ret="cases" for glr using the NegBin distribution
######################################################################

algo.glrnb <- function(disProgObj,
                       control = list(range=range,c.ARL=5,
                         mu0=NULL, alpha=0, Mtilde=1, M=-1, change="intercept",
                         theta=NULL,dir=c("inc","dec"),
                         ret=c("cases","value"),xMax=1e4)) {

  #Small helper function
  either <- function(cond, whenTrue, whenFalse) { if (cond) return(whenTrue) else return(whenFalse) }

  # Set the default values if not yet set
  if(is.null(control[["c.ARL",exact=TRUE]]))
    control$c.ARL <- 5
  if(is.null(control[["change",exact=TRUE]]))
    control$change <- "intercept"
  if(is.null(control[["Mtilde",exact=TRUE]]))
    control$Mtilde <- 1
  if(is.null(control[["M",exact=TRUE]]))
    control$M <- -1
  if(is.null(control[["dir",exact=TRUE]]))
    control$dir <- "inc"
  if(is.null(control[["ret",exact=TRUE]]))
    control$ret <- "value"
  if(is.null(control[["xMax",exact=TRUE]]))
    control$xMax <- 1e4
  if(!is.null(control[["theta",exact=TRUE]])) {
    if(control[["theta",exact=TRUE]] == 1) {
      stop("Error: theta has to be larger than 1!")
    }
  }
  ##Set alpha to null as default. Not necessary, coz it would be taken from
  ##glrnb output.
  ##if(is.null(control[["alpha",exact=TRUE]])) control$alpha <- 0

  #GLM (only filled if estimated)
  m <- NULL

  ################################################
  #Extract the important parts from the arguments
  ################################################
  observed <- disProgObj$observed
  #range is fixed, but t is modified as we iterate the cusum
  t <- control$range ; range <- control$range
  control$mu0Model <- NULL
  control$dir <- match.arg(control$dir, c("inc","dec"))
  dir <- ifelse(control$dir=="inc",1,-1)
  control$ret <- match.arg(control$ret, c("value","cases"))
  ret <- pmatch(control$ret,c("value","cases"))
  mod <- list()

  # Estimate m (the expected number of cases), i.e. parameter lambda of a
  # poisson distribution based on time points 1:t-1
  if (is.null(control[["mu0",exact=TRUE]]) | is.list(control[["mu0",exact=TRUE]])) {
    #Initialize
    if (is.null(control[["mu0",exact=TRUE]])) control$mu0 <- list()
    if (is.null(control[["mu0",exact=TRUE]][["S"]])) control$mu0$S <- 1
    if (is.null(control[["mu0",exact=TRUE]][["trend"]])) control$mu0$trend <- FALSE
    if (is.null(control[["mu0",exact=TRUE]][["refit"]])) control$mu0$refit <- FALSE
    control$mu0Model <- control$mu0

    #Estimate using a hook function (lazy evaluation)
    control$mu0 <- estimateGLRNbHook()$pred

    mod[[1]] <- estimateGLRNbHook()$mod

    # if it is necessary to estimate alpha. Note: glm.nb uses a different
    # parametrization of the negative binomial distribution, i.e. the
    # variance is 'mu + mu^2/size' (?dnbinom).
    # Hence the correct alpha is 1/theta. But now it's the same every time.
    if(is.null(control[["alpha",exact=TRUE]])) control$alpha <- 1/mod[[1]]$theta
  }

  #The counts
  x <- observed[control$range]
  mu0 <- control$mu0

  #Reserve space for the results
  # start with cusum[timePoint -1] = 0, i.e. set cusum[1] = 0
  alarm <- matrix(data = 0, nrow = length(t), ncol = 1)
  upperbound <- matrix(data = 0, nrow = length(t), ncol = 1)


  #Setup counters for the progress
  doneidx <- 0
  N <- 1
  xm10 <- 0
  noofalarms <- 0
  noOfTimePoints <- length(t)
  #Loop as long as we are not through the sequence
  while (doneidx < noOfTimePoints) {
    # cat("Doneidx === ",doneidx,"\n")
    # Call the C-interface -- this should depend on the type
    if (control$change == "intercept") {

      #Generalized likehood ratio vs. ordinary CUSUM
      if (is.null(control[["theta",exact=TRUE]])) {

        if (control$alpha == 0) { #poisson
          if (control$M > 0 ){ # window limited
          	res <- .C("glr_cusum_window",as.integer(x),as.double(mu0),length(x),as.integer(control$M),as.integer(control$Mtilde),as.double(control$c.ARL),N=as.integer(0),val=as.double(numeric(length(x))),cases=as.double(numeric(length(x))),as.integer(dir),as.integer(ret),PACKAGE="surveillance")
          } else { # standard, not window limited
        	res <- .C("glr_cusum",as.integer(x),as.double(mu0),length(x),as.integer(control$Mtilde),as.double(control$c.ARL),N=as.integer(0),val=as.double(numeric(length(x))),cases=as.double(numeric(length(x))),as.integer(dir),as.integer(ret),PACKAGE="surveillance")
          }
        } else { #negbin. This is direcly the window limited version, does M=-1 work here?
          res <- .C("glr_nb_window",x=as.integer(x),mu0=as.double(mu0),alpha=as.double(control$alpha),lx=length(x),Mtilde=as.integer(control$Mtilde),M=as.integer(control$M),c.ARL=as.double(control$c.ARL),N=as.integer(0),val=as.double(numeric(length(x))),dir=as.integer(dir),PACKAGE="surveillance")
          ##hoehle - 2016-01-17. Try out calculating upper bound in terms of cases
          if (control$ret == "cases") {
            ##Warn that this might be slow.
            message("Return of cases is for the GLR detector based on the negative binomial distribution is currently\n only implemented brute force and hence might be very slow!")

###            browser()
            myx <- x
            res$cases <- rep(0,length(res$val))

            for (pos in seq_len(min(length(x),res$N))) {
              myx <- x
              gotAlarm <- (res$N <= pos) #already got an alarm at the position?
              direction <- ifelse(gotAlarm, -1, 1) #go up or down?
              alarmChange <- FALSE #have we suceeded in changing x such that the alarm status changed?

              #Loop over values until one is such that an alarm at (or before!) the time point is given
              while (!alarmChange & (myx[pos] <= control$xMax) & (myx[pos] >=1)) {
                myx[pos] <- myx[pos] + direction
                ##cat("pos=",pos,"x=",myx[pos],"\n")
                tmpRes <- .C("glr_nb_window",x=as.integer(myx),mu0=as.double(mu0),alpha=as.double(control$alpha),lx=length(myx),Mtilde=as.integer(control$Mtilde),M=as.integer(control$M),c.ARL=as.double(control$c.ARL),N=as.integer(0),val=as.double(numeric(length(myx))),dir=as.integer(dir),PACKAGE="surveillance")
                if (!gotAlarm & (tmpRes$N <= pos)) { alarmChange <- TRUE ; res$cases[pos] <- myx[pos]}
                if (gotAlarm  & (tmpRes$N > pos))  { alarmChange <- TRUE ; res$cases[pos] <- myx[pos] + 1}
              }
              if (!alarmChange) { res$cases[pos] <- ifelse(gotAlarm,NA,1e99) } #didn't find alarm before control$xMax
            }
          }
          ##end new 2016 addition to calculate 'cases' for negbin glrnb

        }
      } else { ###################### !is.null(control$theta), i.e. ordinary CUSUM
        if (control$alpha == 0) { #poisson

          res <- .C("lr_cusum",x=as.integer(x),mu0=as.double(mu0),lx=length(x),as.double(control$theta),c.ARL=as.double(control$c.ARL),N=as.integer(0),val=as.double(numeric(length(x))),cases=as.double(numeric(length(x))),as.integer(ret),PACKAGE="surveillance")

        } else { #negbin
          res <- .C("lr_cusum_nb",x=as.integer(x),mu0=as.double(mu0),alpha=as.double(control$alpha),lx=length(x),as.double(control$theta),c.ARL=as.double(control$c.ARL),N=as.integer(0),val=as.double(numeric(length(x))),cases=as.double(numeric(length(x))),as.integer(ret),PACKAGE="surveillance")
        }
      }
    } else { ################### Epidemic chart #######################
      if (control$change == "epi") {
        if (control$alpha == 0) { #pois
          res <- .C("glr_epi_window",as.integer(x),as.double(mu0),length(x),as.integer(control$Mtilde),as.integer(control$M),as.double(xm10),as.double(control$c.ARL),N=as.integer(0),val=as.double(numeric(length(x))),PACKAGE="surveillance")
        } else {
          res <- .C("glr_nbgeneral_window",as.integer(x),as.double(mu0),alpha=as.double(control$alpha),lx=length(x),Mtilde=as.integer(control$Mtilde),M=as.integer(control$M),xm10=as.double(xm10),c.ARL=as.double(control$c.ARL),N=as.integer(0),val=as.double(numeric(length(x))),dir=as.integer(dir),PACKAGE="surveillance")
        }
      }
    }

    ##In case an alarm found log this and reset the chart at res$N+1
    if (res$N <= length(x)) {
      #Put appropriate value in upperbound
      upperbound[1:res$N + doneidx]  <- either(ret == 1, res$val[1:res$N] ,res$cases[1:res$N])
      alarm[res$N + doneidx] <- TRUE

      #Chop & get ready for next round
      xm10 <- x[res$N] #put start value x_0 to last value
      x <- x[-(1:res$N)] ; t <- t[-(1:res$N)]
      #If no refitting is to be done things are easy
      if (!is.list(control$mu0Model) || (control$mu0Model$refit == FALSE)) {
        mu0 <- mu0[-(1:res$N)]
      } else {
        #Update the range (how to change back??)
        range <- range[-(1:res$N)]
        mu0 <- estimateGLRNbHook()$pred
        mod[[noofalarms+2]] <-  estimateGLRNbHook()$mod
        control$mu0[(doneidx + res$N + 1):length(control$mu0)] <- mu0
        #Note: No updating of alpha is currently done.
      }

      noofalarms <- noofalarms + 1

    }
    doneidx <- doneidx + res$N
  }



  #fix of the problem that no upperbound-statistic is returned after
  #last alarm
  upperbound[(doneidx-res$N+1):nrow(upperbound)] <- either(ret == 1, res$val, res$cases)

  #fix of the problem that no upperbound-statistic is returned
  #in case of no alarm
  if (noofalarms == 0) {
    upperbound <- either(ret==1, res$val, res$cases)
  }

  # ensure upper bound is positive and not NaN
  upperbound[is.na(upperbound)] <- 0
  upperbound[upperbound < 0] <- 0


  # Add name and data name to control object
  algoName <- either(control$alpha == 0, "glrpois:", "glrnb:")
  control$name <- paste(algoName, control$change)
  control$data <- paste(deparse(substitute(disProgObj)))
  control$m    <- m
  control$mu0Model$fitted <- mod

  # return alarm and upperbound vectors
  result <- list(alarm = alarm, upperbound = upperbound,
                 disProgObj=disProgObj,control=control)

  class(result) = "survRes" # for surveillance system result
  return(result)
}


#####################################################################
### Function to estimate a Poisson or glm.nb model on the fly - to be
### called within the algo.glrnb function. Experts can customize this
### function.
#####################################################################

estimateGLRNbHook <- function() {
  #Fetch control object from parent
  control <- parent.frame()$control
  #The period
  p <- parent.frame()$disProgObj$freq
  #Current range to perform surveillance on
  range <- parent.frame()$range

  #Define phase1 & phase2 data set (phase2= the rest)
  train <- 1:(range[1]-1)
  test <- range

  #Perform an estimation based on all observations before timePoint
  #Event better - don't do this at all in the algorithm - force
  #user to do it himself - coz its a model selection problem
  data <- data.frame(y=parent.frame()$disProgObj$observed[train],t=train)
  #Build the model equation
  formula <- "y ~ 1 "
  if (control$mu0Model$trend) { formula <- paste(formula," + t",sep="") }
  for (s in seq_len(control$mu0Model$S)) {
    formula <- paste(formula,"+cos(2*",s,"*pi/p*t)+ sin(2*",s,"*pi/p*t)",sep="")
  }

  ##hoehle - 2016-01-16 -- problematic: a full model was fitted, but
  ##this implied a different alpha. Changed now such that a glm
  ##is fitted having the specified alpha (i.e. theta) fixed.
  ##Determine appropriate fitter function
  if (is.null(control[["alpha",exact=TRUE]])) {
    ##Fit while also estimating alpha (if possible!)
    m <- eval(substitute(glm.nb(form,data=data),list(form=as.formula(formula))))
  } else {
    ##Fit the Poisson GLM
    if (control$alpha == 0) {
      message(paste0("glrnb: Fitting Poisson model because alpha == 0"))
      m <- eval(substitute(glm(form,family=poisson(),data=data),list(form=as.formula(formula))))
    } else {
      message(paste0("glrnb: Fitting glm.nb model with alpha=",control$alpha))
      m <- eval(substitute(glm(form,family=negative.binomial(theta=1/control$alpha),data=data),list(form=as.formula(formula))))
    }
  }

  #Predict mu_{0,t}
  pred <- as.numeric(predict(m,newdata=data.frame(t=range),type="response"))

  return(list(mod=m,pred=pred))
}


######################################################################
# simple wrapper for the Poisson case
######################################################################

algo.glrpois <- function(disProgObj,
                         control = list(range=range,c.ARL=5,
                           mu0=NULL, Mtilde=1, M=-1, change="intercept",
                           theta=NULL,dir=c("inc","dec"),
                           ret=c("cases","value"),xMax=1e4)) {
  
  if (is.null(control$alpha)) {
    control$alpha <- 0
  } else if (control$alpha != 0) {
      stop("algo.glrpois has to operate with control$alpha = 0")
  }
  
  algo.glrnb(disProgObj, control)
  
}
