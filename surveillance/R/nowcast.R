######################################################################
# Function to perform nowcast at a specific day "now" using a procedure
# which takes truncation of the available observations into
# account. The full documentation is available in the nowcast.Rd file.
#
# Author: Michael Hoehle <http://www.math.su.se/~hoehle>
#
# Parameters:
#  now  - a Date object representing today
#  when - a vector of Date objects representing the days to do the forecast for.
#         A requirement is that for all elements in when are smaller or equal
#         than "now".
#  data - the Database containing columns dEventCol and dReportCol, which
#      contain the date of the event and of when the report arrives in
#      the database.
#  dEventCol - name of column in data containing time of event occurence
#  dReportCol - name of column in data containing time of reprt arrival
#  method - which method to use

#  D - maximum delay to consider
#  m - moving window for delay estimation
#  control - a list containing the following arguments
#                * gd.prior.kappa - prior for delay is symmetric Dirichlet
#                                   with concentration parameter gd.prior.kappa
#
# Note: As predictions are done simultaneously the entire vector of observations
#       is casted. Then the subset specified in "when" is returned.
#
# Returns:
#  stsNC object with reporting triangle, delay estimate and prediction interval in the appropriate slots.
#
# Todo:
#  * yt.support to N.tInf support  in nowcast??
#  * bayes.notrunc and bayes.notrunc.bnb could become one code segment
#  * Enable user to provide reporting triangle directly.
#  * Function should work for weekly and monthly data as well
######################################################################

nowcast <- function(now,when,data,dEventCol="dHospital",dReportCol="dReport",
        method=c("bayes.notrunc","bayes.notrunc.bnb","lawless","bayes.trunc","unif","bayes.trunc.ddcp"),
        aggregate.by="1 day",
        D=15, m=NULL,
        control=list(
            dRange=NULL,alpha=0.05,nSamples=1e3,
            N.tInf.prior=c("poisgamma","pois","unif"),
            N.tInf.max=300, gd.prior.kappa=0.1,
            ddcp=list(ddChangepoint=NULL,
                logLambda=c("iidLogGa","tps","rw1","rw2"),
                tau.gamma=1,eta.mu=NULL, eta.prec=NULL,
                mcmc=c(burnin=2500,sample=10000,thin=1)),
            score=FALSE,predPMF=FALSE)) {


  #Check if the runjags package is available (required for bayes.trunc.ddcp to work!
  if ("bayes.trunc.ddcp" %in% method) {
    if (!requireNamespace("runjags",quietly=TRUE)) {
      stop("The \"bayes.trunc.ddcp\" method requires the runjags package to be installed, which is available from CRAN.")
    }
  }

  if ((!inherits(now,"Date")) | (length(now)>1)) {
    stop("The parameter 'now' has to be a single Date.")
  }

  #Check if all when_i<= now
  if (!all(when<=now)) {
    stop("Assertion when<=now failed.")
  }

  #Check that specified methods are all valid
  method <- match.arg(method,c("bayes.notrunc","bayes.notrunc.bnb","lawless","bayes.trunc","unif","bayes.trunc.ddcp"),several.ok=TRUE)

  ######################################################################
  # Time aggregation. Make sure it's a valid aggregational level and
  # move all dates to the "first" of this level.
  # @hoehle: Should work for day, weeks and month. Quarter and year not atm.
  ######################################################################
  aggregate.by <- match.arg(aggregate.by,c("1 day","1 week", "1 month"),several.ok=FALSE)
  epochInPeriodStr <- switch(aggregate.by, "1 day"="1","1 week"="%u", "1 month"="%d")

  if (aggregate.by != "1 day") {
      warning("Moving dates to first of each epoch.")


      #Move dates back to first of each epoch unit
      for (colName in c(dEventCol, dReportCol)) {
          data[,colName] <- data[,colName] - as.numeric(format(data[,colName],epochInPeriodStr)) + 1
      }
      #Check now and when
      if (!all( format( c(now,when),epochInPeriodStr) == 1)) {
          stop("The variables 'now' and 'when' needs to be at the first of each epoch")
      }
  }

  #Choose the corect difference function
  if (aggregate.by == "1 day") {
      timeDelay <- function(d1,d2) {as.numeric(d2-d1)}
  }
  if (aggregate.by == "1 week") {
      timeDelay <- function(d1,d2) { floor(as.numeric(difftime(d2,d1,units="weeks")))  } #Count the number of full weeks
  }
  if (aggregate.by == "1 month") {
      timeDelay <- function(d1,d2) {
           #Helper function from http://stackoverflow.com/questions/1995933/number-of-months-between-two-dates
          monnb <- function(d) {
              lt <- as.POSIXlt(as.Date(d, origin="1900-01-01"))
              lt$year*12 + lt$mon
          }
          monnb(d2) - monnb(d1) #count the number of full months
      }
  }

  ######################################################################
  #If there is a specification of dateRange set dMin and dMax accordingly
  #Otherwise use as limits the range of the data
  ######################################################################
  if (is.null(control[["dRange",exact=TRUE]])) {
    dMin <- min(data[,dEventCol],na.rm=TRUE)
    dMax <- max(data[,dEventCol],na.rm=TRUE)
  } else {
    dMin <- control$dRange[1]
    dMax <- control$dRange[length(control$dRange)]
  }
  #@hoehle - check that dRange is proper
  if (!all( format( c(dMin,dMax), epochInPeriodStr) == 1)) {
      stop("The variables in dRange needs to be at the first of each epoch.")
  }
  dateRange <- seq(dMin,dMax,by=aggregate.by)

  ######################################################################
  # Additional manipulation of the control arguments
  ######################################################################

  #Check if alpha is specified
  if (is.null(control[["alpha",exact=TRUE]])) {
    control$alpha <- 0.05
  }
  if (is.null(control[["N.tInf.prior",exact=TRUE]])) {
    control$N.tInf.prior <- "unif"
  }
  if (is.null(control[["N.tInf.max",exact=TRUE]])) {
    control$N.tInf.max <- 300
  }
  if (is.null(control[["gd.prior.kappa",exact=TRUE]])) {
    control$gd.prior.kappa <- 0.1
  }
  if (is.null(control[["nSamples",exact=TRUE]])) {
    control$nSamples <- 1e3
  }
  if (is.null(control[["score",exact=TRUE]])) {
    control$score <- FALSE
  }

  #Checking for the bayes.trun.ddcp procedure. If so make sure params are set up.
  if ("bayes.trunc.ddcp" %in% method) {
      #If no parameters at all set to defaults.
      if (is.null(control[["ddcp",exact=TRUE]])) {
          control$ddcp <- list(ddChangepoint=NULL,
                               logLambda=c("iidLogGa","tps","rw1","rw2"),
                               tau.gamma=1,
                               mcmc=c(burnin=2500,sample=10000,thin=1))
      }
      #Check form og logLambda
      if (is.null(control[["ddcp",exact=TRUE]][["logLambda",exact=TRUE]])) {
          control[["ddcp"]] <- modifyList(control[["ddcp",exact=TRUE]],list(logLambda="iidLogGa"))
      } else {
          control[["ddcp"]]$logLambda <- match.arg(control[["ddcp"]][["logLambda"]],c("iidLogGa","tps","rw1","rw2"))
      }

      #Check breakpoint to use in case of bayes.trunc.ddcp (delay distribution with breakpoint)
      if (is.null(control[["ddcp",exact=TRUE]][["ddChangepoint",exact=TRUE]]) ||
          (!class(control[["ddcp",exact=TRUE]][["ddChangepoint",exact=TRUE]]) == "Date")) {
          stop("Please specify a Date object as changepoint in control$ddChangepoint.")
      } else {
          if (any(control[["ddcp",exact=TRUE]][["ddChangepoint"]] > now)) {
              warning("Some of the elements in ddChangepoint are beyond 'now'. This might be problematic!")
          }
      }

      #Make this an accessible variable
      ddChangepoint <- control$ddcp$ddChangepoint

      #Precision parameter for gamma coefficients for hazard delay distribution
      if (is.null(control[["ddcp",exact=TRUE]][["tau.gamma",exact=TRUE]])) {
          control[["ddcp"]]$tau.gamma <- 1
      }

      if (is.null(control[["ddcp",exact=TRUE]][["eta.mu",exact=TRUE]])) {
          control[["ddcp"]]$eta.mu <- rep(0,length(ddChangepoint))
      } else {
          if (length(control[["ddcp"]]$eta.mu) != length(ddChangepoint)) {
              stop("length of eta.mu is different from the number of change points in 'ddChangepoint'.")
          }
      }
      if (is.null(control[["ddcp",exact=TRUE]][["eta.prec",exact=TRUE]])) {
          control[["ddcp"]]$eta.prec <- rep(1,length(ddChangepoint))
      } else {
          if (length(control[["ddcp"]]$eta.prec) != length(ddChangepoint)) {
              stop("length of eta.prec is different from the number of change points in 'ddChangepoint'.")
          }
      }


      #Check MCMC options
      if (is.null(control[["ddcp",exact=TRUE]][["mcmc",exact=TRUE]])) {
          control[["ddcp"]][["mcmc"]] <- c(burnin=2500,sample=10000,thin=1)
      } else {
          if (!all(names(control[["ddcp",exact=TRUE]][["mcmc",exact=TRUE]]) %in% c("burnin","sample","thin"))) {
              stop("mcmc options need names 'burnin', 'sample' and 'thin'")
          }
      }
  }

  ######################################################################
  # Do preprocessing of the data
  ######################################################################

  #Create a column containing the reporting delay using the timeDelay
  #function
  data$delay <- timeDelay(data[,dEventCol],data[,dReportCol])

  #Handle delays longer than D.
  #@hoehle - handle that the unit might not just be days
  #notThereButDThere <- (data[,dReportCol] > now) & ((data[,dEventCol]) + D <= now)
  notThereButDThere <- (timeDelay(data[,dReportCol],now) < 0) & (timeDelay(data[,dEventCol],now) >= D)
  if (sum(notThereButDThere,na.rm=TRUE)) {
    warning(paste(sum(notThereButDThere,na.rm=TRUE), " observations > \"now\" due to a delay >D. If delay cut to D they would be there."),sep="")
  }

  #Which observations are available at time s
  #@hoehle: data.sub <- data[ na2FALSE(data[,dReportCol] <= now),]
  data.sub <- data[ na2FALSE(timeDelay(data[,dReportCol],now) >= 0),]
  if (nrow(data.sub)==0) {
    stop(paste("No data available at now=",now,"\n"))
  }

  #Create an sts object containing the observed number of counts until s
  sts <- linelist2sts(data.sub,dEventCol,aggregate.by=aggregate.by,dRange=dateRange)
  sts <- as(sts,"stsNC")

  #Create an extra object containing the "truth" based on data
  sts.truth <- linelist2sts(data,dEventCol,aggregate.by=aggregate.by,dRange=dateRange)

  #List of scores to calculate. Can become an argument later on
  scores <- c("logS","RPS","dist.median","outside.ci")

  #Initialize scoring rule results - to be saved in control slot -- dirty
  SR <- array(0,dim=c(nrow(sts),length(method),length(scores)))

  #List for storing the predictive PMFs.
  if (is.null(control[["predPMF",exact=TRUE]])) {
    control$predPMF <- FALSE
  }
  #Prepare a list of different estimated of the delay CDF
  delayCDF <- list()


  ######################################################################
  # Done manipulating the control list with default arguments
  ######################################################################
  sts@control <- control

  #Save truth
  sts@truth <- sts.truth

  #Reserve space for returning the predictive PMFs
  sts@predPMF <- list()

  ######################################################################
  # Consistency checks
  ######################################################################

  #Check if support of N.tInf is large enough
  if (2*control$N.tInf.max < max(observed(sts),na.rm=TRUE)) {
    warning("N.tInf.max appears too small. Largest observed value is more than 50% of N.tInf.max, which -- in case this number is extrapolated -- might cause problems.\n")
  }
  #Create a vector representing the support of N.tInf
  N.tInf.support <- 0:control$N.tInf.max

  #======================================================================
  #======================================================================
  # Build reporting triangle and derived parameters for delay
  #======================================================================
  #======================================================================

  cat("Building reporting triangle...\n")

  #Time origin t_0
  t0 <- min(dateRange)
  #Sequence from time origin until now (per day??)
  #@hoehle
  t02s <- seq(t0,now,by=aggregate.by)
  #Maximum time index
  T <- length(t02s)-1

  #Check if the maximum delay is longer than the available time series
  if (D>T) {
    stop("D>T. Cannot estimate the long delays.")
  }
  #How many observations to take for estimating the delay distribution
  if (is.null(m)) {
    m <- T
  }
  if (m<1) { stop("Assertion m>=1 not fullfilled.") }

  #Define the observation triangle
  n <- matrix(NA,nrow=T+1,ncol=T+1,dimnames=list(as.character(t02s),NULL))

  #Loop over time points. (more efficient that delay and then t)
  for (t in 0:T) {
    #Extract all reports happening at time (index) t.
    #@hoehle: data.att <- data.sub[na2FALSE(data.sub[,dEventCol] == t02s[t+1]), ]
    data.att <- data.sub[na2FALSE(timeDelay(data.sub[,dEventCol], t02s[t+1])) == 0, ]
    #Loop over all delays
    for (x in 0:(T-t)) {
      #Count number with specific delay
      n[t+1,x+1] <- sum(data.att[,"delay"] == x)
    }
  }

  cat("No. cases: ",sum(n,na.rm=TRUE),"\n")

  #Handle delays longer than D
  #@hoehle: Not done!
  nLongDelay <- apply(n[,(D+1)+seq_len(T-D)],1,sum,na.rm=TRUE)
  if (any(nLongDelay>0)) {
    warning(paste(sum(nLongDelay)," cases with a delay longer than D=",D," days forced to have a delay of D days.\n",sep=""))
    n <- n[,1:(D+1)]
    n[,(D+1)] <- n[,(D+1)] + nLongDelay
  } else {
    #No problems. Just extract up to D+1
    n <- n[,1:(D+1)]
  }

  #Calculate n.x and N.x as in (2.7) and (2.8) and Fig.2 of Lawless (1994)
  #Note the different moving window definition as in the Lawless article.
  n.x <- rep(0,times=D+1)
  N.x <- rep(0,times=D+1)
  for (x in 0:D) {
    for (t in max(0,T-m):(T-x)) { #hoehle: Lawless definition is max(0,T-x-x)
      #cat("x=",x,"\tt=",t,":\n")
      n.x[x+1] <- n.x[x+1] + n[t+1,x+1]
      for (y in 0:x) {
        #cat("x=",x,"\tt=",t,"\ty=",y,":\n")
        N.x[x+1] <- N.x[x+1] + n[t+1,y+1]
      }
    }
  }

  cat("No. cases within moving window: ",sum(n.x,na.rm=TRUE),"\n")

  #Available observations at time T, definition of N(t;T) on p.17.
  N.tT <- sapply(0:T, function(t) sum(n[t+1, 0:min(D+1,(T-t)+1)]))

  #Truth - already in another object. Delete??
  N.tInf <- table( factor(as.character(data[,dEventCol]),levels=as.character(t02s)))

  #Store results of the reporting triangle in the control slot together with additional
  #attributes for fast access of, e.g., summaries or defining variables.
  reportingTriangle <- n
  attr(reportingTriangle, "n.x") <- n.x
  attr(reportingTriangle, "N.x") <- N.x
  attr(reportingTriangle, "N.tT") <- N.tT
  attr(reportingTriangle, "N.tInf") <- N.tInf
  attr(reportingTriangle, "T") <- T
  attr(reportingTriangle, "D") <- D
  attr(reportingTriangle, "t02s") <- t02s
  sts@reportingTriangle <- reportingTriangle

  #======================================================================
  # Calculations are jointly for all t values.
  #======================================================================

  #List of casts each containing a table 0..N.tInf.max with the PMF
  Ps <- list()

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #
  # Lawless (1994) method without adjustment for overdispersion
  #
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if ("lawless" %in% method) {
    #Hazard function estimates, i.e. g-function estimate as in (2.9)
    #of Lawless (1994). NAs are set to zero (consequences??)
    g.hat <- ifelse( !is.na(n.x/N.x), n.x/N.x, 0)
    #Force g.hat(0)=1 as stated just below (2.1)
    g.hat[1] <- 1

    #Check how the estimated CDF looks
    #F <- NULL ; for (d in 0:D) { i <- d+seq_len(D-d) ; F[d+1] <- prod(1-g.hat[i+1]) }
    #plot(0:D,F)

    #Compute weights Wt.hat as in eqn. (2.13). Use T1=Inf.
    #Note: Wt.hat estimates F_t(T-t).
    T1 <- Inf
    What.t <- sapply(0:T, function(t) {
      if (t<T-D) {
        1
      } else {
        x = T-t + seq_len(min(T1-t,D)-(T-t)) #(2.13) with modifications. knifflig
        prod(1-g.hat[x+1])
      }
    })

    #Store result of delay distribution estimation for support 0:T
    delayCDF[["lawless"]] <- rev(What.t)

    #Do the prediction as in (2.12)
    Nhat.tT1 <- N.tT / What.t

    #V.Wt as in (2.15)
    Vhat.Wt <- sapply(0:T, function(t) {
      if (t<T-D) {
        0
      } else {
        x = T-t + seq_len(min(T1-t,D)-(T-t))
        What.t[t+1]^2 * sum( g.hat[x+1]/(N.x[x+1]*(1-g.hat[x+1])),na.rm=TRUE)
      }
    })
    #(2.16)
    Vhat.Zt <- sapply(0:T, function(t) {
      if (t<T-D) {
        0
      } else {
        (N.tT[t+1]*(N.tT[t+1]+1))/What.t[t+1]^4*Vhat.Wt[t+1]
        + N.tT[t+1]*(1-What.t[t+1])/What.t[t+1]^2
      }
    })

    #Upper and lower 95% prediction limits based on the asymptotic normal
    #conf.level <- 0.95
    #U <- Nhat.tT1 + qnorm( (1+conf.level)/2)* sqrt(Vhat.Zt)
    #L <- pmax(N.tT, Nhat.tT1 - qnorm( (1+conf.level)/2)* sqrt(Vhat.Zt))

    #Discretize result: all mass below actual observation is put to observation
    PMFs <- matrix(NA, nrow=control$N.tInf.max+1,ncol=T+1)
    #CDF of a left truncated normal, truncated at "at"
    ltruncpnorm <- function(x, mean, sd, at) {
      ifelse( x < at, 0, pnorm(x, mean, sd) / (1-pnorm(at, mean,sd)))
    }
    #Deduce PMF for each time point in the past
    for (i in 1:length(Nhat.tT1)) {
      #Safeguard against NAs
      if (is.na(Vhat.Zt[i])) { warning("Vhat.Zt[i] is NA") ; Vhat.Zt[i] <- -1e99 }
      #If Vhat.Zt = 0 then Nhat.tT1 \equiv 0! Special care needed here
      if (Vhat.Zt[i] > 0) {
        CDF <- c(0,ltruncpnorm(N.tInf.support, mean=Nhat.tT1[i], sd=sqrt(Vhat.Zt[i]),at=N.tT[i]))
        PMFs[,i] <- diff(CDF)
      } else {
        #@hoehle: previous bug: c(1,rep(0,control$N.tInf.max)) ##all mass concentrated in zero, but it should be: Nhat.tT1
        PMFs[,i] <- (N.tInf.support == Nhat.tT1[i])*1
      }
    }

    Ps[["lawless"]] <- PMFs
  } #end lawless procedure


  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #
  # Bayesian method (simple model, clever sampling -> no MCMC)
  #
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  #Part jointly for both bayes and bayes.notrunc
  if (("bayes.trunc" %in% method) | ("bayes.notrunc" %in% method)) {
    cat("bayes prep...\n")

    ######################################################################
    # Prior of N(t,\infty)
    ######################################################################

    N.tInf.prior <- control$N.tInf.prior
    #Extract prior parameters from prior choice
    if (N.tInf.prior == "pois") {
      lambda <- attr(N.tInf.prior,"lambda",exact=TRUE)
    } else {
      if (N.tInf.prior == "poisgamma") {
        #Find size parameters such that mean variance is as target.
        var.prior <- function(size.prior) { mean.prior + mean.prior^2/size.prior }
        #If mean & var specified
        if (all(c("mean.lambda","var.lambda") %in% names(attributes(N.tInf.prior)))) {
          mean.prior <- attr(N.tInf.prior,"mean.lambda",exact=TRUE)
          var.prior.target <- attr(N.tInf.prior,"var.lambda",exact=TRUE)

          size.prior <- uniroot( function(size.prior) { var.prior(size.prior) - var.prior.target},interval=c(1e-12,50))$root
          #Check result
          cat("(E,V) of prior for lambda = (",paste(c(mean.prior,var.prior(size.prior)),collapse=","),")\n")
        } else {
          stop("mean.lambda and var.lambda not part of prior specification")
        }
      } else {
        if (N.tInf.prior == "unif") {
          N.tInf.prior.max <- attr(N.tInf.prior,"N.tInf.prior.max",exact=TRUE)
        } else {
          #No option applied
          stop("Not a valid prior!")
        }
      }
    }

    ######################################################################
    # Define function to generate PMF for max(0,T-D),..,T by sampling.
    #
    # Parameters:
    #  alpha.star, beta.star - vector containing the posterior GD params
    ######################################################################
    pmfBySampling <- function(alpha.star, beta.star) {
      #Sample from posterior distribution, i.e. sample from the reverse distribution
      #and reverse result
      p.sample <- rgd(control$nSamples,alpha.star,beta.star)[,(length(alpha.star)+1):1]

      #All the time points where extrapolation is to be done
      tSet <- max(0,(T-D)):T

      ######################################################################
      # Procedure to generate nowcasts of all time points up to T-D,...,T.
      # This is based on the posterior samples available in p.sample.
      # Current code adds up the PMF tables instead of a pure sample based
      # procedure and also prevents PMF=0 better than tabulating the samples.
      ######################################################################

      N.tT1.pred <- array(NA, dim=c(dim(p.sample)[1],control$N.tInf.max+1,dim(p.sample)[2]),dimnames=list(NULL,seq_len(control$N.tInf.max+1)-1L,tSet))
      for (j in 1:control$nSamples) {
        #Extract delay PMF from sample
        p <- p.sample[j,]
        #Proportion reported up to x, x=0,..,T
        F <- c(rep(1,T-D),rev(cumsum(p)))
        #Guard against numerical instability: ensure that not larger than 1.
        F <- ifelse(F>1,1,F)

        #Loop over all time points to nowcast
        for (i in 1:length(tSet)) {
          t <- tSet[i]
          N.tT1.pred[j,,i] <- switch(N.tInf.prior,
                              "poisgamma"=dpost.bnb(N.tT[t+1],sumpd=F[t+1],mu=mean.prior,size=size.prior,N.tInf.max=control$N.tInf.max))
        }
      }
      #Average the PMFs as in Step (2) of the algorithm
      PMF <- apply(N.tT1.pred,MARGIN=c(2,3),mean)

      #Add part, where no prediction needs to be done
      if (T-D>0) {
        #Empty PMFs
        determined <- matrix(0,nrow=control$N.tInf.max+1,ncol=T-D-1+1)
        #Add "1" entry at the observed
        for (t in 0:(T-D-1)) {
          determined[N.tT[t+1]+1,t+1] <- 1
        }
        PMF <- cbind(determined,PMF)
      }
      return(PMF)
    } #done definition of pmfBySampling
  }

  if ("bayes.trunc" %in% method) {
    cat("bayes.trunc...\n")

    ######################################################################
    #Prior of reporting delay as parameters of generalized Dirichlet prior
    ######################################################################

    #Define symmetric dirichlet as prior, just as in the other case
    alpha.prior <- rep(control$gd.prior.kappa, D)
    beta.prior <- rep(0,D)
    beta.prior[D] <- control$gd.prior.kappa
    for (i in (D-1):1) {
      beta.prior[i] <- alpha.prior[i+1] + beta.prior[i+1]
    }


    ######################################################################
    # Posterior section
    ######################################################################

    #Deduce posterior distribution of delay distribution, i.e. it is again
    #a generalized Dirichlet
    alpha <- beta <- rep(NA,D)
    for (d in 0:(D-1)) {
      alpha[d+1] <- n.x[D-d+1] ##Note: +1 coz index 1 is delay 0.
      beta[d+1] <-  N.x[D-d+1] - n.x[D-d+1]
    }

    #Check if there are any points without data and warn about it.
    if (any(alpha + beta == 0)) {
      warning("The delays ",paste(D-which(alpha+beta==0)-1,collapse=",")," have no observations. Results might be instable and depend all on prior.")
    }

    #Add up. Note: Delay zero (i.e. element D+1) is ignored as this is
    #not modelled explicitily by the GD distribution (sum to 1 constraints)
    alpha.star <- alpha.prior + alpha
    beta.star  <- beta.prior  + beta

    #Compute the expectation of the GD distribution and store this as the delay
    delayCDF[["bayes.trunc"]] <- cumsum(rev(Egd(alpha.star,beta.star)))

    #Save result
    Ps[["bayes.trunc"]] <- pmfBySampling(alpha.star, beta.star)

  } # end "bayes.trunc" %in% method

  #======================================================================
  # Bayesian version which ignores truncation
  #======================================================================

  if ("bayes.notrunc" %in% method) {
    cat("bayes.notrunc...\n")

    ######################################################################
    # Prior section
    ######################################################################
    alpha.prior <- rep(control$gd.prior.kappa, D) #symmetric dirichlet
    beta.prior <- rep(0,D)
    beta.prior[D] <- control$gd.prior.kappa
    for (i in (D-1):1) {
      beta.prior[i] <- alpha.prior[i+1] + beta.prior[i+1]
    }

    ######################################################################
    # Posterior section
    ######################################################################

    #Deduce posterior distribution of delay distribution, i.e. it is again
    #a generalized Dirichlet
    alpha <- beta <- rep(NA,D)
    for (d in 0:(D-1)) {
      alpha[d+1] <- n.x[D-d+1]
      beta[d+1] <-  sum(n.x[D - (d+1):D + 1])
    }

    #Check if there are any points without data and warn about it.
    if (any(alpha + beta == 0)) {
      warning("The delays ",paste(D-which(alpha+beta==0)-1,collapse=",")," have no observations. Results might be instable and depend all on prior.")
    }

    #Posterior parameters.
    alpha.star <- alpha.prior + alpha
    beta.star  <- beta.prior  + beta

    #Check that its a ordinary Dirichlet
    for (i in (D-1):1) {
      if (!all.equal(beta.star[i], alpha.star[i+1] + beta.star[i+1])) {
        warning("Posterior at i=",i," is not an ordinary Dirichlet as it's supposed to be.")
      }
    }

    #Save resulting delay distribution
    delayCDF[["bayes.notrunc"]] <- cumsum(rev(Egd(alpha.star,beta.star)))

    Ps[["bayes.notrunc"]] <- pmfBySampling(alpha.star,beta.star)

  } # end bayes.notrunc

  ######################################################################
  # Unadjusted procedure using beta negative binomial. ToDo:
  # integrate code with other Bayesian procedures
  ######################################################################
  if ("bayes.notrunc.bnb" %in% method) {
    cat("bayes.notrunc.bnb...\n")

    ######################################################################
    # Prior section (same as for all methods)
    ######################################################################
    alpha.prior <- rep(control$gd.prior.kappa, D) #symmetric dirichlet
    beta.prior <- rep(0,D)
    beta.prior[D] <- control$gd.prior.kappa
    for (i in (D-1):1) {
      beta.prior[i] <- alpha.prior[i+1] + beta.prior[i+1]
    }

    ######################################################################
    # Posterior section
    ######################################################################

    #Deduce posterior distribution of delay distribution, i.e. it is again
    #an ordinary Dirichlet
    alpha <- beta <- rep(NA,D)
    for (d in 0:(D-1)) {
      alpha[d+1] <- n.x[D-d+1]
      beta[d+1] <-  sum(n.x[D - (d+1):D + 1])
    }

    #Check if there are any points without data and warn about it.
    if (any(alpha + beta == 0)) {
      warning("The delays ",paste(D-which(alpha+beta==0)-1,collapse=",")," have no observations. Results might be instable and depend all on prior.")
    }

    #Posterior parameters.
    alpha.star <- alpha.prior + alpha
    beta.star  <- beta.prior  + beta

    #Check that its a ordinary Dirichlet
    for (i in (D-1):1) {
      if (!all.equal(beta.star[i], alpha.star[i+1] + beta.star[i+1])) {
        warning("Posterior at i=",i," is not an ordinary Dirichlet as it's supposed to be.")
      }
    }

    #Save resulting delay distribution (i.e. no truncation adjustment)
    delayCDF[["bayes.notrunc"]] <- cumsum(rev(Egd(alpha.star,beta.star)))

    #Allocate PMF to return
    PMF <- matrix(0,nrow=control$N.tInf.max+1,ncol=length(max(0,(T-D)):T))

    #Concentration parameter vector of the ordinary Dirichlet distribution
    #Note. alpha.star vector is reversed (shortest delay last).
    alpha <- rev(c(alpha.star,beta.star[length(beta.star)]))

    #consistency check
    if (!all.equal(rev(Egd(alpha.star,beta.star)),alpha/sum(alpha))) {
      stop("Problem. GD and Dirichlet do not correspond...")
    }

    tSet <- max(0,(T-D)):T
    for (i in 1:length(tSet)) {
      t <- tSet[i]
      alpha.i <- cumsum(alpha)[T-t+1]
      beta.i <-  sum(alpha) - alpha.i
      if (T-t==D) {
        PMF[,i] <- ifelse( N.tInf.support  == N.tT[t+1], 1, 0)
      } else {
        #Calculate PMF knowing the q ~ Beta( , ) by the aggregation
        #property.
        #Note: Vector N.tT starts at time zero, i.e. time T corresponds to T+1
        PMF[,i] <- dbnb( N.tInf.support - N.tT[t+1],n=N.tT[t+1]+1, alpha=alpha.i, beta=beta.i)
      }
    } #done looping over all time points

    #Add part, where no prediction needs to be done
    if (T-D>0) {
      #Empty PMFs
      determined <- matrix(0,nrow=control$N.tInf.max+1,ncol=T-D-1+1)
      #Add "1" entry at the observed
      for (t in 0:(T-D-1)) {
        determined[N.tT[t+1]+1,t+1] <- 1
      }
      PMF <- cbind(determined,PMF)
    }

    Ps[["bayes.notrunc.bnb"]] <- PMF
  } # end bayes.notrunc.bnb


  ######################################################################
  # Fully Bayes version with MCMC
  ######################################################################

  if ("bayes.trunc.ddcp" %in% method) {
      #Allocate result
      PMF <- matrix( 0,ncol=(T+1),nrow=control$N.tInf.max+1)

      #Fix seed value of JAGS RNG for each chain
      n.chains <- 3
      init <- lapply(1:n.chains,function(i) {
          list(.RNG.name="base::Mersenne-Twister",.RNG.seed=i*10)
      })

      #Make design matrix for a quadratic TPS spline in time
      makeTPSDesign <- function(T,degree=2) {
          nbeta=degree + 1
          X <- matrix(NA,ncol=nbeta, nrow=T+1)
          for (t in 0:T) {
          #Form a centered time covariate
              t.centered <- t - T/2
              for(pow in 0:degree) {
                  X[t+1,pow+1]<- t.centered^(pow)
              }
          }

          #Make the knot points evenly spaced between 0,T not including these points
          knots <- seq(0,T,length=min(round(T/6)+2,22))
          knots <- knots[-c(1,length(knots))]
          #Remove knots which are beyond T-maxDelay/2
          knots <- knots[knots <= T-D/2]
          knots <- knots - T/2
          nknots <- length(knots)

          #Penalty as REs - setup design matrix
          Z <- matrix(NA,nrow=T+1,ncol=length(knots))
          for (t in 0:T){
              t.center <- t - T/2
              for(k in 1:nknots){
                  Z[t+1,k]<- pmax((t.center-knots[k]),0)^degree
              }
          }
          return(list(X=X,Z=Z,knots=knots,nknots=nknots,nbeta=nbeta))
      }
      tps <- makeTPSDesign(T=T,degree=2)

      #Design matrix for logistic discrete time hazard model containing
      #changepoints. Could be extended s.t. the user provides W.
      W <- array(NA,dim=c(T+1,length(ddChangepoint),D+1),dimnames=list(as.character(t02s),as.character(ddChangepoint),paste("delay",0:D,sep="")))
      for (t in 0:T){
          for (i in 1:length(ddChangepoint)) {
              W[t+1,i,] <- as.numeric( (t02s[t+1] + 0:D) >= ddChangepoint[i])
          }
      }

      #Priors. Uniform on the delays
      D.prime <- round( D/2-0.4)+1
      p.prior <- rep(1/(D.prime+1), D.prime+1)
      mu.gamma <- qlogis( p.prior[1])
      for (d in 1:(D.prime-1)) {
          mu.gamma <- c(mu.gamma, qlogis( p.prior[d+1] / (1-sum(p.prior[1:d]))))
      }
      tau.gamma <- rep(control$ddcp$tau.gamma,times=D.prime)

      #Prepare data for JAGS
      jagsData <- list(#Data
                       rT=n,T=T+1,m=m+1,maxDelay=D,
                       #Time dependent logistic discrete hazard model
                       W=W, eta.mu=control$ddcp$eta.mu, eta.prec=control$ddcp$eta.prec,
                       mu.gamma=mu.gamma, tau.gamma=tau.gamma,
                       #Epidemic curve
                       alpha.lambda=2500/3000,beta.lambda=50/3000,
                       #Spline related stuff
                       X=tps$X,Z=tps$Z,nknots=tps$nknots,beta.mu=rep(0,tps$nbeta),beta.prec=1e-6*diag(tps$nbeta)
                       )


      #Select appropriate model (change this to be part of the options!!)
      logLambda.method <- control$ddcp$logLambda   #"tps" #"rw2" #"iid" #"rw2" #"rw2" #"iid" #"rw" #"tps"

###      browser()

      #Load the BUGS specification of the Bayesian hierarchical Poisson model
      bugsModel <- readLines(file.path(path.package('surveillance'),'jags',"bhpm.bugs"))

      bugsModel <- gsub(paste("#<",logLambda.method,">",sep=""),"",bugsModel)
      #Problem when eta is scalar (TODO: Improve!!)
      if (length(ddChangepoint) == 1) {
          #Make eta ~ dnorm( , ) instead of eta ~ dmnorm
          bugsModel <- gsub("(^[ ]*eta ~ )(dmnorm)","\\1dnorm",bugsModel)
          #Use eta[1] instead of eta for matrix multiplication
          bugsModel <- gsub("(eta)(.*%\\*%)","eta\\[1\\]\\2",bugsModel)
      }
      #cat(paste(bugsModel,collapse="\n"))
      bugsFile <- tempfile(pattern = "nowcast-")
      writeLines(bugsModel,bugsFile)

      ##browser()
      ## if (FALSE) {
      ##     #Try to compile the model with ordinary rjags to see if there are any problems
      ##     #before doing 3 chains parallelized using runjags.
      ##     model <- jags.model(bugsFile,
      ##                         data = jagsData,
      ##                         init=init, #Fix seed value of JAGS as well
      ##                         n.chains = n.chains, n.adapt = 100)
      ##     list.samplers(model)
      ##     coda.samples(model,variable.names='logLambda',n.iter=100)
      ## }


      ######################################################################
      # runjags way -- ToDo: parametrize using control options!
      ######################################################################

      runjagsMethod <-  'rjparallel' #'rjags'
      monitor <- c('gamma','eta','logLambda','NtInf')
      samples.rj <- runjags::run.jags(bugsFile,#bugsModel,
                             monitor = monitor, data=jagsData, n.chains=3,
                             inits = init,
                             burnin =  control$ddcp$mcmc["burnin"],
                             sample =  control$ddcp$mcmc["sample"],
                             thin =  control$ddcp$mcmc["thin"],
                             adapt=1000,
                           summarise=FALSE,method=runjagsMethod)

      #Extract posterior median of discrete survival time delay distribution model parameters
      dt.surv.samples <- coda::as.mcmc.list(samples.rj, vars = c('gamma','^eta'))
      post.median <- dt.surv.pm <- apply( as.matrix(dt.surv.samples), 2, median)

      #Posterior median of the lambda's
      lambda.post <- exp(apply( as.matrix(coda::as.mcmc.list(samples.rj, vars = c('^logLambda'))), 2,
                               quantile, prob=c(0.025,0.5,0.975)))

      #Extract posterior median of model parameters
      gamma.red <- post.median[grep("gamma",names(post.median))]
      eta <- matrix(post.median[grep("^eta",names(post.median))])
      #Map from reduced set to full set
      gamma <- gamma.red[round( (0:(D-1))/2 - 0.4) + 1]

      #Compute the resulting PMF from the model. Possibly put this in separate function.
      pmf <- matrix(NA, nrow=nrow(W),ncol=D+1,dimnames=list(as.character(t02s),paste("delay",0:D,sep="")))
      #Determine PMF
      for (t in 1:length(t02s)) {
        if (as.character(t02s[t]) %in% rownames(W)) {
          lin.pred <- ( gamma + t(eta) %*% W[t,,0:D])
          pmf[t,] <- haz2pmf(c(plogis(lin.pred),1))
        }
      }

      #Store result as CDF
      delayCDF[["bayes.trunc.ddcp"]] <- t(apply(pmf, 1, cumsum))
      #Store model as attribute
      if(control$ddcp$logLambda != "tps") tps <- NULL
      attr(delayCDF[["bayes.trunc.ddcp"]],"model") <- list(post.median=dt.surv.pm,W=W,lambda.post=lambda.post,tps=tps)

      #Convert to coda compatible output.
      samples <- coda::as.mcmc.list(samples.rj)


      #Extract PMFs
      for (t in 0:T) {
          #Extract samples related to this time point
          vals <- as.matrix(samples[,paste("NtInf[",t+1,"]",sep="")])
          #PMF
          PMF[,t+1] <- prop.table(table(factor(vals,levels=0:control$N.tInf.max)))
      }
      Ps[["bayes.trunc.ddcp"]] <- PMF
  }

  #======================================================================
  #A really bad forecast -- the uniform
  #======================================================================

  if ("unif" %in% method) {
    #Allocate result
    PMF <- matrix( 0,ncol=(T+1),nrow=control$N.tInf.max+1)
    #Loop over all time points to nowcast and put U(N.tT[t],Nmax)
    for (t in 0:T) {
      #How many values are there in N.tT .. Nmax
      noVals <- max(0,control$N.tInf.max - N.tT[t+1]) + 1
      #PMF at t is 0,...0 (N.tT-1 times), 1/noVals,...,1/noVals
      PMF[,t+1] <- c(rep(0,N.tT[t+1]),rep(1/noVals,times=noVals))
    }
    Ps[["unif"]] <- PMF
  }

  ######################################################################
  #Loop over all time points in the vector "when". Only these are
  #returned.
  ######################################################################
  idxt <- which(dateRange %in% when)
  for (i in idxt) {
    #Save PMFs if thats requested
    if (control$predPMF) {
      res <- list()
      for (j in 1:length(method)) {
        res[[method[j]]] <- Ps[[method[j]]][,i]
      }
      sts@predPMF[[as.character(dateRange[i])]] <- res
    }


    #Evaluate scoring rules, if requested
    if (control$score) {
      #Infer the true value
      ytinf <- observed(sts.truth)[i,]
      #Evaluate all scores for all predictive distributions
      for (i.P in 1:length(method)) {
        for (i.score in 1:length(scores)) {
            #cat("i=",i," i.P=",i.P," (",method[i.P],") i.score=",i.score,"\n")
            SR[i,i.P,i.score] <- do.call(scores[i.score],args=list(P=Ps[[method[i.P]]][,i],y=ytinf,alpha=control$alpha))
        }
      }
    } #end if control$score

    #Add first nowcast & ci to stsNC slots
    sts@upperbound[i,] <- median(N.tInf.support[which.max( cumsum(Ps[[method[1]]][,i])>0.5)])
    sts@pi[i,,] <- N.tInf.support[c(which.max(cumsum(Ps[[method[1]]][,i]) > control$alpha/2),which.max(cumsum(Ps[[method[1]]][,i]) > 1-control$alpha/2))]
    dimnames(sts@pi) <- list(as.character(dateRange),NULL,paste( c(control$alpha/2*100,(1-control$alpha/2)*100),"%",sep=""))
  } #end of loop over time points


  #Add scoring rule to output
  if (control$score) {
    dimnames(SR)    <- list(as.character(dateRange),method,scores)
    sts@SR <- SR
  }

  ######################################################################
  #Other arguments to save in control object
  ######################################################################
  sts@control$N.tInf.support <- N.tInf.support
  sts@control$method <- sts@control$name <- method
  #Store variables relevant for the nowcast
  sts@control$D <- D
  sts@control$m <- m
  sts@control$now <- now
  sts@control$when <- when
  sts@control$timeDelay <- timeDelay

  #Store delayCDF object
  sts@delayCDF <- delayCDF

  #For backwards compatibility -- change this in the future TODODODODODO!
  sts@control$yt.support <-  sts@control$N.tInf.support
  sts@control$y.prior.max <- sts@control$N.tInf.max

  #Done
  return(sts)
}



######################################################################
# Helper functions
######################################################################

#Helper function
na2FALSE <- function(x) {x[is.na(x)] <- FALSE ; return(x) }

######################################################################
# Logarithmic score
#
# Parameters:
#  P - predictive distribution, given as a vector containing the PMF
#      with support 0,...,N.prior.max
#  y - the actual observation. Can be a vector.
#
# Returns:
#  -log P(y). If y outside 0,..,N.prior.max then -Inf.
######################################################################

logS <- function(P, y, ...) {
  return(ifelse( y>=0 & y<=length(P)-1, -log(P[y+1]), -Inf))
}

######################################################################
# Ranked probability score
#
# Parameters:
#  P - predictive distribution, given as a vector containing the PMF
#      with support 0,...,N.prior.max
#  y - the actual observation. Can be a vector.
#
# Returns:
#  -log P(y). If y outside 0,..,N.prior.max then -Inf.
######################################################################

RPS <- function(P,y, ...) {
  N.support <- 0:(length(P)-1)
  sum( (cumsum(P) -  (y <= N.support))^2)
}

#Some other scoring rules which are not proper.
dist.median <- function(P,y, ...) {
  point.estimate <- which.max(cumsum(P)>=0.5) - 1
  return(abs(point.estimate - y))
}

#0/1 indicator of observed value outside equal tailed (1-alpha/2) CI
outside.ci <- function(P,y,alpha) {
  N.support <- 0:(length(P)-1)
  ci <- N.support[c(which.max(cumsum(P) > alpha/2),which.max(cumsum(P) >
1-alpha/2))]
  ifelse( y>=ci[1] & y<=ci[2], 0, 1)
}



######################################################################
# Helper functions for sampling the predictive distribution
######################################################################

#Unnormalized in Binomial-Negative-Binomial Hierarchy. Should work for vectors of N.tInf!
#Only kernel parts for N.tInf need to be taken into account
dpost.bnb.unorm <- function(N.tInf, N.tT, sumpd, mu, size) {
  dbinom(N.tT, size=N.tInf, prob=sumpd)*dnbinom(N.tInf, mu=mu,size=size)
  #Direct implementation - appears to be less stable...
  #ifelse(N.tInf >= N.tT,
  #       exp(lgamma(N.tInf+size)-lgamma(N.tInf-N.tT+1) + N.tInf*log( (1-sumpd)*(mu/(mu+size)))),0)
  #Compare the 2
  ## foo.a <-  dbinom(N.tT, size=N.tInf, prob=sumpd)*dnbinom(N.tInf, mu=mu,size=size)
  ## foo.b <- ifelse(N.tInf >= N.tT, #& N.tInf <= size,
  ##                 exp(lgamma(N.tInf+size)-lgamma(N.tInf-N.tT+1) + N.tInf*log( (1-sumpd)*(mu/(mu+size)))),0)
  ## plot(foo.a/sum(foo.a))
  ## points(foo.b/sum(foo.b),col="red")
}

#Sample in binomial-negative-binomial hierarchy
rpost.bnb <- function(n=1, N.tT, sumpd, mu,size, N.tInf.max=1e4) {
  p <- dpost.bnb.unorm(0:N.tInf.max,N.tT=N.tT,sumpd=sumpd, mu=mu,size=size)
  #Set NA values to zero (why would they be NA?)
  #if (is.na(sum(p))) { warning("rpost.bnb: sum is NA") ; browser(p)}
  #Normalize the distribution - safe this for time reasons
  #p <- p/sum(p)
  #Sample
  sample(0:N.tInf.max, size=n, replace=TRUE, prob=p)
}

#PMF for the predictive distribution in binomial-negative-binomial hierarchy.
#Returns entire vector for 0:N.tInf.max
dpost.bnb <- function(N.tT, sumpd, mu,size, N.tInf.max=1e4) {
  p <- dpost.bnb.unorm(0:N.tInf.max,N.tT=N.tT,sumpd=sumpd, mu=mu,size=size)
  #Set NA values to zero (why would they be NA?)
  #if (is.na(sum(p))) { warning("rpost.bnb: sum is NA") ; browser(p)}
  #Normalize the distribution - safe this for time reasons
  return(p/sum(p))
}


######################################################################
# PMF of the beta-negatative binomial distribution
# See Teerapabolarn (2008)
#
# Parameters:
#   k - where to evaluate. can be a vector.
#
# Returns:
# PMF.
######################################################################

dbnb <- function(k,n,alpha,beta) {
  #Check if k's outside the support are requested.
  neg <- k<0
  k[neg] <- 0
  #Calculate the density of the beta-negbin. See Teerapabolarn (2008)
  num <- lgamma(n+alpha)+lgamma(k+beta)+lgamma(n+k)+lgamma(alpha+beta)
  den <- lgamma(n+k+alpha+beta)+lgamma(n)+lgamma(k+1)+lgamma(alpha)+lgamma(beta)
  res <- exp(num-den)
  res[neg] <- 0
  return( res)
}


######################################################################
# Convert discrete time hazard function on 0,...,Dmax to a probability
# mass function.
#
# Parameters:
#  haz - vector with entries for (0,...,Dmax)
# Returns:
#  vector with PMF on 0,...,Dmax.
######################################################################

haz2pmf <- function(haz) {
    PMF <- 0*haz
    for (i in 0:(length(haz)-1)) {
        PMF[i+1] <- haz[i+1] * (1-sum(PMF[seq(i)]))
    }
    return(PMF)
}
