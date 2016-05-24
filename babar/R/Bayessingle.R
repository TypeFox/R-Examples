########################################################################
# Bayessingle.R - code to perform Bayesian analysis for fitting a single
# bacterial growth curve
# (c) Lydia Rickett
########################################################################

.logllSingle <- function(params,par.no,modelfunc,model,dataset,inc.nd,threshold,t.nd,inf.sigma,sigma,transformParams,mumax.prior) {
  # To compute the log likelihood
  
  t <- dataset$t; y <- dataset$y
  ymax <- max(y); ymin <- min(y); tmax <- max(t); ydiff <- (ymax-ymin)
  tParams <- transformParams(params)
  
  if (inf.sigma) {
    sigma <- rep(tParams[par.no+1],length=length(y))
  }
  
  # If using a Cauchy prior for mumax, penalise values outside of a sensible range in the
  # likelihood function (to prevent guesses that are negative or too large)
  if ( ((mumax.prior == "Cauchy") && (model == "linear" || model == "Bar3par") && 
        (tParams[2] < 0 || tParams[2] > 10 * ydiff/tmax)) ||
    ((mumax.prior == "Cauchy") && (model == "logistic" || model == "Bar4par" || model == "Bar6par") 
     && (tParams[3] < 0 || tParams[3] > 10 * ydiff/tmax)) ){
    llhood <- -10^4
  } else {
    fun <- modelfunc(t,tParams[1:par.no])
    res <- sum(((y-fun)/sigma)^2)
    llhood <- log(1/prod(sigma*(2*pi)^(1/2)))-res/2
    
    if (inc.nd && length(t.nd > 0)) {
      fun.nd <- modelfunc(t.nd,tParams[1:par.no])
      # Approximate undetected values using a uniform distribution
      P.nd <- rep(1/threshold,length=length(fun.nd))
      P.d <- numeric(length=length(fun.nd))
      for (k in 1:length(t.nd)){
        if ((fun.nd[k] >= 0) && (fun.nd[k] <= threshold)) {
          P.d[k] <- 0 # In undetected region, uniform distribution
        } else {         
          # In uncensored region, add Gaussian tail with very small variance to avoid negative infinity in log likelihood
          P.d[k] <- ((fun.nd[k]-threshold)/10^(-10))^2 
        }
      }
      llhood <- llhood + log(prod(P.nd)) - sum(P.d)
    }
  }
  
  return (llhood)
}

.generateTransformSingle <- function(dataset,par.no,model,inc.nd,inf.sigma,mumax.prior,mu.mean,mu.sd) {
  # A wrap around to include other input
  
  t <- dataset$t; y <- dataset$y
  ymax <- max(y); ymin <- min(y); ydiff <- (ymax-ymin)
  tmax <- max(t)
  
  transformParams <- function(uParams) {
    # Transform from the unit hypercube to uniform/Jeffreys' prior
    
    tParams = numeric(length = length(uParams))
    
    # When t(0)!=0 or including undetected values, use lower bound y(0)=0 instead, or lower bound might not be low enough
    if (t[1] != 0 | inc.nd){ 
      fact <- 0
    } else {
      fact <- 1
    }
    
    if (par.no == 4) {
      tParams[1] = UniformPrior(uParams[1], fact*(ymin - ydiff/2), ymin + ydiff/2)
      tParams[2] = UniformPrior(uParams[2], ymax - ydiff/2, ymax + ydiff/2)
      if (mumax.prior == "Gaussian") {
        tParams[3] = GaussianPrior(uParams[3], mu.mean, mu.sd)
      } else if (mumax.prior == "Cauchy") {
        tParams[3] = CauchyPrior(uParams[3], mu.mean, mu.sd)
      } else { # mumax.prior = Uniform
        tParams[3] = UniformPrior(uParams[3], 0, 10 * ydiff/tmax)
      }
      tParams[4] = UniformPrior(uParams[4], 0, 9 * ydiff)
    } else if (par.no == 3) {
      if (model == "logistic"){
        tParams[1] = UniformPrior(uParams[1], fact*(ymin - ydiff/2), ymin + ydiff/2)
        tParams[2] = UniformPrior(uParams[2], ymax - ydiff/2, ymax + ydiff/2)
        if (mumax.prior == "Gaussian") {
          tParams[3] = GaussianPrior(uParams[3], mu.mean, mu.sd)
        } else if (mumax.prior == "Cauchy") {
          tParams[3] = CauchyPrior(uParams[3], mu.mean, mu.sd)
        } else { # mumax.prior = Uniform
          tParams[3] = UniformPrior(uParams[3], 0, 10 * ydiff/tmax)
        }
      } else {
        tParams[1] = UniformPrior(uParams[1], fact*(ymin - ydiff/2), ymin + ydiff/2)
        if (mumax.prior == "Gaussian") {
          tParams[2] = GaussianPrior(uParams[2], mu.mean, mu.sd)
        } else if (mumax.prior == "Cauchy") {
          tParams[2] = CauchyPrior(uParams[2], mu.mean, mu.sd)
        } else { # mumax.prior = Uniform
          tParams[2] = UniformPrior(uParams[2], 0, 10 * ydiff/tmax)
        }
        tParams[3] = UniformPrior(uParams[3], 0, 9 * ydiff)
      }
    } else if (par.no == 2){
      tParams[1] = UniformPrior(uParams[1], fact*(ymin - ydiff/2), ymin + ydiff/2)
      if (mumax.prior == "Gaussian") {
        tParams[2] = GaussianPrior(uParams[2], mu.mean, mu.sd)
      } else if (mumax.prior == "Cauchy") {
        tParams[2] = CauchyPrior(uParams[2], mu.mean, mu.sd)
      } else { # mumax.prior = Uniform
        tParams[2] = UniformPrior(uParams[2], 0, 10 * ydiff/tmax)
      }
    } else { # No. parameters = 6
      tParams[1] = UniformPrior(uParams[1], fact*(ymin - ydiff/2), ymin + ydiff/2)
      tParams[2] = UniformPrior(uParams[2], ymax - ydiff/2, ymax + ydiff/2)
      if (mumax.prior == "Gaussian") {
        tParams[3] = GaussianPrior(uParams[3], mu.mean, mu.sd)
      } else if (mumax.prior == "Cauchy") {
        tParams[3] = CauchyPrior(uParams[3], mu.mean, mu.sd)
      } else { # mumax.prior = Uniform
        tParams[3] = UniformPrior(uParams[3], 0, 10 * ydiff/tmax)
      }
      tParams[4] = UniformPrior(uParams[4], 0, 0.9*tmax)
      tParams[5] = UniformPrior(uParams[5], 0, 10/tmax)
      tParams[6] = UniformPrior(uParams[6], 0, 10/ymax)
    }
    
    if (inf.sigma){ # Extra parameter when inferring sigma - bounds between 0.01 and 10
      tParams[length(uParams)] = JeffreysPrior(uParams[length(uParams)], -2, 1)
    }

    return(tParams)
  }
 
  return(transformParams)
}

###
# Functions for growth models:
###

.linear <- function(t,params) {
  # Compute y(t) = ln(x(t)) using the linear model (2-parameters)
  #
  # Arguments: vector t of times and parameters par1 = y0 and par3 = mumax
  
  par1 <- params[1]; par3 <- params[2];
  
  y <- par1 + par3*t
  
  return (y=y)
}

.logistic <- function(t,params) {
  # Compute y(t) = ln(x(t)) using the logistic model (3-parameters, no lag phase)
  #
  # Arguments: vector t of times and parameters par1 = y0, par2 = ymax and par3 = mumax

  par1 <- params[1]; par2 <- params[2]; par3 <- params[3]
  
  Y2 <- rep(1,length(t))
  for (i in 1:length(t)) { 
    if (par3*t[i] > 700) { # Preventing overflow
      Y2[i] <- 0
    }
  }
  y <- par1 + Y2*(par3*t - log(1+(exp(par3*t)-1)/exp(par2-par1)))
  
  return (y=y)
}

.Bar3par <- function(t,params) {
  # Compute y(t) = ln(x(t)) using the 3-parameter Baranyi model (no stationary phase)
  #
  # Arguments: vector t of times and parameters par1 = y0, par3 = mumax and par4 = h0
  
  par1 <- params[1]; par3 <- params[2]; par4 <- params[3]
  
  a1 <- rep(1,length(t)); a2 <- rep(1,length(t)); a3 <- rep(1,length(t))
  for (i in 1:length(t)) {
    if (abs(par3*t[i]-par4) > 700) { # Preventing over/underflow
      if (-(par3*t[i]-par4) > 0) {
        a1[i] <- 0
      } else {
        a3[i] <- 0
      }
    }
    if (-par3*t[i] < -700) {
      a2[i] <- 0
    }
  }
  A <- a1*(t-par4/par3 + log(1-a2*exp(-par3*t)+a3*exp(-(par3*t-par4)))/par3)
  
  Y1 <- rep(1,length(t))
  for (i in 1:length(t)) { 
    if (par3*A[i] < 0) { # Preventing round-off error; mumax*A(t) must be positive
      Y1[i] <- 0
    }
  }
  muA <- Y1*par3*A
  y <- par1 + muA
  
  return (y=y)
}

.Bar4par <- function(t,params) {
  # Compute y(t) = ln(x(t)) using the 4-parameter Baranyi model
  #
  # Arguments: vector t of times and parameters par1 = y0, par2 = ymax, par3 = mumax and par4 = h0
  
  par1 <- params[1]; par2 <- params[2]; par3 <- params[3]; par4 <- params[4]
  
  a1 <- rep(1,length(t)); a2 <- rep(1,length(t)); a3 <- rep(1,length(t))
  for (i in 1:length(t)) {
    if (abs(par3*t[i]-par4) > 700) { # Preventing over/underflow
      if (-(par3*t[i]-par4) > 0) {
        a1[i] <- 0
      } else {
        a3[i] <- 0
      }
    }
    if (-par3*t[i] < -700) {
      a2[i] <- 0
    }
  }
  A <- a1*(t-par4/par3 + log(1-a2*exp(-par3*t)+a3*exp(-(par3*t-par4)))/par3)
  
  Y1 <- rep(1,length(t)); Y2 <- rep(1,length(t))
  for (i in 1:length(t)) { 
    if (par3*A[i] < 0) { # Preventing round-off error; mumax*A(t) must be positive
      Y1[i] <- 0
    }
    if (par3*A[i] > 700) { # Preventing overflow
      Y2[i] <- 0
    }
  }
  muA <- Y1*par3*A
  y <- par1 + Y2*(muA - log(1+(exp(muA)-1)/exp(par2-par1)))
  
  return (y=y)
}

.Bar6par <- function(t,params) {
  # Compute y(t) = ln(x(t)) using the 6-parameter Baranyi model
  #
  # Arguments: vector t of times and parameters par1 = y0, par2 = ymax, par3 = mumax, par4 = lambda, par5 = nu and par6 = m  

  par1 <- params[1]; par2 <- params[2]; par3 <- params[3]; par4 <- params[4]; par5 <- params[5]; par6 <- params[6]
  
  a1 <- rep(1,length(t)); a2 <- rep(1,length(t)); a3 <- rep(1,length(t))
  for (i in 1:length(t)) {
    if (abs(par5*(t[i]-par4)) > 700) { # Preventing over/underflow
      if (-par5*(t[i]-par4) > 0) {
        a1[i] <- 0
      } else {
        a3[i] <- 0
      }
    }
    if (-par5*t[i] < -700) {
      a2[i] <- 0
    }
  }
  A <- a1*(t-par4 + log(1-a2*exp(-par5*t)+a3*exp(-par5*(t-par4)))/par5)
  
  Y1 <- rep(1,length(t)); Y2 <- rep(1,length(t))
  for (i in 1:length(t)) { 
    if (par3*A[i] < 0) { # Preventing round-off error; mumax*A(t) must be positive
      Y1[i] <- 0
    }
    if (par6*par3*A[i] > 700) { # Preventing overflow
      Y2[i] <- 0
    }
  }
  muA <- Y1*par3*A
  y <- par1 + Y2*(muA - log(1+(exp(par6*muA)-1)/exp(par6*(par2-par1)))/par6)
  
  return (y=y)
}

.sortdata <- function(data,inc.nd) {
  # Sort out data for use in nested sampling
  
  t <- numeric(0)
  y <- numeric(0)
  t.nd <- numeric(0) # Times of undetected values
  
  for (i in 1:nrow(data)) {
    if(is.na(data[i,2]) == FALSE){ # Takes care of "NA"/"ND" values
      t <- c(t,data[i,1]); y <- c(y,data[i,2])} else {
        if (inc.nd) {
          t.nd <- c(t.nd,data[i,1])
        } else {
          warning("there are undetected points which are not being used as part of the analysis")
        }
      }
  }
  
  ty <- data.frame(t,y)
  
  return(list(ty=ty,t.nd=t.nd))
  ### Returns: the dataset minus undetected values, (t,y), and list of times of undetected values, t.nd
}

Bayesfit <- structure(function
  ### Perform Bayesian analysis for fitting a single bacterial growth curve using the Baranyi model.
  (data,
  ### A datafile of the curve to be fitted. This should consist of two columns, the first for time and second for logc. 
  ### The bacterial concentration should be given in log_10 cfu and there should be at least 2 data points (the first of which 
  ### may be undetected). Undetected y values should be represented by "NA".
  model,
  ### The growth model to be used. This should be one of "linear", "logistic", "Bar3par", "Bar4par" and "Bar6par".
  inf.sigma=TRUE,
  ### (TRUE/FALSE) Choose whether or not to infer the noise level, sigma, as part of the analysis. If FALSE, sigma should be
  ### specified (or the default value of sigma, 0.3, will be used).
  inc.nd=FALSE,
  ### Choose whether or not to include undetected points as part of the analysis. If TRUE, threshold should be specified.
  sigma=0.3,
  ### The choice of noise level, sigma, in log_10 cfu if it is not inferred as part of the analysis. Default is 0.3.
  threshold=NULL,
  ### Threshold in log_10 cfu below which values are considered as undetected.
  mumax.prior="Uniform",
  ### The type of prior to use for mu_max. This should be one of "Uniform", "Gaussian" or "Cauchy" (or the default "Uniform"
  ### will be used). If "Gaussian" or "Cauchy" are specified, mu.mean and mu.sd should be given. 
  mu.mean=NULL,
  ### The mean to be used when using a Gaussian or Cauchy prior.
  mu.sd=NULL,
  ### The standard deviation to be used when using a Gaussian or Cauchy prior.
  tol=0.1,
  ### The termination tolerance for the nested sampling
  prior.size=100
  ### The number of prior samples to use for nested sampling
) {

  # Setup (converting to natural log scale for use in the model)
  
  if(model!="linear" && model!="logistic" && model!="Bar3par" && model!="Bar4par" && model!="Bar6par") {
    stop("'model' must be one of 'linear', 'logistic', 'Bar3par', 'Bar4par' or 'Bar6par'")
  }
  
  if(mumax.prior!="Uniform" && mumax.prior!="Gaussian" && mumax.prior!="Cauchy") {
    stop("'mumax.prior' must be one of 'Uniform', 'Gaussian' or 'Cauchy'")
  }
  if((mumax.prior=="Gaussian" || mumax.prior=="Cauchy") && (is.null(mu.mean) || is.null(mu.sd))) {
    stop("both 'mu.mean' and 'mu.sd' must be given when using a Gaussian or Cauchy prior for mu_max")
  }
  mu.mean <- mu.mean*log(10); mu.sd <- mu.sd*log(10)
  
  if (model == "linear") {
    par.no <- 2
    modelfunc <- .linear
  } else if (model == "logistic") {
    par.no <- 3
    modelfunc <- .logistic
  } else if (model == "Bar3par") {
    par.no <- 3
    modelfunc <- .Bar3par
  } else if (model == "Bar4par") {
    par.no <- 4
    modelfunc <- .Bar4par
  } else if (model == "Bar6par") {
    par.no <- 6
    modelfunc <- .Bar6par
  }
  
  if (nrow(data) == 1) {
    stop("'data' must include at least 2 data points")
  }
  
  if ((nrow(data) == 2) && (inf.sigma == TRUE)) {
    warning("inferring noise level sigma is not recommended with only 2 data points, sigma has been prescribed as 0.3")
    inf.sigma <- FALSE
  }
  
  sort <- .sortdata(data,inc.nd)
  data <- sort$ty
  t <- data[,1]; y <- log(10)*data[,2] # Convert y from log_10 to natural log for use in model
  if (length(t) ==1) {
    stop("'data' must include at least 2 detected data points")
  }
  dataset <- data.frame(t,y)
  t.nd <- sort$t.nd
  
  if (!inf.sigma) {
    total.pars <- par.no # Calculate total number of parameters
    sigma.save <- sigma
    sigma <- rep(sigma*log(10),length(dataset[,1]))
  } else {
    prior.size <- prior.size + 50
    total.pars <- par.no + 1
  }

  if ((inc.nd) && (is.null(threshold))) {
    stop("'threshold' must be specified when including undetected values")
  } 
  threshold <- threshold*log(10)

  # Perform nested sampling

  #tol <- 0.1 # Set the termination tolerance
  
  # Define transformed priors and log likelihood function
  transformParams <- .generateTransformSingle(dataset,par.no,model,inc.nd,inf.sigma,mumax.prior,mu.mean,mu.sd)
  logllfun <- function(params) {
    return(.logllSingle(params,par.no,modelfunc,model,dataset,inc.nd,threshold,t.nd,inf.sigma,sigma,transformParams,mumax.prior))
  }
  
  # Call the nested sampling routine, which returns the posterior, log evidence & error, logweights and parameter means 
  # and variances
  ret <- nestedSampling(logllfun, total.pars, prior.size, transformParams, exploreFn=ballExplore, tolerance = tol)
  posterior <- ret$posterior
  logevidence <- ret$logevidence
  means <- ret$parameterMeans
  vars <- ret$parameterVariances
  
  # Gather data for posterior plots
  
  # Get equally weighted samples from posterior distribution using staircase sampling
  chosenSamples = getEqualSamples(ret$posterior, n = Inf)
  
  # Fit model with posterior samples
  if (length(t) == 1) {
    tfit <- seq(from=0,by=0.01*t,to=t)
  } else {
    tfit <- seq(from=t[1],by=0.01*t[length(t)],to=t[length(t)])
  }
  posteriorModel = apply(chosenSamples[, -1] , 1, modelfunc, t=tfit)
  
  # Fit model with mean samples
  meanfit <- modelfunc(tfit,means)
  
  # Convert y back to log_10 for output
  
  if (model != "Bar6par") {
    posterior <- cbind(posterior[,1:2],posterior[,3:ncol(posterior)]/log(10))
    means <- means/log(10)
    vars <- vars/log(10)^2
    chosenSamples <- cbind(chosenSamples[,1],chosenSamples[,2:ncol(chosenSamples)]/log(10))
  } else {
    if (inf.sigma == TRUE) {
      posterior <- cbind(posterior[,1:2],posterior[,3:5]/log(10),posterior[,6:7],posterior[,8]*log(10),posterior[,9]/log(10))
      means <- c(means[1:3]/log(10),means[4:5],means[6]*log(10),means[7]/log(10))
      vars <- c(vars[1:3]/log(10)^2,vars[4:5],vars[6]*log(10)^2,vars[7]/log(10)^2)
      chosenSamples <- cbind(chosenSamples[,1],chosenSamples[,2:4]/log(10),chosenSamples[,5:6],chosenSamples[,7]*log(10),
                             chosenSamples[,8]/log(10))
      } else {
      posterior <- cbind(posterior[,1:2],posterior[,3:5]/log(10),posterior[,6:7],posterior[,8]*log(10))
      means <- c(means[1:3]/log(10),means[4:5],means[6]*log(10))
      vars <- c(vars[1:3]/log(10)^2,vars[4:5],vars[6]*log(10)^2)
      chosenSamples <- cbind(chosenSamples[,1],chosenSamples[,2:4]/log(10),chosenSamples[,5:6],chosenSamples[,7]*log(10))
    }
  }
  posteriorModel <- posteriorModel/log(10)
  meanfit <- meanfit/log(10)
  
  # Print out results
  
  cat(rep("#",30),collapse='','\n')
  cat('Model = ', model, '\n')
  cat('mu_max prior type =', mumax.prior, '\n')
  cat(rep("#",30),collapse='',"\n\n")
  cat('log evidence = ', logevidence, '\n')
  cat('Means and standard deviations:', '\n')
  if (model == "linear") {
    cat('Log cell count at time 0, y_0 =',means[1],'+/-',(vars[1])^(1/2),'\n') # Converting back to log_10 scale
    cat('Growth rate, mu_max =',means[2],'+/-',(vars[2])^(1/2),'\n')
  } else if (model == "logistic") {
    cat('Log cell count at time 0, y_0 =',means[1],'+/-',(vars[1])^(1/2),'\n')
    cat('Log final cell count, y_max =',means[2],'+/-',(vars[2])^(1/2),'\n')
    cat('Growth rate, mu_max =',means[3],'+/-',(vars[3])^(1/2),'\n')
  } else if (model == "Bar3par") {
    cat('Log cell count at time 0, y_0 =',means[1],'+/-',(vars[1])^(1/2),'\n')
    cat('Growth rate, mu_max =',means[2],'+/-',(vars[2])^(1/2),'\n')
    cat('h_0 =',means[3],'+/-',(vars[3])^(1/2),'\n')
  } else if (model == "Bar4par") {
    cat('Log cell count at time 0, y_0 =',means[1],'+/-',(vars[1])^(1/2),'\n')
    cat('Log final cell count, y_max =',means[2],'+/-',(vars[2])^(1/2),'\n')
    cat('Growth rate, mu_max =',means[3],'+/-',(vars[3])^(1/2),'\n')
    cat('h_0 =',means[4],'+/-',(vars[4])^(1/2),'\n')
  } else { # model == Bar6par
    cat('Log cell count at time 0, y_0 =',means[1],'+/-',(vars[1])^(1/2),'\n')
    cat('Log final cell count, y_max =',means[2],'+/-',(vars[2])^(1/2),'\n')
    cat('Growth rate, mu_max =',means[3],'+/-',(vars[3])^(1/2),'\n')
    cat('Lag time, lambda =',means[4],'+/-',(vars[4])^(1/2),'\n')
    cat('nu =',means[5],'+/-',(vars[5])^(1/2),'\n')
    cat('m =',means[6],'+/-',(vars[6])^(1/2),'\n')
  }
  if (inf.sigma) {
    cat('Noise level, sigma =',means[par.no+1],'+/-',(vars[par.no+1])^(1/2),'\n\n')
  } else {
    cat('Noise level, sigma = prescribed at',sigma.save,'\n\n')
  }
  
  return(list(posterior=posterior,logevidence=logevidence,means=means,vars=vars,equalposterior=chosenSamples,
              fit.t=tfit,fit.y=posteriorModel,fit.ymean=meanfit))
  ### Returns: 
  ###
  ### posterior: The samples from the posterior, together with their log weights and 
  ###            log likelihoods as a m x n matrix, where m is the
  ###            number of posterior samples and n is the number of
  ###            parameters + 2. The log weights are the first column and the log likelihood values are the
  ###            second column of this matrix. The sum of the log-weights = logZ.
  ###
  ### logevidence: The logarithm of the evidence, a scalar.
  ###
  ### means: A vector of the mean of each parameter, length = no. of parameters.
  ###
  ### vars: A vector of the variance of each parameter, length = no. of parameters.
  ###
  ### equalposterior: Equally weighted posterior samples together with their 
  ###                 log likelihoods as a m x n matrix, where m is the
  ###                 number of posterior samples and n is the number of
  ###                 parameters + 1. The log likelihood values are the
  ###                 first column of this matrix.
  ###
  ### fit.t: A vector of time points at which the model is fitted.
  ###
  ### fit.y: A matrix of fitted model points, y, using posterior parameter samples in the model. Each column represents a different 
  ###        posterior sample.
  ###
  ### fit.ymean: A vector of fitted model points, y, using the mean of the posterior parameter samples in the model.
  
}, ex=function() {
  B092_1.file <- system.file("extdata", "B092_1.csv", package = "babar")
  data <- read.csv(B092_1.file, header=TRUE, sep=",",
                   na.strings=c("ND","NA"))

  # Get a quick approximation of the evidence/model parameters.
  results.linear.short <- Bayesfit(data,model="linear",inf.sigma=FALSE,
                                   tol=10,prior.size=25)

  # Compute a better estimate of the evidence/model parameters (slow so not
  # run as part of the automated examples).
  results.linear.full <- Bayesfit(data,model="linear",inf.sigma=FALSE)
})

