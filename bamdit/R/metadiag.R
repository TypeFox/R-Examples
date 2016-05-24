#' Bayesian Meta-Analysis of diagnostic test data
#'
#' This function performers a Bayesian meta-analysis of diagnostic test data by
#' fitting a bivariate random effects model. The number of true positives and
#' false positives are modeled with two conditional Binomial distributions and
#' the random-effects are based on a bivariate scale mixture of Normals.
#' Computations are done by calling JAGS (Just Another Gibbs Sampler) to perform
#' MCMC (Markov Chain Monte Carlo) sampling and returning an object of the
#' class \emph{mcmc.list}.
#'
#'
#' @param data A data frame with at least 4 columns containing the true positives (tp),
#' number of patients with disease (n1), false positives (fp), number of patients without
#' disease (n2)
#' @param re Random effects distribution for the resulting model. Possible
#' values are \emph{normal} for bivariate random effects and \emph{sm} for scale mixtures
#' @param link The link function used in the model. Possible values are
#' \emph{logit}, \emph{cloglog} \emph{probit}.
#'
#' @param mean.mu.D prior mean of D, default value is 0
#' @param mean.mu.S prior mean of S, default value is 0
#' @param sd.mu.D prior standard deviation of D, default value is 10
#' @param sd.mu.S prior standard deviation of S, default value is 10
#' @param sigma.D.upper upper bound of the uniform prior of sigma.S default value is 10
#' @param sigma.S.upper upper bound of the uniform prior of sigma.S default value is 10
#' @param mean.Fisher.rho mean of rho in the Fisher scale default value is 0
#' @param sd.Fisher.rho   standard deviation of rho in the Fisher scale default value is 1/sqrt(2)
#' @param df              Degrees of freedom for the scale mixture distribution default value is 4
#' @param split.w         Split the w parameter in two indpendentent weights one for each random effect. The defualt value is FALSE.
#' @param n.1.new Number of patients with desease in a predictive study default is 50
#' @param n.2.new Number of patients with non-desease in a predictive study default is 50
#' @param nr.chains       Number of chains for the MCMC computations default 5
#' @param nr.iterations   Number of iterations after adapting the MCMC default is 1000
#' @param nr.adapt        Number of iterations in the adaptation process defualt is 500
#' @param nr.burnin       Number of interation descared for burnin period default is 500
#' @param nr.thin         Thinning rate, it must be a positive integer the defulat value 1
#' @param be.quiet        Do not print warning message if the model does not adapt default value is FALSE
#' @param r2jags          Which interface is used to link R to JAGS (rjags and R2jags) default value is R2Jags TRUE
#' @keywords file
#' @examples
#'
#' ## execute analysis
#' \dontrun{
#' data(mri)
#' mri
#' md <- metadiag(mri)
#' summary(md)
#' }
#'
#' @import R2jags
#' @import rjags
#' @export
#-----
metadiag <- function(
          # Data
           data,
          # Arguments for the model:
          re              = "normal",
          link            = "logit",

					# Hyperpriors parameters............................................
					mean.mu.D       = 0,
					mean.mu.S       = 0,
					sd.mu.D         = 10,
					sd.mu.S         = 10,
					sigma.D.upper   = 10,
					sigma.S.upper   = 10,
					mean.Fisher.rho = 0,
					sd.Fisher.rho   = 1/sqrt(2),
					df              = 4,

          # Split weights
          split.w         = FALSE,

          # Predictions
          n.1.new         = 50,
          n.2.new         = 50,

          # MCMC setup........................................................
					nr.chains       = 2,
					nr.iterations   = 5000,
					nr.adapt        = 1000,
					nr.burnin       = 1000,
          nr.thin         = 1,

          # Further options to link jags and R ...............................
					be.quiet        = FALSE,
          r2jags          = TRUE
          )
{


# Model errors checking-----

if(re=="normal" & split.w==TRUE)stop("Normal random effects and splitting weights are non compatible options")

re.test <- re %in% c("normal", "sm")

if(!re.test)stop("This random effects distribution is not implemented")

link.test <- link %in% c("logit", "cloglog", "probit")

if(!link.test)stop("This link function is not implemented")

	# Setting up hyperparameters ...
	      pre.mu.D <- 1/(sd.mu.D*sd.mu.D)
	      pre.mu.S <- 1/(sd.mu.S*sd.mu.S)
	pre.Fisher.rho <- 1/(sd.Fisher.rho * sd.Fisher.rho)

	# Test if the data is compatible ...
	N <- dim(data)[1]

  # Setting up data nodes ...
	tp <- data[,1]
	fp <- data[,3]
	n1 <- data[,2]
	n2 <- data[,4]

  # Data errors
  if(tp>n1 || fp>n2)stop("the data is inconsistent")
  if(missing(data))stop("NAs are not alow in this function")

#-----
# Data, initial values and parameters .....................................................
	data.model <-
            list(N = N,
	              tp = tp,
	              fp = fp,
	              n1 = n1,
	              n2 = n2,
	       mean.mu.D = mean.mu.D,
	       mean.mu.S = mean.mu.S,
	        pre.mu.D = pre.mu.D,
	        pre.mu.S = pre.mu.S,
	   sigma.D.upper = sigma.D.upper,
	   sigma.S.upper = sigma.S.upper,
	 mean.Fisher.rho = mean.Fisher.rho,
	  pre.Fisher.rho = pre.Fisher.rho,
	 # Predictions
	 n.1.new         = 50,
	 n.2.new         = 50)

  # Parameters to monitor ....................................................................
	parameters.model <-
                c("se.pool",
	                "sp.pool",
                   "se.new",
                   "sp.new",
                   "tp.new",
                   "fp.new",
	                   "mu.D",
                     "mu.S",
                  "sigma.D",
	                "sigma.S",
	                    "rho")

	# This take the weights for the scale mixture random effects model.

        if(re == "sm" & split.w == TRUE)
          parameters.model <- c(parameters.model[], "w1", "w2")
        else
          if(re=="sm") parameters.model <- c(parameters.model, "w")

        if(re=="sm") data.model$df <- df

# Model construction

inits.model <- list(mu.D = 0,
                    mu.S = 0,
                 sigma.D = 1,
                 sigma.S = 1,
                 z       = 0)

# Model BUGS script

#----
blueprint <- function(link = "logit", re = "normal", split.w = FALSE)
{
  if(split.w == TRUE) re <- "sm.split"


  #----
  # Block for data model ......................................................................
  dm <-
    "
  model
{
  for(i in 1:N)
{
  tp[i] ~ dbin(TPR[i], n1[i])
  fp[i] ~ dbin(FPR[i], n2[i])
  "

  #----

  #----
  # Block for the link function................................................................
  link.logit <-
    "
  logit(TPR[i]) <- (D[i] + S[i])/2
  logit(FPR[i]) <- (S[i] - D[i])/2
  "


  link.cloglog <-
    "
  cloglog(TPR[i]) <- (D[i] + S[i])/2
  cloglog(FPR[i]) <- (S[i] - D[i])/2
  "

  link.probit <-
    "
  probit(TPR[i]) <- (D[i] + S[i])/2
  probit(FPR[i]) <- (S[i] - D[i])/2
  "

  #----

  # Block for structural distribution .........................................................
  re.normal <-
    "
  S[i] ~ dnorm(mu.S, pre.S)
  D[i] ~ dnorm(mu.D.S[i], pre.D.S)
  mu.D.S[i] <- mu.D + rho * sigma.D / sigma.S * (S[i] - mu.S)

}

  # Hyper priors
  mu.D ~  dnorm(mean.mu.D, pre.mu.D)
  mu.S ~  dnorm(mean.mu.S, pre.mu.S)

  # Dispersion parameters
  sigma.D ~ dunif(0, sigma.D.upper)
  sigma.S ~ dunif(0, sigma.S.upper)
  pre.D <- 1/(sigma.D*sigma.D)
  pre.S <- 1/(sigma.S*sigma.S)

  #Conditional precision
  pre.D.S <- pre.D/(1-rho*rho)

  # Correlation
  z ~ dnorm(mean.Fisher.rho, pre.Fisher.rho)
  rho <- 2*exp(z)/(1+exp(z)) - 1

  # Predictions ...
  mu.hat[1] <- mu.D
  mu.hat[2] <- mu.S

  Sigma.hat[1, 1] <- pow(sigma.D, 2)
  Sigma.hat[2, 2] <- pow(sigma.S, 2)
  Sigma.hat[1, 2] <- rho * sigma.D * sigma.S
  Sigma.hat[2, 1] <- Sigma.hat[1, 2]

  Omega.hat[1:2, 1:2] <- inverse(Sigma.hat[1:2, 1:2])
  DSnew[1:2] ~ dmnorm(mu.hat[1:2], Omega.hat[1:2 ,1:2])
  "



  #----
  re.sm <-
    "
  S[i] ~ dnorm(mu.S, pre.w.S[i])
  D[i] ~ dnorm(mu.D.S[i], pre.w.D.S[i])
  mu.D.S[i] <- mu.D + rho * sigma.D / sigma.S * (S[i] - mu.S)

  lambda[i] ~ dchisqr(df)
  w[i] <- df/lambda[i]

  pre.w.D[i] <- pre.D / w[i]
  pre.w.S[i] <- pre.S / w[i]

  #Conditional precision
  pre.w.D.S[i] <- pre.w.D[i] / (1-rho*rho)
}

  # Hyper priors
  mu.D ~  dnorm(mean.mu.D, pre.mu.D)
  mu.S ~  dnorm(mean.mu.S, pre.mu.S)

  # Dispersion parameters
  sigma.D ~ dunif(0, sigma.D.upper)
  sigma.S ~ dunif(0, sigma.S.upper)
  pre.D <- 1/(sigma.D*sigma.D)
  pre.S <- 1/(sigma.S*sigma.S)


  # Correlation
  z ~ dnorm(mean.Fisher.rho, pre.Fisher.rho)
  rho <- 2*exp(z)/(1+exp(z)) - 1

  # Predictions ...
  mu.hat[1] <- mu.D
  mu.hat[2] <- mu.S

  Sigma.hat[1, 1] <- pow(sigma.D, 2)
  Sigma.hat[2, 2] <- pow(sigma.S, 2)
  Sigma.hat[1, 2] <- rho * sigma.D * sigma.S
  Sigma.hat[2, 1] <- Sigma.hat[1, 2]

  Omega.hat[1:2,1:2] <- inverse(Sigma.hat[1:2, 1:2])
  DSnew[1:2] ~ dmt(mu.hat[1:2], Omega.hat[1:2 ,1:2], df)
  "

  re.sm.split <-
    "
  S[i] ~ dnorm(mu.S, pre.w.S[i])
  D[i] ~ dnorm(mu.D.S[i], pre.w.D.S[i])
  mu.D.S[i] <- mu.D + rho * sigma.D / sigma.S * (S[i] - mu.S)

  lambda1[i] ~ dchisqr(df)
  w1[i] <- df/lambda1[i]

  lambda2[i] ~ dchisqr(df)
  w2[i] <- df/lambda2[i]

  pre.w.D[i] <- pre.D / w1[i]
  pre.w.S[i] <- pre.S / w2[i]

  #Conditional precision
  pre.w.D.S[i] <- pre.w.D[i] / (1-rho*rho)
}

  # Hyper priors
  mu.D ~  dnorm(mean.mu.D, pre.mu.D)
  mu.S ~  dnorm(mean.mu.S, pre.mu.S)

  # Dispersion parameters
  sigma.D ~ dunif(0, sigma.D.upper)
  sigma.S ~ dunif(0, sigma.S.upper)
  pre.D <- 1/(sigma.D*sigma.D)
  pre.S <- 1/(sigma.S*sigma.S)


  # Correlation
  z ~ dnorm(mean.Fisher.rho, pre.Fisher.rho)
  rho <- 2*exp(z)/(1+exp(z)) - 1

  # Predictions ...
  mu.hat[1] <- mu.D
  mu.hat[2] <- mu.S

  Sigma.hat[1, 1] <- pow(sigma.D, 2)
  Sigma.hat[2, 2] <- pow(sigma.S, 2)
  Sigma.hat[1, 2] <- rho * sigma.D * sigma.S
  Sigma.hat[2, 1] <- Sigma.hat[1, 2]

  Omega.hat[1:2,1:2] <- inverse(Sigma.hat[1:2, 1:2])
  DSnew[1:2] ~ dmt(mu.hat[1:2], Omega.hat[1:2 ,1:2], df)
  "



  # Block of parameters of interest depending on the links .................................
  par.logit <- "
  # Parameters of interest
  # Pooled summaries ...
  x <- (mu.D + mu.S)/2
  y <- (mu.S - mu.D)/2
  se.pool <- ilogit(x)
  sp.pool <- 1 - ilogit(y)

  # Predictive summaries ...
  x.new <- (DSnew[1] + DSnew[2])/2
  y.new <- (DSnew[2] - DSnew[1])/2
  se.new <- ilogit(x.new)
  fpr.new <- ilogit(y.new)
  sp.new <- 1 - fpr.new

  tp.new ~ dbin( se.new, n.1.new)
  fp.new ~ dbin(fpr.new, n.2.new)
}"


  par.cloglog <- "
  # Parameters of interest
  # Pooled summaries ...
  x <- (mu.D + mu.S)/2
  y <- (mu.S - mu.D)/2
  se.pool <- icloglog(x)
  sp.pool <- 1 - icloglog(y)

  # Predictive summaries ...
  x.new <- (DSnew[1] + DSnew[2])/2
  y.new <- (DSnew[2] - DSnew[1])/2
  se.new <- icloglog(x.new)
  fpr.new <- icloglog(y.new)
  sp.new <- 1 - fpr.new

  tp.new ~ dbin( se.new, n.1.new)
  fp.new ~ dbin(fpr.new, n.2.new)
}"

  par.probit <-
  "
  # Parameters of interest
  # Pooled summaries ...
  x <- (mu.D + mu.S)/2
  y <- (mu.S - mu.D)/2
  se.pool <- phi(x)
  sp.pool <- 1 - phi(y)

  # Predictive summaries ...
  x.new <- (DSnew[1] + DSnew[2])/2
  y.new <- (DSnew[2] - DSnew[1])/2
  se.new <- phi(x.new)
  fpr.new <- phi(y.new)
  sp.new <- 1 - fpr.new

  tp.new ~ dbin(se.new, n.1.new)
  fp.new ~ dbin(fpr.new, n.2.new)
  }
  "

#----
# possible models
#normal random effects
m1 <- paste(dm, link.logit,   re.normal, par.logit)
m2 <- paste(dm, link.cloglog, re.normal, par.cloglog)
m3 <- paste(dm, link.probit,  re.normal, par.probit)

# sm random effects
m4 <- paste(dm, link.logit,   re.sm, par.logit)
m5 <- paste(dm, link.cloglog, re.sm, par.cloglog)
m6 <- paste(dm, link.probit,  re.sm, par.probit)


# sm random effects with two t-distributions one for each random effect
m7 <- paste(dm, link.logit,   re.sm.split,  par.logit)
m8 <- paste(dm, link.cloglog, re.sm.split,  par.cloglog)
m9 <- paste(dm, link.probit,  re.sm.split,  par.probit)


#----
  switch(re,
         normal = switch(link,
                           logit = return(m1),
                         cloglog = return(m2),
                          probit = return(m3),
                         stop("The model you requested is not implemented.")
                         ),
             sm = switch(link,
                           logit = return(m4),
                         cloglog = return(m5),
                          probit = return(m6),
                         stop("The model you requested is not implemented.")
                         ),

         sm.split = switch(link,
                           logit = return(m7),
                         cloglog = return(m8),
                          probit = return(m9),
                         stop("The model you requested is not implemented.")
         ),
        stop("The model you requested is not implemented.")
        )

}

model.bugs <- blueprint(link, re, split.w)


if(r2jags == TRUE){
  # Use R2jags as interface for JAGS ...
  results <- jags(              data = data.model,
                  parameters.to.save = parameters.model,
                            #   inits = inits.model,
                          model.file = textConnection(model.bugs),
                            n.chains = nr.chains,
                              n.iter = nr.iterations
                       )
  }
  else {
  # Use rjags as interface for JAGS ...
  # Send the model to JAGS, check syntax, run ...
	jm <- jags.model(file     = textConnection(model.bugs),
	                 data     = data.model,
                   inits    = inits.model,
	                 n.chains = nr.chains,
	                 n.adapt  = nr.adapt,
	                 quiet    = be.quiet)

	results <- coda.samples(jm,
	                        variable.names = parameters.model,
	                        n.iter         = nr.iterations)
  }

if(r2jags == FALSE)
  {cat("You are using the package rjags as interface to JAGS.", "\n")
   cat("The plot functions for output analysis are not implemented in this bamdit version", "\n")
  }
return(results)
}


