#
#    bayesmeta, an R package for Bayesian random-effects meta-analysis.
#    Copyright (C) 2015  Christian Roever
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#


bayesmeta <- function(y,...)
{
  UseMethod("bayesmeta")
}


bayesmeta.default <- function(y, sigma, labels=names(y),
                              tau.prior="uniform",
                              mu.prior=c("mean"=NA,"sd"=NA),
                              mu.prior.mean=mu.prior[1],
                              mu.prior.sd  =mu.prior[2],
                              delta=0.01, epsilon=0.0001,
                              rel.tol.integrate=2*max(c(50*.Machine$double.eps, 0.5e-28)),
                          abs.tol.integrate=rel.tol.integrate,...)
{
  ptm <- proc.time()
  # some preliminary sanity checks:
  stopifnot(is.vector(y), is.vector(sigma),
            all(is.finite(y)), all(is.finite(sigma)),
            all(sigma>0), length(sigma)==length(y),
            length(mu.prior)==2,
            length(mu.prior.mean)==1, length(mu.prior.sd)==1,
            is.na(mu.prior.mean) || is.finite(mu.prior.mean),
            is.na(mu.prior.sd) || (is.finite(mu.prior.mean) && (mu.prior.sd>0)),
            ((is.na(mu.prior.mean) & is.na(mu.prior.sd))
             || (is.finite(mu.prior.mean) & is.finite(mu.prior.sd))),
            (is.function(tau.prior) | (is.character(tau.prior) && (length(tau.prior)==1))))
  
  k <- length(y)
  if (is.null(labels))
    labels <- as.character(1:k)

  if (is.character(tau.prior)) {
    tau.prior <- match.arg(tolower(tau.prior),
                           c("uniform","jeffreys","shrinkage","dumouchel","bergerdeely"))
    tau.prior <- c("uniform"="uniform", "jeffreys"="Jeffreys",
                   "shrinkage"="shrinkage", "dumouchel"="DuMouchel",
                   "bergerdeely"="BergerDeely")[tau.prior]
    stopifnot(is.element(tau.prior, c("uniform", "Jeffreys", "shrinkage", "DuMouchel", "BergerDeely")))
    if (tau.prior=="uniform") {  # uniform prior on tau:
      pdens <- function(t){d<-rep(1,length(t)); d[t<0]<-0; return(d)}
      attr(pdens, "bayesmeta.label") <- "uniform(min=0, max=Inf)"
    } else if (tau.prior=="Jeffreys") {  # Jeffreys prior:
      pdens <- function(t){return(apply(matrix(t,ncol=1),1,function(x){return(sqrt(sum((x/(sigma^2+x^2))^2)))}))}
      attr(pdens, "bayesmeta.label") <- "Jeffreys prior"
    } else if (tau.prior=="BergerDeely") {  # Berger/Deely prior:
      pdens <- function(t){return(apply(matrix(t,ncol=1),1,
                                        function(x){return(exp(log(x)-sum(log(sigma^2+x^2))/k))}))}
      attr(pdens, "bayesmeta.label") <- "Berger/Deely prior"
    } else {
      # harmonic mean of squared standard errors:
      s02 <- k/sum(1/sigma^2)
      if (tau.prior=="shrinkage") {  # "uniform shrinkage" prior:
        pdens <- function(t){return(apply(matrix(t,ncol=1),1,function(x){return(2*x*s02/(s02+x^2)^2)}))}
        attr(pdens, "bayesmeta.label") <- "uniform shrinkage prior"
      } else if (tau.prior=="DuMouchel") {  # DuMouchel prior:
        s0 <- sqrt(s02)
        pdens <- function(t){return(apply(matrix(t,ncol=1),1,function(x){return(s0/(s0+x)^2)}))}
        attr(pdens, "bayesmeta.label") <- "DuMouchel prior"
      } else warning("could not make sense of 'tau.prior' argument")
    } 
    tau.prior <- pdens
    rm("pdens")
  }
  
  dprior <- function(tau=NA, mu=NA, log=FALSE)
  # prior density (marginal or joint)
  {
    if (all(is.na(tau))) { # marginal density for mu:
      if (is.finite(mu.prior.mean))
        result <- dnorm(mu, mean=mu.prior.mean, sd=mu.prior.sd, log=log)
      else
        result <- rep(ifelse(log, 0, 1), length(mu))
    }
    else if (all(is.na(mu))) { # marginal density for tau:
      if (log) result <- log(apply(matrix(tau,ncol=1),1,tau.prior))
      else result <- apply(matrix(tau,ncol=1),1,tau.prior)
    }
    else if (is.finite(mu.prior.mean) & !all(is.na(mu))) {  # joint density:
      result <- (log(tau.prior(tau))
                 + dnorm(mu, mean=mu.prior.mean, sd=mu.prior.sd, log=TRUE))
      if (!log) result <- exp(result)
    }
    else {  # joint density (uniform on mu):
      result <- tau.prior(tau)
      if (log) result <- log(result)
    }
    return(result)
  }

  likelihood <- function(tau=NA, mu=NA, log=FALSE)
  # likelihood function (marginal or joint)
  {
    if (all(is.na(mu)) & all(is.na(tau))) {
      warning("need to supply at least either 'mu' or 'tau'")
      return()
    }
    else if (all(is.na(tau))) { # return marginal likelihood (numerical):
      warning("marginalization over 'tau' not implemented")
      return()
    }
    else if (all(is.na(mu))) {  # return marginal likelihood (analytical):
      logmarglikeli <- function(t)
      {
        if (is.na(mu.prior.mean)) { # uniform prior
          yTerm   <- y
          varTerm <- t^2 + sigma^2
        }
        else {                      # conjugate normal prior
          yTerm   <- c(mu.prior.mean, y)
          varTerm <- c(mu.prior.sd^2, t^2 + sigma^2)
        }
        conditionalmean <- sum(yTerm/varTerm) / sum(1/varTerm)
        return(-0.5*((length(yTerm)-1) * log(2*pi)
                     + sum(log(varTerm))
                     + sum((yTerm-conditionalmean)^2 / varTerm)
                     + log(sum(1/varTerm))))
      }
      result <- apply(matrix(tau, ncol=1), 1, logmarglikeli)
      result[tau<0] <- -Inf
      if (!log) result <- exp(result)
    }
    else {  # return joint likelihood:
      loglikeli <- function(taumu)
      {
        t2s2 <- taumu[1]^2 + sigma^2
        return(-(k/2)*log(2*pi) - 0.5*sum(log(t2s2)) - 0.5*sum((y-taumu[2])^2 / t2s2))
      }
      result <- apply(cbind(tau,mu), 1, loglikeli)
      result[tau<0] <- -Inf
      if (!log) result <- exp(result)
    }
    return(result)
  }
  
  dposterior <- function(tau=NA, mu=NA, theta=mu, log=FALSE, predict=FALSE)#, individual=FALSE)
  # posterior density (marginal or joint)
  {
    individual <- FALSE
    if (all(is.na(mu))) mu <- theta
    stopifnot(length(individual)==1)
    if (! (is.logical(individual) && (!individual))) {
      indiv.logi <- TRUE
      if (is.numeric(individual))   indiv.which <- which(is.element(1:k, individual))
      if (is.character(individual)) indiv.which <- which(is.element(labels, individual))
      if (length(indiv.which)==0) warning("cannot make sense of 'individual' argument: empty subset.")
    }
    else indiv.logi <- FALSE
    if (predict & indiv.logi) warning("need to specify either 'predict' or 'individual' argument, but not both.")
    if (all(is.na(mu)) & all(is.na(tau))) {
      warning("need to supply at least either 'mu' or 'tau'")
      return()
    }
    else if (all(is.na(tau))) { # return marginal posterior (effect mu, numerical):
      if (all(is.na(support))) {
        warning("'support' not initialized.")
        result <- NA
      }
      else {
        result <- numeric(length(mu))
        if (predict)           # posterior predictive distribution
          for (i in 1:length(mu))
            result[i] <- sum(support[,"weight"]
                             * dnorm(mu[i], mean=support[,"mean"], sd=support[,"sd.pred"]))
        else if (indiv.logi) { # individual-effect distribution
          musigma <-   conditionalmoment(tau=support[,"tau"], individual=indiv.which)
          for (i in 1:length(mu))
            result[i] <- sum(support[,"weight"]
                             * dnorm(mu[i], mean=musigma[,"mean"], sd=musigma[,"sd"]))          
        }
        else                   # (marginal) posterior distribution
          for (i in 1:length(mu))
            result[i] <- sum(support[,"weight"]
                             * dnorm(mu[i], mean=support[,"mean"], sd=support[,"sd"]))
        if (log) result <- log(result)
      }
    }
    else {
      if (predict) warning("'predict' argument ignored!")
      if (indiv.logi) warning("'individual' argument ignored!")
      if (all(is.na(mu))) {  # return marginal posterior (heterogeneity tau, analytical):
        logmargpost <- function(t)
        {
          return(ifelse(t<0, -Inf, log(tau.prior(t)) + likelihood(mu=NA, tau=t, log=TRUE)))
        }
        result <- apply(matrix(tau, ncol=1), 1, logmargpost) - log(integral)
        result[is.nan(result)] <- -Inf
        if (!log) result <- exp(result)
      }
      else {  # return joint posterior:
        loglikeli <- function(taumu)
        {
          t2s2 <- taumu[1]^2 + sigma^2
          return(-(k/2)*log(2*pi) - 0.5*sum(log(t2s2)) - 0.5*sum((y-taumu[2])^2 / t2s2))
        }
        logpost <- function(taumu)
        {
          return(dprior(taumu[1], taumu[2], log=TRUE) + loglikeli(taumu))
        }
        result <- apply(cbind(tau,mu), 1, logpost)
        result[tau<0] <- -Inf
        if (!log) result <- exp(result)
      }
    }
    return(result)
  }

  # compute marginal posterior density's normalizing constant:
  integral <- 1
  integral <- integrate(dposterior, lower=0, upper=Inf,
                        rel.tol=rel.tol.integrate, abs.tol=abs.tol.integrate)$value
  if ((!is.finite(integral)) || (integral <= 0))
    warning("failed integrating marginal posterior (tau)")

  pposterior <- function(tau=NA, mu=NA, theta=mu, predict=FALSE)#, individual=FALSE)
  # posterior cumulative distribution function (CDF) of tau
  {
    individual <- FALSE
    if (all(is.na(mu))) mu <- theta
    stopifnot(length(individual)==1)
    if (! (is.logical(individual) && (!individual))) {
      indiv.logi <- TRUE
      if (is.numeric(individual))   indiv.which <- which(is.element(1:k, individual))
      if (is.character(individual)) indiv.which <- which(is.element(labels, individual))
      if (length(indiv.which)==0) warning("cannot make sense of 'individual' argument: empty subset.")
    }
    else indiv.logi <- FALSE
    if (predict & indiv.logi) warning("need to specify either 'predict' or 'individual' argument, but not both.")
    if (all(is.na(mu))) { # numerical integration for tau
      if (predict) warning("'predict' argument ignored!")
      if (indiv.logi) warning("'individual' argument ignored!")
      cdf <- function(x)
      # marginal CDF of tau
      {
        if (x<=0) p <- 0
        else if (x==Inf) p <- 1
        else p <- integrate(dposterior, lower=0, upper=x,
                            rel.tol=rel.tol.integrate, abs.tol=abs.tol.integrate)$value
        return(p)
      }
      result <- apply(matrix(tau, ncol=1), 1, cdf)
    }
    else if (all(is.na(tau))) { # grid approximation for mu
      if (all(is.na(support))) {
        warning("'support' not initialized.")
        result <- NA
      }
      else if (indiv.logi) { # individual-effect distribution
        #musigma <-   conditionalmoment(tau=support[,"tau"], individual=indiv.which)
        musigma <-   conditionalmoment(tau=support[,"tau"])
        cdf <- function(x)
        {
          p <- sum(pnorm(x, mean=musigma[,"mean"], sd=musigma[,"sd"]) * support[,"weight"])
          return(p)
        }
        result <- apply(matrix(mu, ncol=1), 1, cdf)
      }
      else { # posterior or posterior predictive CDF:
        cdf <- function(x)
        {
          p <- sum(pnorm(x, mean=support[,"mean"], sd=support[,ifelse(predict, "sd.pred", "sd")])
                   * support[,"weight"])
          return(p)
        }
        result <- apply(matrix(mu, ncol=1), 1, cdf)
      }
    }
    else {
      warning("need to supply EITHER tau OR mu, not both.")
      result <- NA
    }
    return(result)
  }

  qposterior <- function(tau.p=NA, mu.p=NA, theta.p=mu.p, predict=FALSE)#, individual=FALSE)
  # posterior quantile function of tau or mu
  {
    individual <- FALSE
    if (all(is.na(mu.p))) mu.p <- theta.p
    stopifnot(length(individual)==1)
    if (! (is.logical(individual) && (!individual))) {
      indiv.logi <- TRUE
      if (is.numeric(individual))   indiv.which <- which(is.element(1:k, individual))
      if (is.character(individual)) indiv.which <- which(is.element(labels, individual))
      if (length(indiv.which)==0) warning("cannot make sense of 'individual' argument: empty subset.")
    }
    else indiv.logi <- FALSE
    if (all(is.na(mu.p))) { # compute tau quantile
      stopifnot(all(tau.p>=0), all(tau.p<=1))
      if (predict) warning("'predict' argument ignored!")
      if (indiv.logi) warning("'individual' argument ignored!")
      upper <- 1
      if (any((tau.p<1) & (tau.p>0))) {
        maxp <- max(tau.p[tau.p<1])
        while (pposterior(upper) < maxp) upper <- upper * 2
      }
      qfun <- function(p)
      {
        stopifnot(p>=0, p<=1)
        if (p==0) quant <- 0
        else {
          if (p==1) quant <- Inf
          else
            quant <- uniroot(function(xx){return(pposterior(tau=xx)-p)},
                             interval=c(0,upper))$root
        }
        return(quant)
      }
      result <- apply(matrix(tau.p,ncol=1),1,qfun)
    }
    else if (all(is.na(tau.p))) { # compute mu quantile
      stopifnot(all(mu.p>=0), all(mu.p<=1))
      if (any((mu.p<1) & (mu.p>0))) {
        minp <- min(c(mu.p[mu.p>0], 1-mu.p[mu.p<1]))
        # derive minimum/maximum based on variance and Chebychev inequality:
        if (predict) {
          lower <- sumstats["mean","theta"] - sqrt(1/minp)*sumstats["sd","theta"]*1.1
          upper <- sumstats["mean","theta"] + sqrt(1/minp)*sumstats["sd","theta"]*1.1
        }
        else if (indiv.logi) {
          lower <- y[indiv.which] - sqrt(1/minp)*sigma[indiv.which]*1.1
          upper <- y[indiv.which] + sqrt(1/minp)*sigma[indiv.which]*1.1
        }
        else {
          lower <- sumstats["mean","mu"] - sqrt(1/minp)*sumstats["sd","mu"]*1.1
          upper <- sumstats["mean","mu"] + sqrt(1/minp)*sumstats["sd","mu"]*1.1
        }
        #while (pposterior(mu=lower,predict=predict,individual=individual) > minp)
        while (pposterior(mu=lower,predict=predict) > minp)
            lower <- lower-sumstats["sd","mu"]
        #while (pposterior(mu=upper,predict=predict,individual=individual) < (1-minp))
        while (pposterior(mu=upper,predict=predict) < (1-minp))
            upper <- upper+sumstats["sd","mu"]
      }
      qfun <- function(p)
      {
        stopifnot(p>=0, p<=1)
        if (p==0) quant <- -Inf
        else {
          if (p==1) quant <- Inf
          else
            #quant <- uniroot(function(yy){return(pposterior(mu=yy,predict=predict,individual=individual)-p)},
            #                 interval=c(lower,upper))$root
            quant <- uniroot(function(yy){return(pposterior(mu=yy,predict=predict)-p)},
                             interval=c(lower,upper))$root
        }
        return(quant)
      }
      result <- apply(matrix(mu.p,ncol=1),1,qfun)
    }
    else {
      warning("need to supply EITHER tau.p OR mu.p, not both.")
      result <- NA
    }
    return(result)
  }

  rposterior <- function(n=1, predict=FALSE)#, individual=FALSE)
  # posterior random number generation for tau and mu
  {
    individual <- FALSE
    stopifnot(n>0, n==round(n), length(individual)==1)
    samp <- matrix(NA, nrow=n, ncol=2, dimnames=list(NULL,c("tau","mu")))
    u <- runif(n=n)
    samp[,"tau"] <- apply(matrix(u,ncol=1), 1, function(x){return(qposterior(tau.p=x))})
    cond.sample <- function(t)
    {
      #cm <- conditionalmoment(t, predict=predict, individual=individual)
      cm <- conditionalmoment(t, predict=predict)
      return(rnorm(n=1, mean=cm[,"mean"], sd=cm[,"sd"]))
    }
    samp[,"mu"] <- apply(matrix(samp[,"tau"],ncol=1), 1, cond.sample)
    return(samp)
  }

  post.interval <- function(tau.level=NA, mu.level=NA, theta.level=mu.level,
                            method=c("shortest","central","evidentiary"),
                            predict=FALSE)#, individual=FALSE)
  # determine credibility interval for tau
  {
    individual <- FALSE
    if (is.na(mu.level)) mu.level <- theta.level
    stopifnot((is.finite(tau.level) & ((tau.level>0) & (tau.level<1)))
              | (is.finite(mu.level) & ((mu.level>0) & (mu.level<1))))
    stopifnot(length(individual)==1)
    method <- match.arg(method)
    if (all(is.na(mu.level))) {  # CI for heterogeneity tau:
      if (predict) warning("'predict' argument ignored!")
      if (! (is.logical(individual) && (!individual))) warning("'individual' argument ignored!")
      if (method=="central")         # central interval
        result <- qposterior(tau.p=c((1-tau.level)/2, 1-(1-tau.level)/2))
      else if (method=="shortest") { # shortest interval
        intwidth <- function(left)
        {
          pleft <- pposterior(tau=left)
          right <- qposterior(tau=tau.level+pleft)
          return(right-left)
        }
        opti <- optimize(intwidth, lower=0, upper=qposterior(tau=1-tau.level))$minimum
        # catch marginal minimum:
        if (intwidth(0) < intwidth(opti))
          result <- c(0, qposterior(tau=tau.level))
        else
          result <- c(opti, qposterior(tau=tau.level+pposterior(tau=opti)))
      }
      else {                         # evidentiary interval
        expectedLogLikeli <- function(left)
        {
          pleft <- pposterior(tau=left)
          right <- qposterior(tau=tau.level+pleft)
          expect <- integrate(function(x){return(likelihood(tau=x,log=TRUE))},
                              lower=left, upper=right,
                              rel.tol=rel.tol.integrate, abs.tol=abs.tol.integrate)
          return(expect$value)
        }
        opti <- optimize(expectedLogLikeli, maximum=TRUE,
                         lower=0, upper=qposterior(tau=1-tau.level))$maximum
        # catch marginal minimum:
        if (expectedLogLikeli(0) < expectedLogLikeli(opti))
          result <- c(0, qposterior(tau=tau.level))
        else
          result <- c(opti, qposterior(tau=tau.level+pposterior(tau=opti)))
          
      }
    }
    else if (all(is.na(tau.level))) {  # CI for effect mu:
      if (method=="central")           # central interval
        #result <- qposterior(mu.p=c((1-mu.level)/2, 1-(1-mu.level)/2), predict=predict, individual=individual)
        result <- qposterior(mu.p=c((1-mu.level)/2, 1-(1-mu.level)/2), predict=predict)
      else if (method=="shortest") {   # shortest interval
        intwidth <- function(left)
        {
          #pleft <- pposterior(mu=left, predict=predict, individual=individual)
          pleft <- pposterior(mu=left, predict=predict)
          #right <- qposterior(mu=mu.level+pleft, predict=predict, individual=individual)
          right <- qposterior(mu=mu.level+pleft, predict=predict)
          return(right-left)
        }
        opti <- optimize(intwidth,
                         #lower=qposterior(mu=(1-mu.level)/50, predict=predict, individual=individual),
                         lower=qposterior(mu=(1-mu.level)/50, predict=predict),
                         #upper=qposterior(mu=1-mu.level, predict=predict, individual=individual))$minimum
                         upper=qposterior(mu=1-mu.level, predict=predict))$minimum
        #result <- c(opti, qposterior(mu=mu.level+pposterior(mu=opti, predict=predict, individual=individual), predict=predict, individual=individual))
        result <- c(opti, qposterior(mu=mu.level+pposterior(mu=opti, predict=predict), predict=predict))
      }
      else {                           # evidentiary interval (not yet implemented)
        warning("sorry, evidentiary credible intervals so far only implemented for heterogeneity parameter tau.")
        result <- c(NA,NA)
      }
    }
    else {
      warning("need to supply EITHER tau.level OR mu.level, not both.")
      result <- c(NA,NA)
    }
    return(result)
  }
  
  conditionalmoment <- function(tau, predict=FALSE, #individual=FALSE,
                                simplify=TRUE)
  # compute conditional posterior moments (mean, sd) of mu for given value(s) of tau
  # if (predict==TRUE), moments of the /posterior predictive distribution/ are returned.
  # if (individual==TRUE), moments of the conditional posterior /of the study-specific effects/
  # (theta[i]) are returned.
  {
    individual <- FALSE
    # interpret the "individual" argument, derive a "logical" and "numeric" component:
    if (all(is.logical(individual))) { #  (logical "individual" argument)
      indiv.logi  <- individual[1]
      if (indiv.logi)
        indiv.which <- 1:k
      else
        indiv.which <- NA          
    }
    else { #  (numerical/character "individual" argument)
      indiv.logi  <- TRUE
      if (all(is.numeric(individual))) # number(s) provided
        indiv.which <- which(is.element(1:k, individual))
      else if (all(is.character(individual))) # label(s) provided
        indiv.which <- which(is.element(labels, individual))
      else warning("cannot make sense of 'individual' argument: funny format.")
      if (length(indiv.which)==0)
        warning("cannot make sense of 'individual' argument: empty subset.")
    }
    if (predict & indiv.logi) warning("need to specify either 'predict' or 'individual' argument, but not both.")
    cm <- function(t)
    {
      if (is.na(mu.prior.mean)) { # uniform prior
        yTerm   <- y
        varTerm <- t^2 + sigma^2
      }
      else {                      # conjugate normal prior
        yTerm   <- c(mu.prior.mean, y)
        varTerm <- c(mu.prior.sd^2, t^2 + sigma^2)
      }
      conditionalvar  <- 1 / sum(1/varTerm)
      conditionalmean <- sum(yTerm/varTerm) * conditionalvar        
      return(c("mean"=conditionalmean,"sd"=sqrt(conditionalvar)))
    }
    # compute conditional moments of Mu
    # (or predictive, if requested):
    if (!indiv.logi) {
      result <- t(apply(matrix(tau,ncol=1),1,cm))
      if (predict) result[,"sd"] <- sqrt(result[,"sd"]^2 + tau^2)
    }
    # compute conditional moments of individual, study-specific means:
    else {
      result <- array(NA, dim=c(length(tau), 2, length(indiv.which)),
                      dimnames=list(NULL, c("mean","sd"), labels[indiv.which]))
      musigma <- t(apply(matrix(tau,ncol=1),1,cm))
      # loop over estimates:
      zero <- (tau==0)
      for (i in 1:length(indiv.which)) {
        if (any(!zero)) {
          result[!zero,"mean",i] <- (y[indiv.which[i]]/sigma[indiv.which[i]]^2 + musigma[!zero,"mean"]/tau[!zero]^2) / (1/sigma[indiv.which[i]]^2+1/tau[!zero]^2)
          result[!zero,"sd",i]   <- sqrt(1 / (1/sigma[indiv.which[i]]^2+1/tau[!zero]^2)
                                         + (musigma[!zero,"sd"] * (1/tau[!zero]^2) / (1/sigma[indiv.which[i]]^2+1/tau[!zero]^2))^2)
        }
        if (any(zero)) {
          result[zero,"mean",i] <- musigma[zero,"mean"]
          result[zero,"sd",i]   <- musigma[zero,"sd"]
        }
      }
      # simplify array to matrix (if possible and desired):
      if ((dim(result)[3] == 1) && (simplify)) result <- result[,,1]
    }
    return(result)
  }

  ISquared <- function(tau)
  # compute heterogeneity measure "I-squared" as a function of tau
  {
    I2 <- rep(NA, length(tau))
    s2 <- (k-1)*sum(1/sigma^2) / (sum(1/sigma^2)^2 - sum(1/sigma^4))
    I2[tau>=0] <- tau[tau>=0]^2 / (tau[tau>=0]^2 + s2)
    return(I2)
  }

  discretize <- function(delta=0.01, epsilon=0.0001, alldivs=TRUE)
  {
    symKL <- function(mean1, sd1, mean2, sd2)
    # (general) symmetrized KL-divergence
    {
      stopifnot(sd1>0, sd2>0)
      return((mean1-mean2)^2 * 0.5 * (1/sd1^2 + 1/sd2^2) + (sd1^2-sd2^2)^2 / (2*sd1^2*sd2^2))
    }
    divergence <- function(tau1, tau2)
    # (symmetrized) divergence b/w two conditionals specified through corresponding tau values
    {
      if (!alldivs) { # evaluate divergence based on (conditional) mu posterior only:
        cm <- conditionalmoment(tau=c(tau1, tau2))
        div <- symKL(cm[1,"mean"], cm[1,"sd"], cm[2,"mean"], cm[2,"sd"])
      }
      else { # evaluate divergence based also on (conditional) individual-study
             # and predictive mu posterior, and determine the maximum divergence:
        # determine array of ALL conditional moments:
        #cm <- array(c(conditionalmoment(tau=c(tau1, tau2)),
        #              conditionalmoment(tau=c(tau1, tau2), predict=TRUE),
        #              conditionalmoment(tau=c(tau1, tau2), individual=TRUE)),
        #            dim=c(2, 2, k+2))
        cm <- array(c(conditionalmoment(tau=c(tau1, tau2)),
                      conditionalmoment(tau=c(tau1, tau2), predict=TRUE)),
                    dim=c(2, 2, k+1))
        # determine individual divergences and their maximum:
        #div <- rep(NA, k+2)
        div <- rep(NA, k+1)
        #for (i in 1:(k+2))
        for (i in 1:(k+1))
          div[i] <- symKL(cm[1,1,i], cm[1,2,i], cm[2,1,i], cm[2,2,i])
        div <- max(div)
      }
      return(div)
    }
    # determine range of interest / upper bound on tau:
    maxtau <- qposterior(tau=1-min(c(epsilon, 1-epsilon)))
    # special procedure for 1st bin
    # (want zero to be first reference point):
    tau <- 0
    # search for upper bin margin:
    upper <- 1
    diverg <- divergence(tau, upper)
    while (diverg < delta) {
      upper <- upper*2
      diverg <- divergence(tau, upper)
    }
    ur <- uniroot(function(t){return(divergence(tau, t)-delta)},
                  lower=tau, upper=upper)
    tau <- ur$root
    prob1 <- 0.0
    prob2 <- pposterior(tau=tau)
    # store result for 1st bin:
    result <- matrix(c(0, prob2-prob1),
                     nrow=1, ncol=2,
                     dimnames=list(NULL,c("tau","weight")))
    # determine following bins (2,...):
    bin <- 2
    while ((tau < maxtau) | (bin<=2)) {  # (at least 2 support points)
      result <- rbind(result, rep(0,2))
      # determine bin's reference point:
      diverg <- divergence(tau, upper)
      while (diverg < delta) {
        upper <- upper*2
        diverg <- divergence(tau, upper)
      }
      ur <- uniroot(function(t){return(divergence(tau, t)-delta)},
                    lower=tau, upper=upper)
      tau <- ur$root
      result[bin,"tau"] <- tau
      # determine bin's upper bound:
      diverg <- divergence(tau, upper)
      while (diverg < delta) {
        upper <- upper*2
        diverg <- divergence(tau, upper)
      }
      ur <- uniroot(function(t){return(divergence(tau, t)-delta)},
                    lower=tau, upper=upper)
      tau <- ur$root
      # determine bin's weight:
      prob1 <- prob2
      prob2 <- pposterior(tau=tau)
      result[bin,"weight"] <- max(c(0.0, prob2-prob1))
      bin <- bin+1
    }
    return(result)
  }

  accelerate <- FALSE  # flag to speed up computations at the cost of a few estimates and a little accuracy

  # compute set of tau support points & corresponding weights:
  support <- discretize(delta=delta, epsilon=epsilon, alldivs=ifelse(accelerate, FALSE, TRUE))
  rm(list=c("discretize"))
  support <- cbind(support,
                   conditionalmoment(support[,"tau"]),
                   "sd.pred"=conditionalmoment(support[,"tau"],predict=TRUE)[,"sd"])
  
  # compute tau posterior's summary statistics:
  sumstats <- matrix(NA, nrow=6, ncol=3,
                     dimnames=list(c("mode", "median", "mean","sd", "95% lower", "95% upper"),
                                   c("tau","mu","theta")))
  expectation <- try(integrate(function(x)return(dposterior(x)*x), lower=0, upper=Inf,
                               rel.tol=rel.tol.integrate, abs.tol=abs.tol.integrate)$value, silent=TRUE)
  if (class(expectation)=="try-error") {
    expectation <- NA
    variance <- NA
  }
  else {
    variance <- try(integrate(function(x)return(dposterior(x)*(x-expectation)^2), lower=0, upper=Inf,
                              rel.tol=rel.tol.integrate, abs.tol=abs.tol.integrate)$value,
                    silent=TRUE)
    if (class(variance)=="try-error")
      variance <- NA
  }
  sumstats[c("mean","sd"),"tau"] <- c(expectation, sqrt(variance))
  sumstats["median","tau"] <- qposterior(0.5)
  sumstats[c("95% lower", "95% upper"),"tau"] <- post.interval(tau=0.95, method=ifelse(accelerate,"central","shortest"))
  if (!accelerate) {
    maxi <- optimize(function(x){return(dposterior(tau=x))},
                     lower=0, upper=qposterior(tau=0.9), maximum=TRUE)
    if ((maxi$maximum <= 2 * .Machine$double.eps^0.25)
        && (dposterior(0) > maxi$objective)) maxi <- list("maximum"=0)
    sumstats["mode","tau"] <- maxi$maximum
    rm(list=c("maxi"))
  }
  rm(list=c("expectation", "variance"))
  
  # compute mu posterior's summary statistics:
  if (all(is.finite(support))) {
    sumstats["mean","mu"] <- sum(support[,"mean"] * support[,"weight"])
    sumstats["sd","mu"] <- sqrt(sum(((support[,"mean"]-sumstats["mean","mu"])^2 + support[,"sd"]^2) * support[,"weight"]))
    sumstats["median","mu"] <- qposterior(mu=0.5)
    sumstats[c("95% lower", "95% upper"),"mu"] <- post.interval(mu=0.95)
    if (!accelerate) {
      maxi <- optimize(function(x){return(dposterior(mu=x))},
                       lower=qposterior(mu=0.1), upper=qposterior(mu=0.9), maximum=TRUE)
      sumstats["mode","mu"] <- maxi$maximum
      rm(list=c("maxi"))
    }
  }

  # compute mu posterior-predictive distribution's summary statistics:
  if (all(is.finite(support))) {
    sumstats["mean","theta"] <- sumstats["mean","mu"]
    sumstats["sd","theta"] <- sqrt(sum(((support[,"mean"]-sumstats["mean","theta"])^2 + support[,"sd.pred"]^2) * support[,"weight"]))
    if (!accelerate) {
      sumstats["median","theta"] <- qposterior(mu=0.5, predict=TRUE)
      sumstats[c("95% lower", "95% upper"),"theta"] <- post.interval(mu=0.95, predict=TRUE)
      maxi <- optimize(function(x){return(dposterior(mu=x, predict=TRUE))},
                       lower=qposterior(mu=0.1, predict=TRUE), upper=qposterior(mu=0.9, predict=TRUE), maximum=TRUE)
      sumstats["mode","theta"] <- maxi$maximum
      rm(list=c("maxi"))
    }
  }

  # compute joint & marginal maximum-likelihood (ML) estimates:
  ml.estimate <- rbind("joint"=c("tau"=NA, "mu"=NA), "marginal"=c("tau"=NA, "mu"=NA))
  map.estimate <- rbind("joint"=c("tau"=NA, "mu"=NA),
                        "marginal"=c("tau"=sumstats["mode","tau"],
                                     "mu"=sumstats["mode","mu"]))
  if (!accelerate) {
    opti <- optim(sumstats["median",c("tau","mu")],
                  function(x){return(-likelihood(tau=x[1], mu=x[2], log=TRUE))})
    ml.estimate["joint",] <- opti$par
    maxi <- optimize(function(x){return(likelihood(tau=x))},
                     lower=0, upper=qposterior(tau=0.99), maximum=TRUE)
    if ((maxi$maximum <= 2 * .Machine$double.eps^0.25)
        && (likelihood(0) > maxi$objective)) maxi <- list("maximum"=0)
    ml.estimate["marginal","tau"] <- maxi$maximum
    ml.estimate["marginal","mu"]  <- NA
    # compute joint maximum-a-posteriori (MAP) estimate:
    opti <- optim(sumstats["median",c("tau","mu")],
                  function(x){return(-dposterior(tau=x[1], mu=x[2], log=TRUE))})
    map.estimate["joint",] <- opti$par
    rm(list=c("opti","maxi"))
  }
  
  ## compute "shrinkage" estimates of theta[i]:
  #shrink <- matrix(NA, nrow=8, ncol=k,
  #                 dimnames=list(c("y","sigma","mode", "median", "mean","sd", "95% lower", "95% upper"),
  #                               labels))
  #shrink["y",] <- y
  #shrink["sigma",] <- sigma
  #if (all(is.finite(support)) & (!accelerate)) {
  #  for (i in 1:k) {
  #    musigma <- conditionalmoment(support[,"tau"], individual=i)
  #    shrink["mean",i] <- sum(musigma[,"mean"] * support[,"weight"])
  #    shrink["sd",i]   <- sqrt(sum(((musigma[,"mean"]-shrink["mean",i])^2 + musigma[,"sd"]^2) * support[,"weight"]))
  #    shrink["median",i] <- qposterior(mu=0.5, individual=i)
  #    shrink[c("95% lower", "95% upper"),i] <- post.interval(mu=0.95, individual=i)
  #    maxi <- optimize(function(x){return(dposterior(mu=x, individual=i))},
  #                     lower=qposterior(mu=0.1), upper=qposterior(mu=0.9), maximum=TRUE)
  #    shrink["mode",i] <- maxi$maximum
  #  }
  #}
  ptm <- proc.time()[1] - ptm[1]
  
  muprior <- c(mu.prior.mean, mu.prior.sd)
  names(muprior) <- c("mean","sd")
  
  # assemble eventual result to be returned:
  result <- list("y"             = y,
                 "sigma"         = sigma,
                 "labels"        = labels,
                 "k"             = k,
                 "tau.prior"     = tau.prior,
                 "mu.prior"      = muprior,
                 "dprior"        = dprior,
                 "likelihood"    = likelihood,
                 "dposterior"    = dposterior,
                 "pposterior"    = pposterior,
                 "qposterior"    = qposterior,
                 "rposterior"    = rposterior,
                 "post.interval" = post.interval,
                 "cond.moment"   = conditionalmoment,
                 "I2"            = ISquared,
                 "summary"       = sumstats,
                 "ML"            = ml.estimate,
                 "MAP"           = map.estimate,
                 #"theta"         = shrink,
                 "support"       = support,
                 "call"          = match.call(expand.dots=FALSE),
                 "init.time"     = c("seconds"=unname(ptm)))
  class(result) <- "bayesmeta"
  return(result)
}


bayesmeta.escalc <- function(y, labels=NULL, ...)
# apply "bayesmeta()" to an "escalc" object
{
  attri <- attributes(y)
  if (all(is.element(c("yi.names", "vi.names"), names(attri)))) { # (for recent "metafor" versions)
    var.names <- c(attri$yi.names, attri$vi.names)
  } else if (is.element("var.names", names(attri))) {             # (for older "metafor" versions)
    var.names <- attri$var.names
  } else {
    stop(paste("Cannont extract \"yi.names\" and \"vi.names\" (or \"var.names\") attribute(s) from escalc object ",
               "(check use of the \"var.names\" option).", sep=""))
  }
  stopifnot(length(var.names)==2, all(is.character(var.names)))
  if (!all(is.element(var.names, names(y)))) {
    stop(paste("Cannont find columns \"",
               var.names[1],"\" and/or \"",
               var.names[2],"\" in escalc object ",
               "(check use of the \"var.names\" option).", sep=""))
  }
  if (is.null(labels)) {
    if (is.element("slab", names(attributes(y[,var.names[1]]))))
      labels <- as.character(attr(y[,var.names[1]], "slab"))
  }
  return(bayesmeta.default(y=as.vector(y[,var.names[1]]),
                           sigma=sqrt(as.vector(y[,var.names[2]])),
                           labels=labels, ...))
}


print.bayesmeta <- function(x,...)
# print a short summary
{
  cat(" 'bayesmeta' object.\n")
  cat(paste("\n",x$k," estimates:\n", sep=""))
  if (length(x$labels)>10)
    cat(paste(c(x$labels[1:10], "..."), collapse=", "))
  else
    cat(paste(x$labels[1:length(x$labels)], collapse=", "))
  cat("\n")
  cat("\ntau prior:\n")
  if (is.element("bayesmeta.label", names(attributes(x$tau.prior))))
    cat(paste(attributes(x$tau.prior)[["bayesmeta.label"]], "\n"))
  else
    print(x$tau.prior)
  cat("\nmu prior:\n")
  if (is.finite(x$mu.prior["mean"]))
    cat(paste("normal(mean=",x$mu.prior["mean"],", sd=",x$mu.prior["sd"],")\n",sep=""))
  else
    cat("uniform(min=-Inf, max=Inf)\n")
  cat("\nML and MAP estimates:\n")
  print(rbind("ML joint"=x$ML["joint",],
              "ML marginal"=x$ML["marginal",],
              "MAP joint"=x$MAP["joint",],
              "MAP marginal"=x$MAP["marginal",]))
  cat("\nmarginal posterior summary:\n")
  print(x$summary[,c("tau","mu")])
  invisible(x)
}


summary.bayesmeta <- function(object,...)
# print a longer summary
{
  cat(" 'bayesmeta' object.\n")
  data <- matrix(c(object$y, object$sigma), ncol=2, dimnames=list(object$labels, c("y","sigma")))
  cat(paste("data (", object$k, " estimates):\n",sep=""))
  if (nrow(data)<=20)
    print(data)
  else {
    print(data[1:20,])
    cat(paste(" [...]  (truncated;", object$k, "estimates total)\n"))
  }
  cat("\ntau prior:\n")
  if (is.element("bayesmeta.label", names(attributes(object$tau.prior))))
    cat(paste(attributes(object$tau.prior)[["bayesmeta.label"]], "\n"))
  else
    print(object$tau.prior)
  cat("\nmu prior:\n")
  if (is.finite(object$mu.prior["mean"]))
    cat(paste("normal(mean=",object$mu.prior["mean"],", sd=",object$mu.prior["sd"],")\n",sep=""))
  else
    cat("uniform(min=-Inf, max=Inf)\n")
  cat("\nML and MAP estimates:\n")
  print(rbind("ML joint"=object$ML["joint",],
              "ML marginal"=object$ML["marginal",],
              "MAP joint"=object$MAP["joint",],
              "MAP marginal"=object$MAP["marginal",]))
  cat("\nmarginal posterior summary:\n")
  print(object$summary)
  cat("\nrelative heterogeneity I^2 (posterior median):", object$I2(tau=object$summary["median","tau"]), "\n")
  invisible(object)
}


plot.bayesmeta <- function(x, main=deparse(substitute(x)),
                           which=1:4, prior=FALSE, bs=FALSE,...)
# generate forest plot and joint and marginal density plots.
{
  q975 <- qnorm(0.975)
  maxtau <- x$qposterior(tau=0.995)*1.1
  murange <- x$qposterior(mu=c(0.005,0.995))
  murange <- murange + c(-1,1)*diff(murange)*0.05

  forestplot <- function(x, main="", battleship=bs)
  {
    # determine x-axis range:
    xrange <- range(c(x$y-q975*x$sigma, x$y+q975*x$sigma,
                      x$summary[c("95% lower", "95% upper"),c("mu","theta")]))
    # empty plot:
    plot(xrange, c(-2, x$k)+c(-1,1)*0.5,
         type="n", axes=FALSE,
         xlab="", ylab="", main=main)
    # add horizontal line dividing data and estimates:
    abline(h=0, col="lightgrey")
    # add vertical "posterior median effect" line:
    abline(v=x$summary["median","mu"], col="black", lty="12")
    if (battleship) {
      ticklength <- 0.45 # (here: maximum height of gaussian)
      maxdens <- c(dnorm(0, sd=x$sigma),
                   x$dposterior(mu=x$summary["mode","mu"]),
                   x$dposterior(theta=x$summary["mode","theta"], predict=TRUE))
      relmaxdens <- maxdens / max(maxdens)
      # density for central 95%
      narg1 <- seq(-q975, q975, le=51)
      ndens1 <- dnorm(narg1) / dnorm(0)
      # density for tails out to +/- 6 sigma
      narg2 <- seq(q975, 6.0, le=40)
      ndens2 <- dnorm(narg2) / dnorm(0)
      
      relsigma <- x$sigma / min(x$sigma)
      for (i in 1:x$k) { # loop over estimates
        # right tail:
        polygon(x$y[i] + c(narg2,rev(narg2))*x$sigma[i],
                x$k-(i-1) + c(ndens2,-rev(ndens2)) * relmaxdens[i] * ticklength,
                col="grey", border=NA)
        # left tail:
        polygon(x$y[i] - c(narg2,rev(narg2))*x$sigma[i],
                x$k-(i-1) + c(ndens2,-rev(ndens2)) * relmaxdens[i] * ticklength,
                col="grey", border=NA)
        # central chunk:
        polygon(x$y[i] + c(narg1,rev(narg1))*x$sigma[i],
                x$k-(i-1) + c(ndens1,-ndens1) * relmaxdens[i] * ticklength,
                col="black", border=NA)
        # central vertical line:
        lines(c(x$y[i], x$y[i]), x$k-(i-1)+ c(-1,1) * relmaxdens[i] * ticklength, col="grey")
      }
    }
    else {
      ticklength <- 0.3
      # draw horizontal lines for individual-study confidence intervals:
      matlines(rbind(x$y-q975*x$sigma, x$y+q975*x$sigma),
               rbind(x$k:1, x$k:1), col=grey(0.3), lty="solid")
      # draw vertical lines for individual-study effect estimates:
      matlines(rbind(x$y, x$y),
               rbind((x$k:1)+ticklength, (x$k:1)-ticklength), col="black", lty="solid")
    }
    if (battleship) {
      # draw blob for mean effect estimate:
      ticklength <- 0.45 # (here: maximum height of gaussian)
      quant <- x$summary[c("95% lower","95% upper"),"mu"]
      quant <- c(quant[1]-2*(x$summary["median","mu"]-quant[1]),
                 quant,
                 quant[2]+2*(quant[2]-x$summary["median","mu"]))
      arg1 <- seq(quant[1], quant[2], le=50)  # left tail
      arg2 <- seq(quant[2], quant[3], le=51)  # central chunk
      arg3 <- seq(quant[3], quant[4], le=50)  # right tail
      dmode <- x$dposterior(mu=x$summary["mode","mu"])
      dens1 <- x$dposterior(mu=arg1) / dmode * relmaxdens[x$k+1]
      dens2 <- x$dposterior(mu=arg2) / dmode * relmaxdens[x$k+1]
      dens3 <- x$dposterior(mu=arg3) / dmode * relmaxdens[x$k+1]
      polygon(c(arg1, rev(arg1)), -1+c(dens1, -rev(dens1))*ticklength, col="grey", border=NA)
      polygon(c(arg3, rev(arg3)), -1+c(dens3, -rev(dens3))*ticklength, col="grey", border=NA)
      polygon(c(arg2, rev(arg2)), -1+c(dens2, -rev(dens2))*ticklength, col="black", border=NA)
      dm <- x$dposterior(mu=x$summary["median","mu"]) / dmode * relmaxdens[x$k+1]
      lines(rep(x$summary["median","mu"],2), -1+c(-1,1)*dm*ticklength, col="grey")
      
      # draw blob for prediction interval:
      quant <- x$summary[c("95% lower","95% upper"),"theta"]
      quant <- c(quant[1]-2*(x$summary["median","theta"]-quant[1]),
                 quant,
                 quant[2]+2*(quant[2]-x$summary["median","theta"]))
      arg1 <- seq(quant[1], quant[2], le=50)
      arg2 <- seq(quant[2], quant[3], le=51)
      arg3 <- seq(quant[3], quant[4], le=50)
      dmode <- x$dposterior(theta=x$summary["mode","theta"], predict=TRUE)
      dens1 <- x$dposterior(theta=arg1, predict=TRUE) / dmode * relmaxdens[x$k+2]
      dens2 <- x$dposterior(theta=arg2, predict=TRUE) / dmode * relmaxdens[x$k+2]
      dens3 <- x$dposterior(theta=arg3, predict=TRUE) / dmode * relmaxdens[x$k+2]
      polygon(c(arg1, rev(arg1)), -2+c(dens1, -rev(dens1))*ticklength, col="grey", border=NA)
      polygon(c(arg3, rev(arg3)), -2+c(dens3, -rev(dens3))*ticklength, col="grey", border=NA)
      polygon(c(arg2, rev(arg2)), -2+c(dens2, -rev(dens2))*ticklength, col="black", border=NA)
      dm <- x$dposterior(theta=x$summary["median","theta"], predict=TRUE) / dmode * relmaxdens[x$k+2]
      lines(rep(x$summary["median","theta"],2), -2+c(-1,1)*dm*ticklength, col="grey")
    }
    else {
      # draw diamond for mean effect estimate:
      ticklength <- 0.4
      polygon(x$summary[c("95% lower", "median", "95% upper", "median"),"mu"],
              rep(-1,4)+c(0,1,0,-1)*ticklength,
              border=NA, col=grey(0.4))
      # draw rectangle for effect prediction interval:
      ticklength <- 0.2
      polygon(x$summary[c("95% lower", "95% lower", "95% upper", "95% upper"),"theta"],
              rep(-2,4)+c(-1,1,1,-1)*ticklength,
              border=NA, col=grey(0.4))
    }
    # add axes & bounding box:
    axis(2, at=x$k:1, labels=x$labels, las=1)
    axis(2, at= c(-1,-2), labels=c(expression(mu), expression(theta[pred.])), las=1)
    axis(1); box()
    invisible()
  }
  
  jointdensity <- function(x, main="")
  {
    # range of tau values:
    tau <- seq(0,maxtau,le=50)
    # range of mu values:
    mu <- seq(murange[1], murange[2], le=50)
    # grid of tau/mu value combinations:
    taumu <- expand.grid(tau,mu)
    # evaluate posterior density at grid points:
    post <- matrix(x$dposterior(tau=taumu[,1], mu=taumu[,2], log=TRUE),
                   nrow=length(tau), ncol=length(mu))
    # determine MAP value:
    map.value <- x$dposterior(x$MAP["joint",1], x$MAP["joint",2], log=TRUE)
    # draw greyscale image:
    image(tau, mu, exp(post), axes=FALSE,
          col=grey((seq(1,0,le=128))^2),
          breaks=seq(0,exp(map.value),length=129),
          xlab="", ylab="", main=main, sub="(joint posterior density)", cex.sub=0.8)
    # add the blue lines for conditional mean & 95% confidence bounds:
    tau2 <- seq(0,maxtau,le=200)
    cm <- x$cond.moment(tau=tau2)
    lines(tau2, cm[,"mean"], col="blue")  
    lines(tau2, cm[,"mean"]-q975*cm[,"sd"], col="blue", lty="dashed")  
    lines(tau2, cm[,"mean"]+q975*cm[,"sd"], col="blue", lty="dashed")
    # add green lines for marginal means & 95% confidence bounds:
    abline(v=x$summary[c("95% lower", "median", "95% upper"),"tau"],
           h=x$summary[c("95% lower", "median", "95% upper"),"mu"],
           col="green2", lty=c("16", "44", "16"))
    # draw ML estimate:
    points(x$ML["joint",1], x$ML["joint",2], col="magenta", pch=4)
    # draw MAP estimate:
    points(x$MAP["joint",1], x$MAP["joint",2], col="red", pch=3)
    # add contour lines:
    contour(tau, mu, post-map.value, add=TRUE, col="red",
            levels=-0.5*qchisq(p=c(0.5, 0.9, 0.95, 0.99), df=2),
            labels=paste(c(50, 90, 95, 99),"%",sep=""))
    # add axes, bounding box, labels, ...
    axis(1); axis(2); box()
    mtext(side=1, line=par("mgp")[1], expression("heterogeneity "*tau))
    mtext(side=2, line=par("mgp")[1], expression("effect "*mu))
    invisible()
  }

  mumarginal <- function(x, main="", priorline=prior)
  {
    # range of mu values:
    mu <- seq(murange[1]-diff(murange)*0.05, murange[2]+diff(murange)*0.05, le=200)
    # corresponding posterior density:
    dens <- x$dposterior(mu=mu)
    # empty plot:
    plot(murange, c(0,max(dens,na.rm=TRUE)), type="n", axes=FALSE,
         xlab="", ylab="", main=main)
    # light grey shaded contour for density across whole range:
    polygon(c(min(mu), mu, max(mu)), c(0,dens,0), border=NA, col=grey(0.90))
    # dark grey shaded contour for density within 95% bounds:
    indi <- ((mu>=x$summary["95% lower","mu"]) & (mu<=x$summary["95% upper","mu"]))
    polygon(c(rep(x$summary["95% lower","mu"],2), mu[indi], rep(x$summary["95% upper","mu"],2)),
            c(0, x$dposterior(mu=x$summary["95% lower","mu"]),
              dens[indi], x$dposterior(mu=x$summary["95% upper","mu"]), 0),
            border=NA, col=grey(0.80))
    # vertical line for posterior median:
    lines(rep(x$summary["median","mu"],2),
          c(0,x$dposterior(mu=x$summary["median","mu"])), col=grey(0.6))
    # actual density line:
    lines(mu, dens, col="black")
    # x-axis:
    abline(h=0, col=grey(0.40))
    # prior density (if requested):
    if (priorline) lines(mu, x$dprior(mu=mu), col="black", lty="dashed")
    # add axes, labels, bounding box, ...
    mtext(side=1, line=par("mgp")[1], expression("effect "*mu))
    mtext(side=2, line=par("mgp")[2], expression("marginal posterior density"))
    axis(1); box()
    invisible()
  }

  taumarginal <- function(x, main="", priorline=prior)
  {
    # range of tau values:
    tau <- seq(0, maxtau*1.1, le=200)
    # corresponding posterior density:
    dens <- x$dposterior(tau=tau)
    # empty plot:
    plot(c(0,maxtau), c(0,max(dens,na.rm=TRUE)),
         type="n", axes=FALSE, xlab="", ylab="", main=main)
    # light grey shaded contour for density across whole range:
    polygon(c(0,tau,max(tau)), c(0,dens,0), border=NA, col=grey(0.90))
    # dark grey shaded contour for density within 95% bounds:
    indi <- ((tau>=x$summary["95% lower","tau"]) & (tau<=x$summary["95% upper","tau"]))
    polygon(c(rep(x$summary["95% lower","tau"],2), tau[indi], rep(x$summary["95% upper","tau"],2)),
            c(0, x$dposterior(tau=x$summary["95% lower","tau"]),
              dens[indi], x$dposterior(tau=x$summary["95% upper","tau"]), 0),
            border=NA, col=grey(0.80))
    # vertical line at posterior median:
    lines(rep(x$summary["median","tau"],2), c(0,x$dposterior(tau=x$summary["median","tau"])), col=grey(0.6))
    # actual density line:
    lines(tau, dens, col="black")
    # x-axis:
    abline(h=0, v=0, col=grey(0.40))
    # prior density (if requested):
    if (priorline) lines(tau, x$dprior(tau=tau), col="black", lty="dashed")
    # add axes, labels, bounding box, ...
    mtext(side=1, line=par("mgp")[1], expression("heterogeneity "*tau))
    mtext(side=2, line=par("mgp")[2], expression("marginal posterior density"))
    axis(1); box()
    invisible()
  }

  # main function:
  stopifnot(all(is.element(which, 1:4)),
            length(which)==length(unique(which)))  
  par.ask <- par("ask")
  for (i in 1:length(which)) {
    if (which[i]==1) forestplot(x, main, battleship=bs)
    else if (which[i]==2) jointdensity(x, main)
    else if (which[i]==3) mumarginal(x, main, priorline=prior)
    else if (which[i]==4) taumarginal(x, main, priorline=prior)
    if (i==1) {
      on.exit(par(ask=par.ask))
      par(ask=TRUE)
    }
  }
  par(ask=par.ask)
  invisible(x)
}


dlomax <- function(x, scale=1, shape=1, log=FALSE)
# probability density function for Lomax distribution
{
  stopifnot(shape>0, scale>0)
  result <- rep(-Inf, length(x))
  result[x>=0] <- log(shape)-log(scale)-(shape+1)*log(1+x[x>=0]/scale)
  if (!log) result <- exp(result)
  return(result)
}


dhalfnormal <- function(x, scale=1, log=FALSE)
# probability density function for half-normal distribution
{
  stopifnot(scale>0)
  result <- rep(-Inf, length(x))
  result[x>=0] <- log(2) + dnorm(x[x>=0], mean=0, sd=scale, log=TRUE)
  if (!log) result <- exp(result)
  return(result)  
}


dhalft <- function(x, scale=1, df, log=FALSE)
# probability density function for half-Student-t distribution
{
  stopifnot(scale>0, df>0)
  result <- rep(-Inf, length(x))
  result[x>=0] <- log(2) - log(scale) + dt(x[x>=0]/scale, df=df, log=TRUE)
  if (!log) result <- exp(result)
  return(result)  
}


dhalfcauchy <- function(x, scale=1, log=FALSE)
# probability density function for half-Cauchy distribution
{
  stopifnot(scale>0)
  return(dhalft(x, scale=scale, df=1, log=log))
}


dinvchisq <- function(x, scale=1, df=1, log=FALSE)
# probability density function for scaled inverse chi-squared distribution
{
  halfdf <- df/2
  result <- rep(-Inf, length(x))
  if (df>0)
    result[x>=0] <- halfdf*log(halfdf) - lgamma(halfdf) + halfdf*log(scale) - (halfdf+1)*log(x[x>=0]) - ((df*scale) / (2*x[x>=0]))
  else if(df==0)    # Jeffreys prior
    result[x>=0] <- -log(x[x>=0])
  else result[x>=0] <- 0 # uniform prior
  if (!log) result <- exp(result)
  return(result)
}


drayleigh <- function(x, scale=1, log=FALSE)
# probability density function for Rayleigh distribution
{
  result <- rep(-Inf, length(x))
  result[x>0] <- log(x[x>0])-2*log(scale)-0.5*(x[x>0]/scale)^2
  if (!log) result <- exp(result)
  return(result)
}


# Turner & al. prior data:
# ========================
# list of 5 possible intervention comparison types:
ctypes <- c("pharmacological vs. placebo / control",
            "pharmacological vs. pharmacological",
            "non-pharmacological vs. placebo / control",
            "non-pharmacological vs. pharmacological",
            "non-pharmacological vs. non-pharmacological")
# list of 16 possible outcome types:
otypes <- c("all-cause mortality",
            "obstetric outcomes",
            "cause-specific mortality / major morbidity event / composite (mortality or morbidity)",
            "resource use / hospital stay / process",
            "surgical / device related success / failure",
            "withdrawals / drop-outs",
            "internal / structure-related outcomes",
            "general physical health indicators",
            "adverse events",
            "infection / onset of new disease",
            "signs / symptoms reflecting continuation / end of condition",
            "pain",
            "quality of life / functioning (dichotomized)",
            "mental health indicators",
            "biological markers (dichotomized)",
            "subjective outcomes (various)")  
# matrix of mu values (see Tab.IV):
meanmat <- matrix(-c(3.95, 3.52, 3.71, 2.34, 2.14, 2.99, 2.71, 2.29, 1.87, 2.49, 2.06, 1.83, 2.54, 2.12, 1.77, 2.70,
                     4.18, 3.75, 3.95, 2.58, 2.37, 3.23, 2.94, 2.53, 2.10, 2.73, 2.29, 2.06, 2.78, 2.35, 2.00, 2.93,
                     4.17, 3.74, 3.93, 2.56, 2.36, 3.21, 2.93, 2.51, 2.10, 2.71, 2.28, 2.05, 2.77, 2.34, 1.99, 2.92,
                     2.92, 2.49, 2.68, 1.31, 1.11, 1.96, 1.67, 1.26, 0.84, 1.46, 1.03, 0.80, 1.51, 1.09, 0.74, 1.67,
                     3.50, 3.08, 3.27, 1.90, 1.69, 2.55, 2.26, 1.85, 1.43, 2.05, 1.61, 1.38, 2.10, 1.67, 1.33, 2.26),
                  nrow=16, ncol=5,
                 dimnames=list(otypes, ctypes))
# matrix of sigma values (see Tab.IV):
sdmat <- matrix(c(1.34, 1.74, 1.74, 1.74, 1.74, 1.74, 1.74, 1.53, 1.52, 1.52, 1.51, 1.52, 1.54, 1.53, 1.52, 1.52,
                  1.41, 1.79, 1.79, 1.79, 1.79, 1.79, 1.79, 1.58, 1.58, 1.58, 1.58, 1.58, 1.60, 1.60, 1.58, 1.58,
                  1.55, 1.91, 1.91, 1.91, 1.91, 1.91, 1.92, 1.72, 1.71, 1.71, 1.71, 1.71, 1.73, 1.72, 1.71, 1.71,
                  1.02, 1.50, 1.51, 1.50, 1.50, 1.51, 1.51, 1.25, 1.24, 1.24, 1.24, 1.25, 1.27, 1.27, 1.24, 1.25,
                  1.26, 1.68, 1.68, 1.68, 1.68, 1.68, 1.68, 1.46, 1.45, 1.45, 1.45, 1.45, 1.47, 1.47, 1.45, 1.45),
                nrow=16, ncol=5,
                dimnames=list(otypes, ctypes))
# array of mu/sigma values:
TurnerEtAlParameters <- array(c(as.vector(meanmat), as.vector(sdmat)),
                              dim=c(16,5,2),
                              dimnames=list("outcome"=otypes, "comparison"=ctypes, "parameter"=c("mu","sigma")))
# remove obsolete objects:
rm(list=c("ctypes", "otypes", "meanmat", "sdmat"))


TurnerEtAlPrior <- function(outcome=c("all-cause mortality",
                                      "obstetric outcomes",
                                      "cause-specific mortality / major morbidity event / composite (mortality or morbidity)",
                                      "resource use / hospital stay / process",
                                      "surgical / device related success / failure",
                                      "withdrawals / drop-outs",
                                      "internal / structure-related outcomes",
                                      "general physical health indicators",
                                      "adverse events",
                                      "infection / onset of new disease",
                                      "signs / symptoms reflecting continuation / end of condition",
                                      "pain",
                                      "quality of life / functioning (dichotomized)",
                                      "mental health indicators",
                                      "biological markers (dichotomized)",
                                      "subjective outcomes (various)"),
                            comparator1=c("pharmacological", "non-pharmacological", "placebo / control"),
                            comparator2=c("pharmacological", "non-pharmacological", "placebo / control"))
#
# Function to return prior parameters, densities, etc. as proposed in
#
#   Turner et al.  Predictive distributions for between-study heterogeneity
#   and simple methods for their application in Bayesian meta-analysis.
#   Statisics in Medicine 4(6):984-998, 2015.
#
# (see Table IV).
#
{
  # match the provided arguments:
  outcome     <- match.arg(outcome)
  comparator1 <- match.arg(comparator1)
  comparator2 <- match.arg(comparator2)
  # list of 3 possible comparators:
  clist <- c("pharmacological", "non-pharmacological", "placebo / control")
  # index matrix to match pairs of comparators to one of five possible scenarios:
  cmatrix <- matrix(c(2,4,1,4,5,3,1,3,NA), nrow=3, ncol=3,
                    dimnames=list(clist, clist))
  # list of 5 possible intervention comparison types:
  ctypes <- dimnames(TurnerEtAlParameters)[["comparison"]]
  # figure out current comparison scenario:
  comparisontype <- ctypes[cmatrix[comparator1, comparator2]]
  # assemble function output:
  mu    <- TurnerEtAlParameters[outcome, comparisontype, "mu"]
  sigma <- TurnerEtAlParameters[outcome, comparisontype, "sigma"]
  dprior <- function(tau, logarithmic=FALSE)
  {
    result <- rep(-Inf, length(tau))
    proper <- is.finite(tau)
    proper[proper] <- (tau[proper]>0)
    result[proper] <- log(2) + log(tau[proper]) + dlnorm(tau[proper]^2, meanlog=mu, sdlog=sigma, log=TRUE)
    if (!logarithmic) result <- exp(result)
    return(result)
  }
  attr(dprior, "bayesmeta.label") <- paste("log-normal(mu=",sprintf("%1.2f",mu),", sigma=",sprintf("%1.2f",sigma),")",sep="")
  result <- list("parameters"      = TurnerEtAlParameters[outcome, comparisontype, ],
                 "outcome.type"    = outcome,
                 "comparison.type" = comparisontype,
                 "dprior"          = dprior,
                 "pprior"          = function(tau) {return(plnorm(tau^2, meanlog=mu, sdlog=sigma))},
                 "qprior"          = function(p) {return(sqrt(qlnorm(p, meanlog=mu, sdlog=sigma)))})
  return(result)
}


# Rhodes & al. prior data:
# ========================
# list of 5 possible intervention comparison types:
ctypes <- c("pharmacological vs. placebo / control",
            "pharmacological vs. pharmacological",
            "non-pharmacological (any)")

# list of 16 possible outcome types:
otypes <- c("obstetric outcome",
            "resource use and hospital stay / process",
            "internal and external structure-related outcome",            
            "general physical health and adverse event and pain and quality of life / functioning",
            "signs / symptoms reflecting continuation / end of condition and infection / onset of new acute / chronic disease",
            "mental health outcome",
            "biological marker",
            "various subjectively measured outcomes")

# matrix of location parameter values (see Tab.3):
locationmat <- matrix(-c(4.13, 2.55, 2.43, 3.16, 3.00, 2.99, 3.41, 2.76,
                         4.40, 2.83, 2.70, 3.44, 3.27, 3.27, 3.68, 3.03,
                         3.99, 2.41, 2.29, 3.02, 2.86, 3.85, 3.27, 2.62),
                      nrow=8, ncol=3,
                      dimnames=list(otypes, ctypes))

# matrix of location parameter values for RESPIRATORY DISEASES (see Tab.A.3.1):
locationmat.r <- matrix(-c(6.03, 4.46, 4.33, 5.07, 4.90, 4.90, 5.31, 4.66,
                           6.31, 4.73, 4.61, 5.34, 5.18, 5.17, 5.59, 4.94,
                           5.89, 4.32, 4.19, 4.93, 4.76, 4.76, 5.17, 4.52),
                        nrow=8, ncol=3,
                        dimnames=list(otypes, ctypes))

# matrix of location parameter values for CANCER (see Tab.A.3.2):
locationmat.c <- matrix(-c(1.57,-0.01, 0.13, 0.60, 0.44, 0.43, 0.85, 0.20,
                           1.85, 0.27, 0.14, 0.88, 0.71, 0.71, 1.13, 0.48,
                           1.43,-0.15,-0.27, 0.46, 0.30, 0.29, 0.71, 0.06),
                        nrow=8, ncol=3,
                        dimnames=list(otypes, ctypes))

# matrix of scale parameter values (see Tab.3):
scalemat <- matrix(c(2.34, 2.73, 2.50, 2.50, 2.50, 2.16, 2.83, 2.58,
                     2.31, 2.70, 2.46, 2.44, 2.47, 2.14, 2.78, 2.59,
                     2.11, 2.57, 2.32, 2.27, 2.33, 1.93, 2.66, 2.41),
                   nrow=8, ncol=3,
                   dimnames=list(otypes, ctypes))

# matrix of scale parameter values for RESPIRATORY DISEASES (see Tab.A.3.1):
scalemat.r <- matrix(c(2.36, 2.74, 2.51, 2.51, 2.50, 2.17, 2.83, 2.59,
                       2.31, 2.70, 2.46, 2.45, 2.47, 2.14, 2.78, 2.59,
                       2.21, 2.57, 2.33, 2.28, 2.33, 1.94, 2.66, 2.41),
                     nrow=8, ncol=3,
                     dimnames=list(otypes, ctypes))

# matrix of scale parameter values for CANCER (see Tab.A.3.2):
scalemat.c <- matrix(c(2.45, 2.83, 2.61, 2.61, 2.60, 2.28, 2.93, 2.68,
                       2.41, 2.79, 2.56, 2.55, 2.57, 2.25, 2.87, 2.68,
                       2.24, 2.68, 2.45, 2.40, 2.46, 2.08, 2.78, 2.53),
                     nrow=8, ncol=3,
                     dimnames=list(otypes, ctypes))

# array of mu/sigma values:
RhodesEtAlParameters <- array(c(as.vector(locationmat), as.vector(scalemat),
                                as.vector(locationmat.r), as.vector(scalemat.r),
                                as.vector(locationmat.c), as.vector(scalemat.c)),
                              dim=c(8,3,2,3),
                              dimnames=list("outcome"=otypes, "comparison"=ctypes,
                                            "parameter"=c("location","scale"),
                                            "medical area"=c("other","respiratory","cancer")))
# remove obsolete objects:
rm(list=c("ctypes", "otypes", "locationmat", "scalemat", "locationmat.r", "scalemat.r", "locationmat.c", "scalemat.c"))

RhodesEtAlPrior <- function(outcome=c(NA,
                                      "obstetric outcome",
                                      "resource use and hospital stay / process",
                                      "internal and external structure-related outcome",            
                                      "general physical health and adverse event and pain and quality of life / functioning",
                                      paste("signs / symptoms reflecting continuation / end of condition and infection",
                                            "/ onset of new acute / chronic disease"),
                                      "mental health outcome",
                                      "biological marker",
                                      "various subjectively measured outcomes"),
                            comparator1=c("pharmacological", "non-pharmacological", "placebo / control"),
                            comparator2=c("pharmacological", "non-pharmacological", "placebo / control"),
                            area=c("other","respiratory","cancer"))
#
# Function to return prior parameters as proposed in
#
#   Rhodes et al.
#   Predictive distributions were developed for the extent of heterogeneity in meta-analyses of continuous outcome data.
#   Journal of Clinical Epidemiology 68(1):52-60, 2015.
#
# (see Table 3).
#
# Technically, this function mostly retrieves the parameters from a pre-defined array,
# the 3-dimensional "RhodesEtAlParameters" array. (This is more convenient than having
# to deal with the array itself, mostly due to partial argument matching, etc.)
#
{
  # match the provided arguments:
  outcome     <- match.arg(outcome)
  if (!is.na(outcome)) {
    comparator1 <- match.arg(comparator1)
    comparator2 <- match.arg(comparator2)
    area        <- match.arg(area)
    # list of 3 possible comparators:
    clist <- c("pharmacological", "non-pharmacological", "placebo / control")
    # index matrix to match pairs of comparators to one of five possible scenarios:
    cmatrix <- matrix(c(2,3,1,3,3,3,1,3,NA), nrow=3, ncol=3,
                      dimnames=list(clist, clist))
    # list of 5 possible intervention comparison types:
    ctypes <- dimnames(RhodesEtAlParameters)[["comparison"]]
    # figure out current comparison scenario:
    comparisontype <- ctypes[cmatrix[comparator1, comparator2]]
    # assemble function output:
    location <- RhodesEtAlParameters[outcome, comparisontype, "location", area]
    scale    <- RhodesEtAlParameters[outcome, comparisontype, "scale", area]
  } else {
    outcome <- comparisontype <- area <- "any"
    location <- -3.44
    scale <- 2.59
  }
  dlt5 <- function(x, location=0, scale=1)
  # log-density of log-t distribution with 5 d.f.
  {
    stopifnot(is.finite(location), is.finite(scale), scale>0)
    return(dt((log(x)-location)/scale, df=5, log=TRUE) - log(scale) - log(x))
  }
  dprior <- function(tau, logarithmic=FALSE)
  {
    result <- rep(-Inf, length(tau))
    proper <- is.finite(tau)
    proper[proper] <- (tau[proper]>0)
    #result[proper] <- log(2) + log(tau[proper]) + dlnorm(tau[proper]^2, meanlog=mu, sdlog=sigma, log=TRUE)
    result[proper] <- log(2) + log(tau[proper]) + dlt5(tau[proper]^2, location=location, scale=scale)
    if (!logarithmic) result <- exp(result)
    return(result)
  }
  attr(dprior, "bayesmeta.label") <- paste("log-Student-t(location=",sprintf("%1.2f",location),", scale=",sprintf("%1.2f",scale),", d.f.=5)",sep="")
  plt5 <- function(x, location=0, scale=1)
  # cumulative distribution function (CDF) of log-t distribution with 5 d.f.
  {
    stopifnot(is.finite(location), is.finite(scale), scale>0)
    return(pt((log(x)-location)/scale, df=5))
  }
  qlt5 <- function(p, location=0, scale=1)
  # quantile function (inverse CDF) of log-t distribution with 5 d.f.
  {
    stopifnot(is.finite(location), is.finite(scale), scale>0)
    return(exp((qt(p, df=5)*scale)+location))
  }
  #result <- list("parameters"      = RhodesEtAlParameters[outcome, comparisontype, , area],
  result <- list("parameters"      = c("location"=location, "scale"=scale),
                 "outcome.type"    = outcome,
                 "comparison.type" = comparisontype,
                 "medical.area"    = area,
                 "dprior"          = dprior,
                 "pprior"          = function(tau) {return(plt5(tau^2, location=location, scale=scale))},
                 "qprior"          = function(p) {return(sqrt(qlt5(p, location=location, scale=scale)))})
  return(result)
}


forest.bayesmeta <- function(x, xlab="effect size", refline=0, cex=1,...)
# forest plot for a "bayesmeta" object
# based on the "metafor" package's plotting functions
{
  if (!requireNamespace("metafor"))
    stop("required 'metafor' package not available!")
  metafor::forest.default(x=x$y, sei=x$sigma,
                          showweight=FALSE,  # (IV-weights don't make sense here)
                          ylim=c(-x$k-4, 1),
                          level=95,          # (95% level is intentionally hard-coded)
                          refline=refline,
                          xlab=xlab,
                          slab=x$labels,
                          rows=seq(-2, -x$k - 1, by = -1),
                          cex=cex, ...)
  metafor::addpoly(x$summary["median","mu"], ci.lb=x$summary["95% lower","mu"], ci.ub=x$summary["95% upper","mu"],
                   rows = -x$k-2.5, mlab=expression("mean effect ("*mu*")"), level=95, cex=cex, ...)
  metafor::addpoly(x$summary["median","theta"], ci.lb=x$summary["95% lower","theta"], ci.ub=x$summary["95% upper","theta"],
                   rows = -x$k-3.5, mlab=expression("prediction ("*vartheta[k+1]*")"), level=95, cex=cex, ...)
  
  plotdata <- cbind("95% lower"=x$y-qnorm(0.975)*x$sigma, "estimate"=x$y, "95% upper"=x$y+qnorm(0.975)*x$sigma)
  plotdata <- rbind(plotdata, t(x$summary[c("95% lower","median","95% upper"),c("mu","theta")]))
  rownames(plotdata) <- c(x$labels, c("mean effect(mu)", "prediction (theta)"))
  invisible(plotdata)
}
