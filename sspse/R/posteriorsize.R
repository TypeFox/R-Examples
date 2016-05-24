#' Estimating hidden population size using RDS data
#' 
#' \code{\link{posteriorsize}} computes the posterior distribution of the
#' population size based on data collected by Respondent Driven Sampling. The
#' approach approximates the RDS via the Sequential Sampling model of Gile
#' (2008). As such, it is referred to as the Sequential Sampling - Population Size Estimate (SS-PSE).
#' It uses the order of selection of the sample to provide information
#' on the distribution of network sizes over the population members.
#' 
#' @param s vector of integers; the vector of degrees from the RDS in order
#' they are recorded.
#' @param median.prior.size scalar; A hyperparameter being the mode of the
#' prior distribution on the population size.
#' @param interval count; the number of proposals between sampled statistics.
#' @param burnin count; the number of proposals before any MCMC sampling is
#' done. It typically is set to a fairly large number.
#' @param maxN integer; maximum possible population size. By default this is
#' determined from an upper quantile of the prior distribution.
#' @param K count; the maximum degree for an individual. This is usually
#' calculated as \code{round(quantile(s,0.80))}.
#' @param samplesize count; the number of Monte-Carlo samples to draw to
#' compute the posterior. This is the number returned by the
#' Metropolis-Hastings algorithm.The default is 1000.
#' @param quartiles.prior.size vector of length 2; A pair of hyperparameters
#' being the lower and upper quartiles of the prior distribution on the
#' population size. For example, \cr \code{quartiles.prior.size=c(1000,4000)}
#' corresponds to a prior where the lower quartile (25\%) is 1000 and the upper
#' (75\%) is 4000.
#' @param mean.prior.size scalar; A hyperparameter being the mean of the prior
#' distribution on the population size.
#' @param mode.prior.size scalar; A hyperparameter being the mode of the prior
#' distribution on the population size.
#' @param priorsizedistribution character; the type of parametric distribution
#' to use for the prior on population size. The options are \code{beta} (for a
#' Beta prior on the sample proportion (i.e. \eqn{n/N})), \code{flat}
#' (uniform), \code{nbinom} (Negative-Binomial), and \code{pln}
#' (Poisson-log-normal). The default is \code{beta}.
#' @param effective.prior.df scalar; A hyperparameter being the effective
#' number of samples worth of information represented in the prior distribution
#' on the population size. By default this is 1, but it can be greater (or
#' less!) to allow for different levels of uncertainty.
#' @param sd.prior.size scalar; A hyperparameter being the standard deviation
#' of the prior distribution on the population size.
#' @param mode.prior.sample.proportion scalar; A hyperparameter being the mode
#' of the prior distribution on the sample proportion \eqn{n/N}.
#' @param alpha scalar; A hyperparameter being the first parameter of the beta
#' prior model for the sample proportion. By default this is NULL, meaning that
#' 1 is chosen. it can be any value at least 1 to allow for different levels of
#' uncertainty.
#' @param degreedistribution count; the parametric distribution to use for the
#' individual network sizes (i.e., degrees). The options are \code{cmp},
#' \code{nbinom}, and \code{pln}.  These correspond to the
#' Conway-Maxwell-Poisson, Negative-Binomial, and Poisson-log-normal. The
#' default is \code{cmp}.
#' @param mean.prior.degree scalar; A hyper parameter being the mean degree for
#' the prior distribution for a randomly chosen person. The prior has this
#' mean.
#' @param sd.prior.degree scalar; A hyper parameter being the standard
#' deviation of the degree for a randomly chosen person.  The prior has this
#' standard deviation.
#' @param max.sd.prior.degree scalar; The maximum allowed value of \code{sd.prior.degree}.
#' If the passed or computed value is higher, it is reduced to this value.
#' This is done for numerical stability reasons.
#' @param df.mean.prior scalar; A hyper parameter being the degrees-of-freedom
#' of the prior for the mean. This gives the equivalent sample size that would
#' contain the same amount of information inherent in the prior.
#' @param df.sd.prior scalar; A hyper parameter being the degrees-of-freedom of
#' the prior for the standard deviation. This gives the equivalent sample size
#' that would contain the same amount of information inherent in the prior for
#' the standard deviation.
#' @param Np integer; The overall degree distribution is a mixture of the
#' \code{Np} rates for \code{1:Np} and a parametric degree distribution model
#' truncated below \code{Np}. Thus the model fits the proportions of the
#' population with degree \code{1:Np} each with a separate parameter. This
#' should adjust for an lack-of-fit of the parametric degree distribution model
#' at lower degrees, although it also changes the model away from the
#' parametric degree distribution model.
#' @param nk vector; the vector of counts for the number of people in the
#' sample with degree k. This is usually computed from \eqn{s} automatically as
#' \code{tabulate(s,nbins=K)} and not usually specified by the user.
#' @param n vector; the vector of counts for the number of people in the sample
#' with degree k. This is usually computed from \eqn{s} automatically and not
#' usually specified by the user.
#' @param muproposal scalar; The standard deviation of the proposal
#' distribution for the mean degree.
#' @param sigmaproposal scalar; The standard deviation of the proposal
#' distribution for the standard deviation of the degree.
#' @param burnintheta count; the number of proposals in the Metropolis-Hastings
#' sub-step for the degree distribution parameters (\eqn{\theta}) before any
#' MCMC sampling is done. It typically is set to a modestly large number.
#' @param parallel count; the number of parallel processes to run for the
#' Monte-Carlo sample.  This uses PVM or MPI. The default is 1, that is not to
#' use parallel processing.
#' @param parallel.type The type of parallel processing to use. The options are
#' "PVM" or "MPI". This requires the corresponding type to be installed.
#' @param seed integer; random number integer seed.  Defaults to \code{NULL} to
#' use whatever the state of the random number generator is at the time of the
#' call.
#' @param maxbeta scalar; The maximum allowed value of the \code{beta} parameter.
#' If the implied or computed value is higher, it is reduced to this value.
#' This is done for numerical stability reasons.
#' @param dispersion scalar; dispersion to use in the reported network size
#' compared to the actual network size.
#' @param supplied list; If supplied, is a list with components \code{maxN} and
#' \code{sample}. In this case \code{supplied} is a matrix with a column named
#' \code{N} being a sample from a prior distribution for the population size.
#' The value \code{maxN} specifies the maximum value of the population size, a
#' priori.
#' @param verbose logical; if this is \code{TRUE}, the program will print out
#' additional information, including goodness of fit statistics.
#' @return \code{\link{posteriorsize}} returns a list consisting of the
#' following elements: \item{pop}{vector; The final posterior draw for the
#' degrees of the population. The first \eqn{n} are the sample in sequence and
#' the reminder are non-sequenced.} \item{K}{count; the maximum degree for an
#' individual. This is usually calculated as twice the maximum observed
#' degree.} \item{n}{count; the sample size.} \item{samplesize}{count; the
#' number of Monte-Carlo samples to draw to compute the posterior. This is the
#' number returned by the Metropolis-Hastings algorithm.The default is 1000.}
#' \item{burnin}{count; the number of proposals before any MCMC sampling is
#' done. It typically is set to a fairly large number.} \item{interval}{count;
#' the number of proposals between sampled statistics.} \item{mu}{scalar; The
#' hyper parameter \code{mean.prior.degree} being the mean degree for the prior
#' distribution for a randomly chosen person. The prior has this mean.}
#' \item{sigma}{scalar; The hyper parameter \code{sd.prior.degree} being the
#' standard deviation of the degree for a randomly chosen person. The prior has
#' this standard deviation.} \item{df.mean.prior}{scalar; A hyper parameter
#' being the degrees-of-freedom of the prior for the mean. This gives the
#' equivalent sample size that would contain the same amount of information
#' inherent in the prior.} \item{df.sd.prior}{scalar; A hyper parameter being
#' the degrees-of-freedom of the prior for the standard deviation. This gives
#' the equivalent sample size that would contain the same amount of information
#' inherent in the prior for the standard deviation.} \item{Np}{integer; The
#' overall degree distribution is a mixture of the \code{1:Np} rates and a
#' parametric degree distribution model truncated below Np. Thus the model fits
#' the proportions of the population with degree \code{1:Np} each with a
#' separate parameter. This should adjust for an lack-of-fit of the parametric
#' degree distribution model at lower degrees, although it also changes the
#' model away from the parametric degree distribution model.}
#' \item{muproposal}{scalar; The standard deviation of the proposal
#' distribution for the mean degree.} \item{sigmaproposal}{scalar; The standard
#' deviation of the proposal distribution for the standard deviation of the
#' degree.} \item{N}{vector of length 5; summary statistics for the posterior
#' population size.  \describe{ \item{MAP}{maximum aposteriori value of
#' N} \item{Mean AP}{mean aposteriori value of N} \item{Median
#' AP}{median aposteriori value of N} \item{P025}{the 2.5th
#' percentile of the (posterior) distribution for the N. That is, the lower
#' point on a 95\% probability interval.} \item{P975}{the 97.5th
#' percentile of the (posterior) distribution for the N. That is, the upper
#' point on a 95\% probability interval.} } } \item{maxN}{integer; maximum
#' possible population size. By default this is determined from an upper
#' quantile of the prior distribution.} \item{sample}{matrix of dimension
#' \code{samplesize}\eqn{\times} \code{10} matrix of summary statistics from
#' the posterior. this is also an object of class \code{mcmc} so it can be
#' plotted and summarized via the \code{mcmc.diagnostics} function in the
#' \code{ergm} package (and also the \code{coda} package). The statistics are:
#' \describe{ \item{N}{population size.} \item{mu}{scalar; The mean
#' degree for the prior distribution for a randomly chosen person. The prior
#' has this mean.} \item{sigma}{scalar; The standard deviation of the degree
#' for a randomly chosen person. The prior has this standard deviation.}
#' \item{degree1}{scalar; the number of nodes of degree 1 in the population (it
#' is assumed all nodes have degree 1 or more).} \item{lambda}{scalar; This is
#' only present for the \code{cmp} model. It is the \eqn{\lambda} parameter in
#' the standard parametrization of the Conway-Maxwell-Poisson model for the
#' degree distribution.} \item{nu}{scalar; This is only present for the
#' \code{cmp} model. It is the \eqn{\nu} parameter in the standard
#' parametrization of the Conway-Maxwell-Poisson model for the degree
#' distribution.} } } \item{lpriorm}{vector; the vector of (log) prior
#' probabilities on each value of \eqn{m=N-n} - that is, the number of
#' unobserved members of the population. The values are
#' \code{n:(length(lpriorm)-1+n)}.} \item{burnintheta}{count; the number of
#' proposals in the Metropolis-Hastings sub-step for the degree distribution
#' parameters (\eqn{\theta}) before any MCMC sampling is done. It typically is
#' set to a modestly large number.} \item{verbose}{logical; if this is
#' \code{TRUE}, the program printed out additional information, including
#' goodness of fit statistics.} \item{predictive.degree.count}{vector; a vector
#' of length the maximum degree (\code{K}) (by default \cr \code{K=2*max(sample
#' degree)}).  The \code{k}th entry is the posterior predictive number persons
#' with degree \code{k}.  That is, it is the posterior predictive distribution
#' of the number of people with each degree in the population.}
#' \item{predictive.degree}{vector; a vector of length the maximum degree
#' (\code{K}) (by default \cr \code{K=2*max(sample degree)}).  The \code{k}th entry
#' is the posterior predictive proportion of persons with degree \code{k}.
#' That is, it is the posterior predictive distribution of the proportion of
#' people with each degree in the population.} \item{MAP}{vector of length 6
#' of MAP estimates corresponding to the output \code{sample}. These are:
#' \describe{ \item{N}{population size.} \item{mu}{scalar; The mean
#' degree for the prior distribution for a randomly chosen person. The prior
#' has this mean.} \item{sigma}{scalar; The standard deviation of the degree
#' for a randomly chosen person. The prior has this standard deviation.}
#' \item{degree1}{scalar; the number of nodes of degree 1 in the population (it
#' is assumed all nodes have degree 1 or more).} \item{lambda}{scalar; This is
#' only present for the \code{cmp} model. It is the \eqn{\lambda} parameter in
#' the standard parametrization of the Conway-Maxwell-Poisson model for the
#' degree distribution.} \item{nu}{scalar; This is only present for the
#' \code{cmp} model. It is the \eqn{\nu} parameter in the standard
#' parametrization of the Conway-Maxwell-Poisson model for the degree
#' distribution.} } } \item{mode.prior.sample.proportion}{scalar; A
#' hyperparameter being the mode of the prior distribution on the sample
#' proportion \eqn{n/N}.} \item{median.prior.size}{scalar; A hyperparameter
#' being the mode of the prior distribution on the population size.}
#' \item{mode.prior.size}{scalar; A hyperparameter being the mode of the prior
#' distribution on the population size.} \item{mean.prior.size}{scalar; A
#' hyperparameter being the mean of the prior distribution on the population
#' size.} \item{quartiles.prior.size}{vector of length 2; A pair of
#' hyperparameters being the lower and upper quartiles of the prior
#' distribution on the population size.} \item{degreedistribution}{count; the
#' parametric distribution to use for the individual network sizes (i.e.,
#' degrees). The options are \code{cmp}, \code{nbinom}, and \code{pln}.  These
#' correspond to the Conway-Maxwell-Poisson, Negative-Binomial, and
#' Poisson-log-normal. The default is \code{cmp}.}
#' \item{priorsizedistribution}{character; the type of parametric distribution
#' to use for the prior on population size. The options are \code{beta} (for a
#' Beta prior on the sample proportion (i.e. \eqn{n/N}), \code{nbinom}
#' (Negative-Binomial), \code{pln} (Poisson-log-normal), \code{flat} (uniform),
#' and \code{continuous} (the continuous version of the Beta prior on the
#' sample proportion. The default is \code{beta}. }
#' @section Details on priors: The best way to specify the prior is via the
#' hyperparameter \code{mode.prior.size} which specifies the mode of the prior
#' distribution on the population size. You can alternatively specify the
#' hyperparameter \code{median.prior.size} which specifies the median of the
#' prior distribution on the population size, or \code{mean.prior.sample
#' proportion} which specifies the mean of the prior distribution on the
#' proportion of the population size in the sample or \code{mode.prior.sample
#' proportion} which specifies the mode of the prior distribution on the
#' proportion of the population size in the sample. Finally, you can specify
#' \code{quartiles.prior.size} as a vector of length 2 being the pair of lower
#' and upper quartiles of the prior distribution on the population size.
#' @seealso network, statnet, degreenet
#' @references
#' 
#' Gile, Krista J. (2008) \emph{Inference from Partially-Observed Network
#' Data}, Ph.D. Thesis, Department of Statistics, University of Washington.
#' 
#' Gile, Krista J. and Handcock, Mark S. (2010) \emph{Respondent-Driven
#' Sampling: An Assessment of Current Methodology}, Sociological Methodology
#' 40, 285-327.
#' 
#' Gile, Krista J. and Handcock, Mark S. (2014) \pkg{sspse}: Estimating Hidden 
#' Population Size using Respondent Driven Sampling Data
#' R package, Los Angeles, CA.  Version 0.5, \url{http://hpmrg.org}.
#' 
#' Handcock MS (2003).  \pkg{degreenet}: Models for Skewed Count Distributions
#' Relevant to Networks.  Statnet Project, Seattle, WA.  Version 1.2,
#' \url{http://statnetproject.org}.
#' 
#' Handcock, Mark S., Gile, Krista J. and Mar, Corinne M. (2014)
#' \emph{Estimating Hidden Population Size using Respondent-Driven Sampling
#' Data}, Electronic Journal of Statistics, 8, 1, 1491-1521
#' 
#' Handcock, Mark S., Gile, Krista J. and Mar, Corinne M. (2015)
#' \emph{Estimating the Size of Populations at High Risk for HIV using Respondent-Driven 
#' Sampling Data}, Biometrics.
#' @keywords models
#' @examples
#' 
#' \dontrun{
#' N0 <- 200
#' n <- 100
#' K <- 10
#' 
#' # Create probabilities for a Waring distribution 
#' # with scaling parameter 3 and mean 5, but truncated at K=10.
#' probs <- c(0.33333333,0.19047619,0.11904762,0.07936508,0.05555556,
#'            0.04040404,0.03030303,0.02331002,0.01831502,0.01465201)
#' probs <- probs / sum(probs)
#' 
#' # Look at the degree distribution for the prior
#' # Plot these if you want
#' # plot(x=1:K,y=probs,type="l")
#' # points(x=1:K,y=probs)
#' #
#' # Create a sample
#' #
#' set.seed(1)
#' pop<-sample(1:K, size=N0, replace = TRUE, prob = probs)
#' s<-sample(pop, size=n, replace = FALSE, prob = pop)
#'  
#' out <- posteriorsize(s=s,interval=10)
#' plot(out, HPD.level=0.9,data=pop[s])
#' summary(out, HPD.level=0.9)
#' # Let's look at some MCMC diagnostics
#' plot(out, HPD.level=0.9,mcmc=TRUE)
#' }
#' @export posteriorsize
posteriorsize<-function(s,
		  median.prior.size=NULL,
                  interval=10,
                  burnin=5000,
                  maxN=NULL,
                  K=max(s,na.rm=TRUE),
                  samplesize=1000,
		  quartiles.prior.size=NULL,
		  mean.prior.size=NULL,
		  mode.prior.size=NULL,
		  priorsizedistribution=c("beta","flat","nbinom","pln","supplied"),
		  effective.prior.df=1,
                  sd.prior.size=NULL,
		  mode.prior.sample.proportion=NULL,
		  alpha=NULL,
		  degreedistribution=c("cmp","nbinom","pln"),
                  mean.prior.degree=NULL, sd.prior.degree=NULL, max.sd.prior.degree=4,
                  df.mean.prior=1,df.sd.prior=3,
                  Np=0,
                  nk=NULL,
                  n=length(s),
                  muproposal=0.1, 
                  sigmaproposal=0.15, 
                  burnintheta=500,
                  parallel=1, parallel.type="MPI", seed=NULL, 
                  maxbeta=120, dispersion=0,
                  supplied=list(maxN=maxN),
                  verbose=TRUE){
#
  degreedistribution=match.arg(degreedistribution)
  posfn <- switch(degreedistribution,
                  nbinom=posnbinom,
                  pln=pospln,
		  cmp=poscmp,
		  poscmp)
  remvalues <- !is.na(s)
  if(sum(remvalues) < length(s)){
   warning(paste(length(s)-sum(remvalues),"of",length(s),
  	"sizes values were missing and were removed."), call. = FALSE)
   s <- s[remvalues]
   n <- length(s)
  }
  priorsizedistribution=match.arg(priorsizedistribution)
  if(priorsizedistribution=="nbinom" && missing(mean.prior.size)){
    stop("You need to specify 'mean.prior.size', and possibly 'sd.prior.size' if you use the 'nbinom' prior.") 
  }
  if(is.null(K)){
    K=round(quantile(s,0.90))+1
    degs <- s
    degs[degs>K] <- K
    degs[degs==0]<-1
    ds<-degs
    isnas <- is.na(degs)
    degs <- sum(!isnas)*(degs)/sum(degs,na.rm=TRUE)
    weights <- (1/degs)
    weights[is.na(weights)] <- 0
    mean.pd <- sum(ds*weights)/sum(weights)
    sd.pd <- sum(ds*ds*weights)/sum(weights)
    sd.pd <- sqrt(sd.pd - mean.pd^2)
    if(sd.pd > max.sd.prior.degree*mean.pd){
     sd.pd <- min(max.sd.prior.degree*mean.pd, sd.pd)
    }
    xv <- ds
#   xp <- weights*ds
    xp <- weights
    xp <- length(xp)*xp/sum(xp)
    txp <- tapply(xp,xv,sum)
    txv <- tapply(xv,xv,median)
    fit <- cmpmle(txv,txp,cutoff=1,cutabove=K-1,guess=c(mean.pd, sd.pd))
    y=dcmp.natural(v=fit,x=(0:max(s)))
    K=(0:max(s))[which.max(cumsum(y)>0.99)]
#   K=round(quantile(s,0.99))
  }
  cat(sprintf("The cap on influence of the personal network size is K = %d.\n",K))
  if(is.null(mean.prior.degree)){
    degs <- s
    degs[degs>K] <- K
    degs[degs==0]<-1
    ds<-degs
    isnas <- is.na(degs)
    degs <- sum(!isnas)*(degs)/sum(degs,na.rm=TRUE)
    weights <- (1/degs)
    weights[is.na(weights)] <- 0
    mean.prior.degree <- sum(ds*weights)/sum(weights)
    if(is.null(sd.prior.degree)){
     sd.prior.degree <- sum(ds*ds*weights)/sum(weights)
     sd.prior.degree <- sqrt(sd.prior.degree - mean.prior.degree^2)
    }
    xv <- ds
#   xp <- weights*ds
    xp <- weights
    xp[is.na(xp)] <- 0
    xp <- length(xp)*xp/sum(xp)
    txp <- tapply(xp,xv,sum)
    txv <- tapply(xv,xv,median)
    fit <- cmpmle(txv,txp,cutoff=1,cutabove=K-1,
            guess=c(mean.prior.degree,sd.prior.degree))
    fit <- cmp.mu(fit,max.mu=5*mean.prior.degree)
    if(verbose){
      cat(sprintf("The preliminary empirical value of the mean of the prior distribution for degree is %f.\n",mean.prior.degree))
      cat(sprintf("The preliminary empirical value of the s.d. of the prior distribution for degree is %f.\n",sd.prior.degree))
    }
    mean.prior.degree = fit[1]
    sd.prior.degree = fit[2]
  }
  if(verbose){
    cat(sprintf("The mean of the prior distribution for degree is %f.\n",mean.prior.degree))
    cat(sprintf("The s.d. of the prior distribution for degree is %f.\n",sd.prior.degree))
  }
  if(sd.prior.degree > max.sd.prior.degree*mean.prior.degree){
    sd.prior.degree <- min(max.sd.prior.degree*mean.prior.degree, sd.prior.degree)
    cat(sprintf("The suggested s.d. of the prior distribution for degree is too large and has been reduced to the more reasonable %f.\n",sd.prior.degree))
  }
  ### are we running the job in parallel (parallel > 1), if not just 
  #   call the degree specific function
  if(parallel==1){
      Cret <- posfn(s=s,K=K,nk=nk,n=n,maxN=maxN,
                    mean.prior.degree=mean.prior.degree,df.mean.prior=df.mean.prior,
                    sd.prior.degree=sd.prior.degree,df.sd.prior=df.sd.prior,
                    muproposal=muproposal, sigmaproposal=sigmaproposal, 
		    Np=Np,
                    samplesize=samplesize,burnin=burnin,interval=interval,
		    burnintheta=burnintheta,
		    priorsizedistribution=priorsizedistribution,
		    mean.prior.size=mean.prior.size, sd.prior.size=sd.prior.size,
		    mode.prior.sample.proportion=mode.prior.sample.proportion,
		    median.prior.size=median.prior.size,
		    mode.prior.size=mode.prior.size,
		    quartiles.prior.size=quartiles.prior.size,
                    effective.prior.df=effective.prior.df,
                    alpha=alpha,
                    seed=seed,
                    supplied=supplied,
                    maxbeta=maxbeta,
                    dispersion=dispersion)
  }
  ### since running job in parallel, start pvm (if not already running)
  else{
    cl <- sspse::beginparallel(parallel,type=parallel.type)
    ### divide the samplesize by the number of parallel runs (number of MCMC samples)
    samplesize.parallel=round(samplesize/parallel)
    ### cluster call, send following to each of the virtual machines, posnbinom function
    ### with it's arguments
    outlist <- parallel::clusterCall(cl, posfn,
      s=s,K=K,nk=nk,n=n,maxN=maxN,
      mean.prior.degree=mean.prior.degree,df.mean.prior=df.mean.prior,
      sd.prior.degree=sd.prior.degree,df.sd.prior=df.sd.prior,
      muproposal=muproposal, sigmaproposal=sigmaproposal, 
      Np=Np,
      samplesize=samplesize.parallel,burnin=burnin,interval=interval,
      burnintheta=burnintheta,
      priorsizedistribution=priorsizedistribution,
      mean.prior.size=mean.prior.size, sd.prior.size=sd.prior.size,
      mode.prior.sample.proportion=mode.prior.sample.proportion,
      median.prior.size=median.prior.size,
      mode.prior.size=mode.prior.size,
      quartiles.prior.size=quartiles.prior.size,
      effective.prior.df=effective.prior.df,
      alpha=alpha,
      supplied=supplied,
      dispersion=dispersion,
      maxbeta=maxbeta)
#
#   Process the results
#
    ### Snow returns a list of length parallel where each element is the return of each posfn
    ### Following loops combines the separate MCMC samples into 1 using rbind, creating a matrix
    Cret <- outlist[[1]]
    Nparallel <- length(outlist)
    Cret$samplesize <- samplesize
    for(i in (2 : Nparallel)){
     z <- outlist[[i]]
     Cret$sample <- rbind(Cret$sample,z$sample)
     Cret$predictive.degree.count<-Cret$predictive.degree.count+z$predictive.degree.count
     Cret$predictive.degree<-Cret$predictive.degree+z$predictive.degree
    }
    Cret$predictive.degree.count<-Cret$predictive.degree.count/Nparallel
    Cret$predictive.degree<-Cret$predictive.degree/Nparallel
    #
    degnames <- NULL
    if(Np>0){degnames <- c(degnames,paste("pdeg",1:Np,sep=""))}
    colnamessample <- c("N","mu","sigma","degree1","totalsize")
    if(length(degnames)>0){
     colnamessample <- c(colnamessample, degnames)
    }
    if(degreedistribution=="cmp"){
     colnamessample <- c(colnamessample, c("lambda","nu"))
    }
    colnames(Cret$sample) <- colnamessample
    m <- apply(Cret$sample,2,median,na.rm=TRUE)
    Cret$sample[is.na(Cret$sample[,"mu"]),"mu"] <- m["mu"]
    Cret$sample[is.na(Cret$sample[,"sigma"]),"sigma"] <- m["sigma"]
#   Any NA and NaN are typically in pdeg and so should be 0.
    Cret$sample[is.na(Cret$sample)] <- 0
    Cret$sample[is.nan(Cret$sample)] <- 0
    
    ### Coda package which does MCMC diagnostics, requires certain attributes of MCMC sample
    endrun <- burnin+interval*(samplesize)
    attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    
#   ### Remove the padding from the last draws from the populations of degrees
#   Nlastpos <- Cret$sample[nrow(Cret$sample),"N"]
#   Cret$pop<-Cret$pop[1:Nlastpop]

    ### compute modes of posterior samples,Maximum A Posterior (MAP) values N, mu, sigma, degree1
    Cret$MAP <- apply(Cret$sample,2,mode.density)
    Cret$MAP["N"] <- mode.density(Cret$sample[,"N"],lbound=n,ubound=Cret$maxN)
    if(verbose){
     cat("parallel samplesize=", parallel,"by", samplesize.parallel,"\n")
    }
    
    ### stop cluster
    sspse::endparallel(cl,type=parallel.type)
  }
  Cret$N <- c(Cret$MAP["N"], 
              mean(Cret$sample[,"N"]),
              median(Cret$sample[,"N"]),
	      quantile(Cret$sample[,"N"],c(0.025,0.975)))
  names(Cret$N) <- c("MAP","Mean AP","Median AP","P025","P975")
  #
  if(Cret$predictive.degree[length(Cret$predictive.degree)] > 0.01){
   warning("There is a non-trivial proportion of the posterior mass on very high degrees. This may indicate convergence problems in the MCMC.")
  }
  Cret$degreedistribution <- degreedistribution
  Cret$priorsizedistribution <- priorsizedistribution
# Cret$mean.prior.size <- mean.prior.size
  ### return result
  class(Cret) <- "sspse"
  Cret
}

posize.warning <- "POSTERIOR SIZE CALCULATION FAILED" #added for posteriorsize dialog to hide error message unless needed
