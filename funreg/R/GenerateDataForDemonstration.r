#' @title Generate data for some demonstration examples
#' @description Simulates a dataset with two functional covariates, four 
#' subject-level scalar covariates, and a binary outcome.
#' @param nsub The number of subjects in the simulated dataset.
#' @param b0.true The true value of the intercept.
#' @param b1.true The true value of the first covariate.
#' @param b2.true The true value of the second covariate.
#' @param b3.true The true value of the third covariate.
#' @param b4.true The true value of the fourth covariate.
#' @param nobs The total number of possible observation times.
#' @param observe.rate The average proportion of those possible
#' times at which any given subject is observed.
#' @return Returns a \code{data.frame} representing \code{nobs}
#'  measurements for each subject. The rows of this \code{data.frame}
#'  tell the values of two time-varying covariates on a dense grid
#'  of \code{nobs} observation times. It also contains an
#'  \code{id} variable, four subject-level covariates
#'   (\code{s1}, ..., \code{s4}) and one subject-level
#' response (\code{y}), which are replicated for each observation.
#'  For each observation, there is also its observation 
#'  time \code{time}, there are both the smooth latent value of the covariates 
#' (\code{true.x1} and \code{true.x2}) and 
#' versions observed with error (\code{x1} and \code{x2}), and there are 
#' also the local values of the functional regression coefficients 
#' (\code{true.betafn1} and \code{true.betafn2}).  Lastly, 
#' each row has a random value for \code{include.in.subsample},
#' telling whether it should be considered as an observed data
#' point (versus an unobserved moment in the simulated subject's life).
#' \code{include.in.subsample} is simply generated as a Bernoulli random variable with 
#' success probability \code{observe.rate}.  
#' @note \code{nobs} is the number of simulated data rows per
#' simulated subject.  It  
#' should be selected to be large because \code{x} covariates are conceptually
#'  supposed to be smooth functions of time. However, in the
#' simulated data analyses we actually only use a small random
#' subset of the generated time points, because this is more
#' realistic for many behavioral and medical science datasets.
#' Thus, the number of possible observation times per subject
#' is \code{nobs}, and the mean number of actual observation
#' times per subject is \code{nobs} times \code{observe.rate}. 
#' This smaller 'observed' dataset can be obtained by 
#' deleting from the dataset those observations having 
#' \code{include.in.subsample==FALSE}.
#' @importFrom mvtnorm rmvnorm
#'@export
generate.data.for.demonstration <- function(nsub=400,
                      b0.true=-5,
                    
                                            b1.true=0,
                      b2.true=+1,
                      b3.true=-1,
                      b4.true=+1,
                      nobs = 500,
                      observe.rate=.1) {
  nbins <- 35;
  min.time.grid <- 0;
  max.time.grid <- 10;
  cov.scale <- 2;
  x.scale <- 1;
  x.error.sd <- 1;
  time.grid <- seq(min.time.grid,max.time.grid,length=nobs);
  time.bins <- seq(min.time.grid,max.time.grid,length=nbins);
  betafn1.true <- cos(time.grid*pi/2/(max.time.grid-min.time.grid));
  betafn1.true.by.bins <- cos(time.bins*pi/2/(max.time.grid-min.time.grid));
  betafn2.true <- rep(0,length(betafn1.true));
  betafn2.true.by.bins <- rep(0,length(betafn1.true.by.bins));

  nobs <- length(time.grid);
  times <- matrix(rep(time.grid,times=nsub),nrow=nsub,byrow=TRUE);
  expected.crossings <- (1/(2*pi*cov.scale)*(max(time.grid)-min(time.grid)));
          # see Rasmussen & Williams, p. 83;
  gaussian.cov.mat <- matrix(NA,nobs,nobs);
  for (i in 1:nobs) {
    for (j in 1:nobs) {
      gaussian.cov.mat[i,j] <- x.scale*exp(-(1/(2*cov.scale^2))*
                       (time.grid[i]-time.grid[j])^2);
    }
  }
  gaussian.cov.mat.by.bins <- matrix(NA,nbins,nbins);
  for (i in 1:length(time.bins)) {
    for (j in 1:length(time.bins)) {
      gaussian.cov.mat.by.bins[i,j] <- x.scale*exp(-(1/(2*cov.scale^2))*
                       (time.bins[i]-time.bins[j])^2);
    }
  }
  for (i in 1:nobs) {
    for (j in 1:nobs) {
      gaussian.cov.mat[i,j] <- x.scale*exp(-(1/(2*cov.scale^2))*
                       (time.grid[i]-time.grid[j])^2);
    }
  }
  start.time <- Sys.time();
  true.x1 <- rmvnorm(n=nsub,sigma=gaussian.cov.mat,method="svd")+1;
  x1 <- true.x1 + matrix(rnorm(nsub*nobs,0,x.error.sd),nrow=nsub,byrow=TRUE);
  true.x2 <- rmvnorm(n=nsub,sigma=gaussian.cov.mat,method="svd")+2;
  x2 <- true.x2 + matrix(rnorm(nsub*nobs,0,x.error.sd),nrow=nsub,byrow=TRUE);
  #######################
  # Generate the y values
  step.size <- time.grid[2]-time.grid[1];
  product.x1.betafn1 <- rep(NA,nsub);
  product.x2.betafn2 <- rep(NA,nsub);
  for (i in 1:nsub) {
    product.x1.betafn1[i] <- sum(true.x1[i,] * betafn1.true * step.size);
    product.x2.betafn2[i] <- sum(true.x2[i,] * betafn2.true * step.size);
  }
  s1 <- rbinom(nsub,1,.5);
  s2 <- rbinom(nsub,1,.5);
  s3 <- rnorm(nsub);
  s4 <- rnorm(nsub);
  linear.predictor <- b0.true +
            s1*b1.true +
            s2*b2.true +
            s3*b3.true +
            s4*b4.true +
            product.x1.betafn1+
            product.x2.betafn2;
  inv.logit.linear.predictor <- exp(linear.predictor)/(1+exp(linear.predictor));
  y <- rbinom(nsub,1,inv.logit.linear.predictor);
  subject.level.data <- cbind(id=1:nsub,
                s1=s1,
                s2=s2,
                s3=s3,
                s4=s4,
                y=y);
  assessment.level.data <- cbind(id=as.vector(t(row(x1))),
                  time=as.vector(t(times)),
                  true.x1=as.vector(t(true.x1)),
                  true.x2=as.vector(t(true.x2)),
                  true.betafn1=rep(betafn1.true,times=nsub),
                  true.betafn2=rep(betafn2.true,times=nsub),
                  x1=as.vector(t(x1)),
                  x2=as.vector(t(x2)),
                  include.in.subsample=rbinom(nsub*nobs,1,observe.rate));
  full.data <- merge(subject.level.data,
            assessment.level.data,
            by="id");
  the.data <- full.data[which(full.data$include.in.subsample==TRUE),1:13];
  return(the.data);
}