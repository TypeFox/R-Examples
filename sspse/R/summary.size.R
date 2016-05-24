#' Summarizing Population Size Estimation Model Fits
#' 
#' This is the \code{summary} method for class \code{"sspse"} objects.
#' These objects encapsulate an estimate of the posterior distribution of
#' the population size based on data collected by Respondent Driven Sampling.
#' The approach approximates the RDS via the Sequential Sampling model of Gile
#' (2008).  As such, it is referred to as the Sequential Sampling - Population Size Estimate (SS-PSE).
#' It uses the order of selection of the sample to provide information
#' on the distribution of network sizes over the population members.
#' \code{summary} method for class \code{"sspse"}. posterior distribution of
#' the population size based on data collected by Respondent Driven Sampling.
#' The approach approximates the RDS via the Sequential Sampling model of Gile
#' (2008). As such, it is referred to as the Sequential Sampling - Population Size Estimate (SS-PSE).
#' It uses the order of selection of the sample to provide information
#' on the distribution of network sizes over the population members.
#' 
#' \code{print.summary.sspse} tries to be smart about formatting the
#' coefficients, standard errors, etc. and additionally gives
#' \sQuote{significance stars} if \code{signif.stars} is \code{TRUE}.
#' 
#' Aliased coefficients are omitted in the returned object but restored by the
#' \code{print} method.
#' 
#' Correlations are printed to two decimal places (or symbolically): to see the
#' actual correlations print \code{summary(object)$correlation} directly.
#' 
#' @aliases summary.sspse
#' @param object an object of class \code{"sspse"}, usually, a result of a call
#' to \code{\link{posteriorsize}}.
#' @param support the number of equally-spaced points to use for the support of
#' the estimated posterior density function.
#' @param HPD.level numeric; probability level of the highest probability
#' density interval determined from the estimated posterior.
#' @param \dots further arguments passed to or from other methods.
#' @return The function \code{summary.sspse} computes and returns a two row matrix of
#' summary statistics of the prior and estimated posterior distributions. The rows correspond to the \code{Prior} and the
#' \code{Posterior}, respectively. 
#' The rows names are \code{Mean}, \code{Median}, \code{Mode}, \code{25\%}, \code{75\%}, and \code{90\%}.
#' These correspond to the distributional mean, median, mode, lower quartile, upper quartile and 90\% quantile, respectively. 
#' @seealso The model fitting function \code{\link{posteriorsize}},
#' \code{\link{summary}}.
#' 
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
#' 
#' @method summary sspse
#' @export
summary.sspse <- function(object, support=1000, HPD.level=0.95,...){
#summary.sspse <- function(object, ...){
  p.args <- as.list( sys.call() )[-c(1,2)]
  formal.args<-formals(sys.function())[-1]
# control <- list(support=1000,HPD.level=0.95)

  control<-list()
  names.formal.args <- names(formal.args)
  names.formal.args <- names.formal.args[-match("...",names.formal.args)]
  for(arg in names.formal.args){ control[arg]<-list(get(arg)) }
  for(arg in names(p.args)){ control[arg]<-list(get(arg)) }

#suppressMessages(require(locfit, quietly=TRUE))
out <- object$sample
if(is.null(out)){
  object$n <- min(object$x)
  object$lpriorm <- log(object$lprior)
}
if(!is.null(out)){
  outN <- out[,"N"]
  #a=locfit( ~ lp(outN, nn=0.35, h=0, maxk=500))
  xp <- seq(object$n,object$maxN, length=control$support)
# a=locfit::locfit( ~ lp(outN,nn=0.5))
# posdensN <- predict(a, newdata=xp)
  a=bgk_kde(outN,n=2^(ceiling(log((object$maxN-object$n))/log(2))),MIN=object$n,MAX=object$maxN)
  posdensN <- spline(x=a[1,],y=a[2,],xout=xp)$y
  posdensN <- control$support*posdensN / ((object$maxN-object$n)*sum(posdensN))
  # Next from coda
  # hpd <- HPDinterval(object$sample[,"N"])[1:2]
  # MSH using locfit
  cy <- cumsum(posdensN/sum(posdensN))
  hpd <- c(xp[which.max(cy>((1-control$HPD.level)/2))],
           xp[which.max(cy>((1+control$HPD.level)/2))])
  if(is.na(hpd[1])) hpd[1] <- xp[1]
  if(is.na(hpd[2])) hpd[2] <- xp[length(xp)]
#
  map <- xp[which.max(posdensN)]
  mp <- sum(xp*posdensN)/sum(posdensN)
  l90 <- xp[which.max(cy>0.9)]
  l50 <- xp[which.max(cy>0.5)]
  l25 <- xp[which.max(cy>0.25)]
  l75 <- xp[which.max(cy>0.75)]
}
#
lpriorm <- exp(object$lpriorm-max(object$lpriorm))
lpriorm <- lpriorm[object$n+(1:length(lpriorm)) > object$n & object$n+(1:length(lpriorm)) < object$maxN]
lpriorm <- lpriorm / sum(lpriorm)
cy <- cumsum(lpriorm)
xp <- seq(object$n,object$maxN)
pl025 <- xp[which.max(cy>((1-control$HPD.level)/2))]
pl95  <- xp[which.max(cy>((1+control$HPD.level)/2))]
pl90  <- xp[which.max(cy>0.9)]
#
#cat(sprintf("Prior:\nMean = %d, Median = %d, Mode = %d, 90%% = %d, 25%% = %d, 75%% = %d.\n",
# round(object$mean.prior.size), round(object$median.prior.size), round(object$mode.prior.size), round(pl90), round(object$quartiles.prior.size[1]), round(object$quartiles.prior.size[2])))
#cat(sprintf("Posterior:\nMean = %d, Median = %d, MAP = %d, 90%% = %d, HPD = (%d, %d).\n",
# round(mp),round(l50),round(map),round(l90),round(hpd[1]),round(hpd[2])))
#
if(!is.null(out)){
 res <- matrix(c(
  round(object$mean.prior.size), round(object$median.prior.size), round(object$mode.prior.size),
  round(object$quartiles.prior.size[1]),round(object$quartiles.prior.size[2]),
  round(pl90),round(pl025),round(pl95),
  round(mp),round(l50),round(map),round(l25),round(l75),round(l90),round(hpd[1]),round(hpd[2])),byrow=TRUE,nrow=2)
  rownames(res) <- c("Prior","Posterior")
  colnames(res) <- c("Mean","Median","Mode","25%","75%","90%",
    paste(round(100*(1-control$HPD.level)/2,1),"%",sep=""),
    paste(round(100*(1+control$HPD.level)/2,1),"%",sep=""))
}else{
 res <- matrix(c(
  round(object$mean.prior.size), round(object$median.prior.size), round(object$mode.prior.size), round(object$quartiles.prior.size[1]), round(object$quartiles.prior.size[2]), round(pl90), round(pl025), round(pl95)
  ),byrow=TRUE,nrow=1)
  rownames(res) <- c("Prior")
  colnames(res) <- c("Mean","Median","Mode","25%","75%","90%",
    paste(round(100*(1-control$HPD.level)/2,1),"%",sep=""),
    paste(round(100*(1+control$HPD.level)/2,1),"%",sep=""))
}
res <- as.data.frame(res)
if(!is.null(out)){
  attr(res, "heading") <- "Summary of Population Size Estimation"
}else{
  attr(res, "heading") <- "Summary of Population Size Prior"
}
class(res) <- c("Anova","data.frame")
res
}
