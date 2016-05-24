#' plot.mcmcComposite plots the result of a MCMC search
#' @title Plot the result of a mcmcComposite object
#' @author Marc Girondot
#' @return None
#' @param x A mcmcComposite object
#' @param chain The chain to use
#' @param parameters Name of parameters or their number (see description)
#' @param legend If FALSE, the legend is not shown
#' @param scale.prior If TRUE, the prior is scaled at the same size as posterior
#' @param ... Graphical parameters to be send to hist()
#' @family mcmcComposite functions
#' @description Plot the results within a mcmcCompositve object.\cr
#' The parameters to use can be called by:\cr
#' \code{parameters="all"}\cr
#' \code{parameters=1:4}\cr
#' \code{parameters=c("PAR1", "PAR2", "PAR5")}\cr
#' \code{parameters=c(TRUE, TRUE, FALSE, TRUE)}\cr
#' @examples 
#' \dontrun{
#' library(HelpersMG)
#' require(coda)
#' x <- rnorm(30, 10, 2)
#' dnormx <- function(x, par) return(-sum(dnorm(x, mean=par['mean'], sd=par['sd'], log=TRUE)))
#' parameters_mcmc <- data.frame(Density=c('dnorm', 'dlnorm'), 
#' Prior1=c(10, 0.5), Prior2=c(2, 0.5), SDProp=c(1, 1), 
#' Min=c(-3, 0), Max=c(100, 10), Init=c(10, 2), stringsAsFactors = FALSE, 
#' row.names=c('mean', 'sd'))
#' mcmc_run <- MHalgoGen(n.iter=1000, parameters=parameters_mcmc, data=x, 
#' likelihood=dnormx, n.chains=1, n.adapt=100, thin=1, trace=1)
#' plot(mcmc_run, xlim=c(0, 20))
#' plot(mcmc_run, xlim=c(0, 10), parameters="sd")
#' mcmcforcoda <- as.mcmc(mcmc_run)
#' #' heidel.diag(mcmcforcoda)
#' raftery.diag(mcmcforcoda)
#' autocorr.diag(mcmcforcoda)
#' acf(mcmcforcoda[[1]][,"mean"], lag.max=20, bty="n", las=1)
#' acf(mcmcforcoda[[1]][,"sd"], lag.max=20, bty="n", las=1)
#' batchSE(mcmcforcoda, batchSize=100)
#' # The batch standard error procedure is usually thought to 
#' # be not as accurate as the time series methods used in summary
#' summary(mcmcforcoda)$statistics[,"Time-series SE"]
#' summary(mcmc_run)
#' as.parameters(mcmc_run)
#' lastp <- as.parameters(mcmc_run, index="last")
#' parameters_mcmc[,"Init"] <- lastp
#' # The n.adapt set to 1 is used to not record the first set of parameters
#' # then it is not duplicated (as it is also the last one for 
#' # the object mcmc_run)
#' mcmc_run2 <- MHalgoGen(n.iter=1000, parameters=parameters_mcmc, data=x, 
#' likelihood=dnormx, n.chains=1, n.adapt=1, thin=1, trace=1)
#' mcmc_run3 <- merge(mcmc_run, mcmc_run2)
#' ####### no adaptation, n.adapt must be 0
#' parameters_mcmc[,"Init"] <- c(mean(x), sd(x))
#' mcmc_run3 <- MHalgoGen(n.iter=1000, parameters=parameters_mcmc, data=x, 
#' likelihood=dnormx, n.chains=1, n.adapt=0, thin=1, trace=1)
#' }
#' @method plot mcmcComposite
#' @export

plot.mcmcComposite <- function(x, ... , chain=1, parameters=1, scale.prior=FALSE, legend=TRUE) {

resultMCMC <- x

mcmc <- resultMCMC[["resultMCMC"]]
Parameters_L <- resultMCMC[["parametersMCMC"]]
Parameters <- Parameters_L[[1]]
n.iter <- Parameters_L[[2]]
n.chains <- Parameters_L[[3]]
n.adapt <- Parameters_L[[4]]
thin <- Parameters_L[[5]]

possible <- rownames(Parameters)
NbTS <- length(possible)

if (parameters[[1]]=="all") {
	series<-rep(TRUE, NbTS)
} else {
	if (any(!is.logical(parameters))) {
		if (is.numeric(parameters)) {
# Même si un nombre de paramètres plus important est indiqué, ne provoque pas d'erreur
			series <- rep(FALSE, max(NbTS, length(parameters)))
			series[parameters] <- TRUE
			series <- series[1:NbTS]
		} else {
			series <- (possible==parameters)
		}
	} else {
# c'est des valeurs logiques, je verifie si le bon nombre, sinon je recycle
		if (length(parameters)!=NbTS) {
			series <- rep(parameters, NbTS)
			series <- series[1:NbTS]
		}
	}
}

seriesx <- series

nbseries <- length(seriesx[seriesx==TRUE])

for(variable in 1:NbTS) {

if (seriesx[variable]) {

nitercorrige <- floor(n.iter/thin)
vals <- mcmc[[chain]][,variable]
# c'est quoi ça ??? 6/10/2012
# vals <- vals[(length(vals)-nitercorrige):length(vals)]

	x <- vals

if (Parameters[variable, "Density"]=="dunif") {
  # 22/8/2014
  # je suis en uniforme. Il faut que je limite les breaks
  # 9/3/2015 Si l'un est négatif ou nul, je dois faire un shift d'origine
  shift <- as.numeric(Parameters[variable,"Prior1"])
  
  
  mx <- as.numeric(Parameters[variable,"Prior2"])-shift
  mn <- as.numeric(Parameters[variable,"Prior1"])-shift
  mxlog <- 10^(floor(log10(mx))-1)
  # C'est quoi ???
  xl <- as.numeric(c(mn, mx))+c(-mxlog, mxlog)
  br1 <- mn
  br2 <- mx
  interval <- (mx-mn)/20
  decalage <- br1 %% interval
  
  br <- seq(from=decalage+shift, to=mx+shift, by=interval)
  
  # tpt <- list(las=1, xlim=c(0,30), breaks=c(0, 1.00095, 2.0009, 3.00085, 4.0008, 5.00075, 6.0007, 7.00065, 8.0006, 9.00055, 10.0005, 11.00045, 12.0004, 13.00035, 14.0003, 15.00025, 16.0002, 17.00015, 18.0001, 19.00005, 20))
  tpt <- list(...)
  
	L <- modifyList(list(ylab="Density", xlab=rownames(Parameters)[[variable]], main="", freq=FALSE, 
	xlim=xl, breaks=br), modifyList(list(x=x), tpt)) 
} else {
	L <- modifyList(list(ylab="Density", xlab=rownames(Parameters)[[variable]], main="", freq=FALSE), modifyList(list(x=x), tpt)) 
}

	do.call(hist, L) 
	
	scl <- ScalePreviousPlot()

  sequence <- seq(from=scl$xlim[1], to=scl$xlim[2], length=200)
  p1 <- as.numeric(Parameters[variable,"Prior1"])
  p2 <- as.numeric(Parameters[variable,"Prior2"])
  y <- get(Parameters[variable, "Density"])(sequence, p1, p2)
  yl <- c(0, max(y))
  
if (Parameters[variable, "Density"]!="dunif") {
  par(new=TRUE)
  if (scale.prior) {
    plot(sequence, y, type="l", col="red", axes=FALSE, bty="n", xlab="", ylab="", main="", ylim=yl)
  } else {
    plot(sequence, y, type="l", col="red", axes=FALSE, bty="n", xlab="", ylab="", main="", ylim=scl$ylim[1:2])
  }
} else {
  if (scale.prior) {
  segments(x0=scl$xlim[1], x1=p1, y0=0, y1=0, col="red")
  segments(x0=p1, x1=p1, y0=0, y1=scl$ylim[2], col="red")
  segments(x0=p1, x1=p2, y0=scl$ylim[2], y1=scl$ylim[2], col="red")
  segments(x0=p2, x1=p2, y0=scl$ylim[2], y1=0, col="red")
  segments(x0=p2, x1=scl$xlim[2], y0=0, y1=0, col="red")
  } else {
    segments(x0=scl$xlim[1], x1=p1, y0=0, y1=0, col="red")
    segments(x0=p1, x1=p1, y0=0, y1=yl[2], col="red")
    segments(x0=p1, x1=p2, y0=yl[2], y1=yl[2], col="red")
    segments(x0=p2, x1=p2, y0=yl[2], y1=0, col="red")
    segments(x0=p2, x1=scl$xlim[2], y0=0, y1=0, col="red")
  }
}

if (legend) {
	legend("topright", c("Prior", "Posterior"), lty=1, col=c('red', 'black'), bty = "n")
}
    }
  }
}
