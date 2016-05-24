#' Fit a Spatial Occupancy Model via Gibbs Sampling
#' 
#' This package contains functions that fit a spatial occupancy model where the
#' true occupancy is a function of a spatial process.  An efficient Gibbs
#' sampling algorithm is used by formulating the detection and occupancy
#' process models with a probit model instead of the traditional logit based
#' model.
#' 
#' \tabular{ll}{ Package: \tab stocc\cr Type: \tab Package\cr Version: \tab
#' 1.30\cr Date: \tab August 21, 2015\cr License: \tab file LICENSE\cr LazyLoad: \tab
#' yes\cr }
#' 
#' @name stocc-package
#' @aliases stocc-package stocc
#' @docType package
#' @author Devin S. Johnson
#' 
#' Maintainer: Devin S. Johnson <devin.johnson@@noaa.gov>
#' 
#' @importFrom graphics boxplot 
#' @importFrom stats density dist knots model.matrix 
#'           pnorm rbinom rgamma rnorm sd terms
NULL


#' A simulated data set of environmental covariates
#' 
#' This data represents a simulated study area. The study area is a 40 x 40
#' grid of pixels. There are two variables, a factor variable (e.g., a habitat
#' layer), as well as, a continuous covariate.
#' 
#' 
#' @name habData
#' @docType data
#' @format A data frame with 1600 observations on the following 5 variables.
#' \describe{ \item{site}{Site labels} 
#' \item{x}{Longitude
#' coordinate} \item{y}{Latitude coordinate} \item{habCov1}{a
#' factor with levels \code{1} \code{2} \code{3}} \item{habCov2}{a numeric
#' vector} }
#' @keywords datasets
#' @examples
#' 
#' data(habData)
#' image(x=seq(0.5,39.5,1), y=seq(0.5,39.5,1), 
#' 	z=t(matrix(as.numeric(habData$habCov1),40)), main="habData: Factor environmental covariate",
#' 	xlab="x", ylab="y", col=rainbow(3))
#' 
#' 
#' dev.new()
#' image(x=seq(0.5,39.5,1), y=seq(0.5,39.5,1), 
#' 	z=t(matrix(habData$habCov2,40)), main="habData: Continuous environmental covariate",
#' 	xlab="x", ylab="y", col=terrain.colors(50))
#' 
#' 
NULL





#' Simulated occupancy for the 40 x 40 study area.
#' 
#' This data represnts truth with regards to occupancy in the simulated study
#' area. The probability of occupancy was simulated as \code{pnorm(0, X%*%gamma
#' + K alpha, 1, lower=FALSE)}, where \code{K} and \code{alpha} were constructed 
#' from a reduced rank is an ICAR process with precision
#' (\code{tau}) = 0.3 and \code{gamma = c(-1, 0, 0, 1)}
#' 
#' 
#' @name occupancyData
#' @docType data
#' @format A data frame with 1600 observations on the following 5 variables.
#' \describe{ 
#' \item{site}{Site labels} 
#' \item{x}{Longitude coordinate} 
#' \item{y}{Latitude coordinate} 
#' \item{psi}{True probability of occupancy} 
#' \item{psi.fix}{The fixed effects portion of the occupancy process map}
#' \item{occ}{True realized occupancy} 
#' }
#' @examples
#' 
#' data(occupancyData)
#' ##
#' ## Blue points represent realized occupancy.
#' ##
#' image(x=seq(0.5,39.5,1), y=seq(0.5,39.5,1), z=t(matrix(occupancyData$psi,40)), 
#' 	xlab="x", ylab="y", main="Occupancy process with realized occupancy")
#' points(occupancyData$x[occupancyData$occ==1], occupancyData$y[occupancyData$occ==1], 
#'  pch=20, cex=0.25, col="blue")
#' 
NULL



#' Simulated occupancy survey data
#' 
#' Data set representing a simulated survey of the 40 x 40 study area.
#' Approximately 1/3 of the 1600 sites were visited at least once. Those sites
#' that were surveyed were visited a random number of times with an average of
#' 2.5 visits. Detection was simulated as a function of 2 covariates, a
#' continuous one (cov1) and a factor (cov2). These are NOT the same as the
#' cov1 and cov2 of the habData data frame. The coefficients used were
#' \code{beta = c(1, 0, 0.5, 1, 0)}. Thus detection given occupancy of site i
#' at time j = \code{pnorm(0,X\%*\%beta,lower=FALSE)}.
#' 
#' 
#' @name visitData
#' @docType data
#' @format A data frame with 1340 observations on the following 6 variables.
#' \describe{ \item{site}{Site labels} \item{x}{Longitude
#' coordinate} \item{y}{Latitude coordinate} \item{detCov1}{a
#' numeric vector} \item{detCov2}{a factor with levels \code{0} \code{1}
#' \code{2} \code{3}} \item{obs}{a numeric vector} }
#' @keywords datasets
#' @examples
#' 
#' data(visitData)
#' data(occupancyData)
#' ##
#' ## Blue points represent visited sites and green circles represent confirmed occupancy.
#' ##
#' image(x=seq(0.5,39.5,1), y=seq(0.5,39.5,1), z=t(matrix(occupancyData$psi,40)), 
#' 	xlab="x", ylab="y", main="Occupancy process with visits")
#' points(visitData$x[visitData$obs==1], visitData$y[visitData$obs==1], col="green")
#' points(visitData$x, visitData$y, col="blue", pch=20, cex=0.25)
#' 
NULL



.onAttach <- function(library, pkgname)
{
  info <-utils::packageDescription(pkgname)
  package <- info$Package
  version <- info$Version
  date <- info$Date
  packageStartupMessage(
    paste(paste(package, version, paste("(",date, ")", sep=""), "\n"), 
          "Type 'demo(package='stocc')' to see a list of demos for this package.\n",
          "BE CAREFUL! The MCMC code can take a while to run if you start the demo.\n",
          "The raw code for the demos can be found by typing 'system.file('demo', package='stocc')'")
                      )

}


###
### Some Misc. functions for the future...
###

postMode <- function(x){
  out <- boxplot(x, plot=FALSE)$out
  x <- x[!x%in%out]
  dd <- density(x)
  return(dd$x[dd$y==max(dd$y)])
}

###
### Random effects formula...
###

### testing
# f <- ~ x*y + (1|s1) + (w|s1:s2)
# dat <- data.frame(x=rnorm(10),y=rgamma(10,1), w=rpois(10,1), s1=factor(rep(c(1:5),each=2)), s2=factor(rep(c(1,2),each=5)))

proc.form <- function(f){
  tms <- terms(f)
  tms.lab <- attr(tms, "term.labels")
  tms.lst <- strsplit(tms.lab, rep(" | ",length(tms.lab)), fixed=TRUE)
  fix.var <- attr(tms, "term.labels")[sapply(tms.lst, "length")==1]
  fix.model <- paste("~ ",paste(fix.var, collapse=" + "))
  re.lst <- tms.lst[sapply(tms.lst, "length")==2]
  if(length(re.lst)==0) re.model <- NULL
  else{
    re.model <- lapply(re.lst, function(x){list(model=paste("~",x[1]), sub=paste("~",x[2],"-1", collapse=""))})
  }
  return(list(fix.model=fix.model, re.model=re.model))
}

make.re.mat <- function(x){
  
}

