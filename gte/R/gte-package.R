#' Generalized Turnbull's Estimator
#' 
#' The \code{gte} function computes the generalized Turnbull's estimator proposed by Dehghan and Duchesne (2011).
#' It is a nonparametric estimator of a conditional survival function given a vector of continuous covariates 
#' that can handle interval-censored lifetimes.   
#' 
#' \tabular{ll}{ 
#' Package: \tab gte\cr 
#' Type: \tab Package\cr 
#' Version: \tab 1.2-2\cr 
#' Date: \tab 2015-02-25\cr 
#' License: \tab GPL-2\cr 
#' }
#' 
#' @name gte-package
#' @aliases gte-package
#' @docType package
#' @author Mohammad Hossein Dehghan, Thierry Duchesne and Sophie Baillargeon
#' 
#' Maintainer: Thierry Duchesne <thierry.duchesne@@mat.ulaval.ca>
#' @references Dehghan, M. H. and Duchesne, T. (2011). A generalization of Turnbull's estimator for 
#' nonparametric estimation of the conditional survival function with interval-censored data.
#' \emph{Lifetime Data Analysis}, \bold{17}, 234-255.
#' @keywords package
NULL

#' Simulated Data
#' 
#' Simulated Interval-censored data
#' 
#' The value \code{R = NA} means that the observation is right censored (occurs 2 times). 
#' If \code{L = NA}, then the observation is left censored (occurs 26 times).
#' An observation with \code{R = L} means that the time of occurence of the event
#' is known exactly (occurs 3 times).  
#' 
#' @name simul
#' @docType data
#' @format A data frame with 100 observations on the following 3 variables.
#' \describe{ 
#'   \item{L}{the left endpoints of the censoring interval}
#'   \item{R}{the right endpoints of the censoring interval} 
#'   \item{Z}{a continuous covariate} 
#' }
#' @references Dehghan, M. H. and Duchesne, T. (2011). A generalization of Turnbull's estimator for 
#' nonparametric estimation of the conditional survival function with interval-censored data.
#' \emph{Lifetime Data Analysis}, \bold{17}, 234-255.
#' @keywords datasets
#' @examples
#' data(simul)
#' Fit <- gte(Surv(L, R, type="interval2") ~ Z, data=simul, z=15)
#' plot(Fit)
NULL
