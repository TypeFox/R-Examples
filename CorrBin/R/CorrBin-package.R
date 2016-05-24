

#'Nonparameterics for Correlated Binary and Multinomial Data
#'
#'This package implements nonparametric methods for analysing exchangeable
#'binary and multinomial data with variable cluster sizes with emphasis on trend testing. The
#'input should specify the treatment group, cluster-size, and the number of
#'responses (i.e. the number of cluster elements with the outcome of interest)
#'for each cluster.
#'
#'\tabular{ll}{ 
#' Package: \tab CorrBin\cr 
#' Type: \tab Package\cr 
#' Version: \tab 1.5\cr
#' Date: \tab 2014-12-18\cr 
#' License: \tab GPL 2\cr 
#' LazyLoad: \tab yes\cr
#'} 
#'\itemize{ \item The \code{\link{CBData}/\link{CMData}} and \code{\link{read.CBData}/\link{read.CMData}}
#'functions create a `CBData' or `CMData' object used by the analysis functions.  
#'\item \code{\link{ran.CBData}} and \code{\link{ran.CMData}} can be used to generate random 
#' binary or multinomial data using a variety of distributions.  
#'\item \code{\link{mc.test.chisq}} tests the assumption of marginal compatibility
#'underlying all the methods, while \code{\link{mc.est}} estimates the
#'distribution of the number of responses under marginal compatibility.  
#'\item Finally, \code{\link{trend.test}} performs three different tests for trend
#'along the treatment groups for binomial data. }
#'
#'@name CorrBin-package
#'@aliases CorrBin-package CorrBin
#'@docType package
#'@author Aniko Szabo
#'
#'Maintainer: Aniko Szabo <aszabo@@mcw.edu>
#'@references Szabo A, George EO. (2009) On the Use of Stochastic Ordering to
#'Test for Trend with Clustered Binary Data. \emph{Biometrika}
#'
#'Stefanescu, C. & Turnbull, B. W. (2003) Likelihood inference for exchangeable
#'binary data with varying cluster sizes. \emph{Biometrics}, 59, 18-24
#'
#'Pang, Z. & Kuk, A. (2007) Test of marginal compatibility and smoothing
#'methods for exchangeable binary data with unequal cluster sizes.
#'\emph{Biometrics}, 63, 218-227
#'@keywords package nonparametric
NULL
