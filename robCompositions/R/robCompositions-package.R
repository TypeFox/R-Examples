#' Robust Estimation for Compositional Data.
#' 
#' The package contains methods for imputation of compositional data including
#' robust methods, (robust) outlier detection for compositional data, (robust)
#' principal component analysis for compositional data, (robust) factor
#' analysis for compositional data, (robust) discriminant analysis (Fisher
#' rule) and (robust) Anderson-Darling normality tests for compositional data
#' as well as popular log-ratio transformations (alr, clr, ilr, and their
#' inverse transformations).
#' 
#' \tabular{ll}{ Package: \tab robCompositions\cr Type: \tab Package\cr
#' Version: \tab 1.3.3\cr Date: \tab 2009-11-28\cr License: \tab GPL 2\cr
#' LazyLoad: \tab yes\cr }
#' 
#' @name robCompositions-package
#' @aliases robCompositions-package robCompositions
#' @docType package
#' @author Matthias Templ, Peter Filzmoser, Karel Hron,
#' 
#' Maintainer: Matthias Templ <templ@@tuwien.ac.at>
#' @references Aitchison, J. (1986) \emph{The Statistical Analysis of
#' Compositional Data} Monographs on Statistics and Applied Probability.
#' Chapman \& Hall Ltd., London (UK). 416p. \
#' 
#' Filzmoser, P., and Hron, K. (2008) Outlier detection for compositional data
#' using robust methods. \emph{Math. Geosciences}, \bold{40} 233-248.
#' 
#' Filzmoser, P., Hron, K., Reimann, C. (2009) Principal Component Analysis for
#' Compositional Data with Outliers. \emph{Environmetrics}, \bold{20} (6),
#' 621--632.
#' 
#' P. Filzmoser, K. Hron, C. Reimann, R. Garrett (2009): Robust Factor Analysis
#' for Compositional Data.  \emph{Computers and Geosciences}, \bold{35} (9),
#' 1854--1861.
#' 
#' Hron, K. and Templ, M. and Filzmoser, P. (2010) Imputation of missing values
#' for compositional data using classical and robust methods
#' \emph{Computational Statistics and Data Analysis}, \bold{54} (12),
#' 3095--3107.
#' 
#' C. Reimann, P. Filzmoser, R.G. Garrett, and R. Dutter (2008): Statistical
#' Data Analysis Explained.  \emph{Applied Environmental Statistics with R}.
#' John Wiley and Sons, Chichester, 2008.
#' @keywords package
#' @examples
#' 
#' ## k nearest neighbor imputation
#' data(expenditures)
#' expenditures[1,3]
#' expenditures[1,3] <- NA
#' impKNNa(expenditures)$xImp[1,3]
#' 
#' ## iterative model based imputation
#' data(expenditures)
#' x <- expenditures
#' x[1,3]
#' x[1,3] <- NA
#' xi <- impCoda(x)$xImp
#' xi[1,3]
#' s1 <- sum(x[1,-3])
#' impS <- sum(xi[1,-3])
#' xi[,3] * s1/impS
#' 
#' xi <- impKNNa(expenditures)
#' xi
#' summary(xi)
#' \dontrun{plot(xi, which=1)}
#' plot(xi, which=2)
#' plot(xi, which=3)
#' 
#' ## pca
#' data(expenditures)
#' p1 <- pcaCoDa(expenditures)
#' p1
#' plot(p1)
#' 
#' ## outlier detection
#' data(expenditures)
#' oD <- outCoDa(expenditures)
#' oD
#' plot(oD)
#' 
#' ## transformations
#' data(arcticLake)
#' x <- arcticLake
#' x.alr <- addLR(x, 2)
#' y <- addLRinv(x.alr)
#' addLRinv(addLR(x, 3))
#' data(expenditures)
#' x <- expenditures
#' y <- addLRinv(addLR(x, 5))
#' head(x)
#' head(y)
#' addLRinv(x.alr, ivar=2, useClassInfo=FALSE)
#' 
#' data(expenditures)
#' eclr <- cenLR(expenditures)
#' inveclr <- cenLRinv(eclr)
#' head(expenditures)
#' head(inveclr)
#' head(cenLRinv(eclr$x.clr))
#' 
#' require(MASS)
#' Sigma <- matrix(c(5.05,4.95,4.95,5.05), ncol=2, byrow=TRUE)
#' z <- isomLRinv(mvrnorm(100, mu=c(0,2), Sigma=Sigma))
#' 
NULL







