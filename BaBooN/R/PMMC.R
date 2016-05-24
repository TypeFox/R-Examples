# Wrapper for Rcpp invoke of Predictive Mean Matching function of BBPMMrow
# Version:             0.2
# Date:         2014-12-17
#	Author:             T.S.
# Note:     Relies on Rcpp / RcppArmadillo
# Further infos, references and credits:
#  See for Rcpp: Eddelbuettel, D. and Francois, R. (2011) Rcpp: Seamless R and C++ Integration.
#                Journal of Statistical Software, Vol. 40, No. 8, pp. 1--18. URL http://www.jstatsoft.org/v40/i08/.
#
#  See for RcppArmadillo: Eddelbuettel, D. and Sanderson, C. (2014) RcppArmadillo: Accelerating R with high-performance C++ linear algebra.
#                         Computational Statistics and Data Analysis, Vol. 71, March 2014, pp. 1054--1063.
# License:  (GPL >= 2)

PMMC <- function(...){
  .Call( "PMMC", ..., PACKAGE = "BaBooN" )
}