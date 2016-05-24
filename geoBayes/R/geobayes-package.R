##' Analysis of geostatistical data using Bayes and Empirical Bayes
##' methods.
##'  
##' This package provides functions to fit geostatistical data. The
##' data can be continuous, binary or count data and the models
##' implemented are flexible. Conjugate priors are assumed on some
##' parameters while inference on the other parameters can be done
##' through a full Bayesian analysis of by empirical Bayes methods.
##'
##' Some demonstration examples are provided. Type \code{demo(package
##' = "geoBayes")} to examine them.
##' @title The \code{geoBayes} package
##' @name geoBayes
##' @docType package
##' @useDynLib geoBayes
##' @author Evangelos Evangelou <e.evangelou@@maths.bath.ac.uk> and
##' Vivekananda Roy <vroy@@iastate.edu>
##' @seealso \code{geoR}, \code{geoRglm}
##' @examples \dontrun{
##' demo(package = "geoBayes")
##' demo(rhizoctonia3, package = "geoBayes")
##' }
##' @references Roy, V., Evangelou, E. and Zhu, Z. (2014). Empirical
##' Bayes methods for the transformed Gaussian random fields model
##' with additive measurement errors. In Upadhyay, S. K., Singh, U.,
##' Dey, D. K., and Loganathan, A., editors, \emph{Current Trends in
##' Bayesian Methodology with Applications}, Boca Raton, FL, USA, CRC
##' Press.
##'
##' Roy, V., Evangelou, E., and Zhu, Z. (2015). Efficient estimation
##' and prediction for the Bayesian spatial generalized linear mixed
##' model with flexible link functions. \emph{Biometrics}.
##' \url{http://dx.doi.org/10.1111/biom.12371}
##' @keywords package
NULL
