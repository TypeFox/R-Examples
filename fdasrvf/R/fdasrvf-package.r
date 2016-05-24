#' Elastic Functional Data Analysis
#'
#' A library for functional data analysis using the square root
#' velocity framework which performs pair-wise and group-wise
#' alignment as well as modeling using functional component
#' analysis
#'
#' @name fdasrvf
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2 [math.ST].
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @references J. D. Tucker, W. Wu, and A. Srivastava, ``Phase-Amplitude Separation of Proteomics Data Using Extended Fisher-Rao Metric," Electronic Journal of Statistics, Vol 8, no. 2. pp 1724-1733, 2014.
#' @references J. D. Tucker, W. Wu, and A. Srivastava, ``Analysis of signals under compositional noise With applications to SONAR data," IEEE Journal of Oceanic Engineering, Vol 29, no. 2. pp 318-330, Apr 2014.
#' @references Tucker, J. D. 2014, Functional Component Analysis and Regression using Elastic Methods. Ph.D. Thesis, Florida State University.
#' @references Robinson, D. T. 2012, Function Data Analysis and Partial Shape Matching in the Square Root Velocity Framework. Ph.D. Thesis, Florida State University.
#' @references Huang, W. 2014, Optimization Algorithms on Riemannian Manifolds with Applications. Ph.D. Thesis, Florida State University.
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @docType package
#' @useDynLib fdasrvf
#' @import foreach mvtnorm matrixcalc splines parallel doParallel Rcpp fields
#' @importFrom graphics layout legend matplot plot title
#' @importFrom stats approx cov optim predict quantile rnorm runif sd smooth.spline var spline
#' @aliases fdasrvf fdasrvf-package
NULL
#' Simulated two Gaussian Dataset
#'
#' A functional dataset where the individual functions are given by:
#' \eqn{y_i(t) = z_{i,1} e^{-(t-1.5)^2/2} + z_{i,2}e^{-(t+1.5)^2/2}}, \eqn{t \in [-3, 3], ~i=1,2,\dots, 21},
#' where \eqn{z_{i,1}} and \eqn{z_{i,2}} are \emph{i.i.d.} normal with mean one and standard deviation
#' 0.25. Each of these functions is then warped according to: \eqn{\gamma_i(t) = 6({e^{a_i(t+3)/6} -1 \over e^{a_i} - 1}) - 3}
#' if  \eqn{a_i \neq 0}, otherwise \eqn{\gamma_i = \gamma_{id}} (\eqn{gamma_{id}(t) = t})
#' is the identity warping). The variables are as follows: f containing the
#' 21 functions of 101 samples and time which describes the sampling
#'
#'
#' @docType data
#' @keywords datasets
#' @name simu_data
#' @usage data("simu_data")
#' @format A list which contains f and time
NULL
#' Distributed Gaussian Peak Dataset
#'
#' A functional dataset where the individual functions are given by a Gaussian
#' peak with locations along the \eqn{x}-axis. The variables are as follows:
#' f containing the 29 functions of 101 samples and time which describes the
#' sampling
#'
#'
#' @docType data
#' @keywords datasets
#' @name toy_data
#' @usage data("toy_data")
#' @format A list which contains f and time
NULL
#' Aligned Distributed Gaussian Peak Dataset
#'
#' A functional dataset where the individual functions are given by a Gaussian
#' peak with locations along the \eqn{x}-axis. The variables are as follows:
#' f containing the 29 functions of 101 samples and time which describes the
#' sampling which as been aligned
#'
#'
#' @docType data
#' @keywords datasets
#' @name toy_warp
#' @usage data("toy_warp")
#' @format A list which contains the outputs of the time_warping function
NULL
#' Aligned Simulated two Gaussian Dataset
#'
#' A functional dataset where the individual functions are given by:
#' \eqn{y_i(t) = z_{i,1} e^{-(t-1.5)^2/2} + z_{i,2}e^{-(t+1.5)^2/2}}, \eqn{t \in [-3, 3], ~i=1,2,\dots, 21},
#' where \eqn{z_{i,1}} and \eqn{z_{i,2}} are \emph{i.i.d.} normal with mean one and standard deviation
#' 0.25. Each of these functions is then warped according to: \eqn{\gamma_i(t) = 6({e^{a_i(t+3)/6} -1 \over e^{a_i} - 1}) - 3}
#' if  \eqn{a_i \neq 0}, otherwise \eqn{\gamma_i = \gamma_{id}} (\eqn{gamma_{id}(t) = t})
#' is the identity warping). The variables are as follows: f containing the
#' 21 functions of 101 samples and time which describes the sampling which has been aligned
#'
#'
#' @docType data
#' @keywords datasets
#' @name simu_warp
#' @usage data("simu_warp")
#' @format A list which contains the outputs of the time_warping function
NULL
#' MPEG7 Curve Dataset
#'
#' Contains the MPEG7 curve data set which is 20 curves in 65 classes. The array
#' is structured with dimension (2,100,65,20)
#'
#'
#' @docType data
#' @keywords datasets
#' @name beta
#' @usage data("mpeg7")
#' @format an array of shape (2,100,65,20)
NULL
