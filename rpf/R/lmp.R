##' Create logistic function of a monotonic polynomial (LMP) model
##'
##' This model is a dichotomous response model originally proposed by
##' Liang (2007) and is implemented using the parameterization by
##' Falk & Cai (in press).
##'
##' The LMP model replaces the linear predictor part of the
##' two-parameter logistic function with a monotonic polynomial,
##' \eqn{m(\theta,\omega,\xi,\mathbf{\alpha},\mathbf{\tau})}{m(theta, omega, alpha, tau)},
##'
##' \deqn{\mathrm P(\mathrm{pick}=1|\omega,\xi,\mathbf{\alpha},\mathbf{\tau},\theta)
##' = \frac{1}{1+\exp(-(\xi + m(\theta,\omega,\xi,\mathbf{\alpha},\mathbf{\tau})))}
##' }{P(pick=1|omega,xi,alpha,tau,th) = 1/(1+exp(-(xi + m(theta,omega,xi,alpha,tau))))}
##'
##' where \eqn{\mathbf{\alpha}}{alpha} and \eqn{\mathbf{\tau}}{tau} are vectors
##' of length k.
##'
##' The order of the polynomial is always odd and is controlled by
##' the user specified non-negative integer, k. The model contains
##' 2+2*k parameters and are used as input to the \code{rpf.prob}
##' function in the following order:
##' \eqn{\omega}{omega} - the natural log of the slope of the item model when k=0,
##' \eqn{\xi}{xi} - the intercept,
##' \eqn{\alpha}{alpha} and \eqn{\tau}{tau} - two parameters that control bends in
##' the polynomial. These latter parameters are repeated in the same order for
##' models with k>1. For example, a k=2 polynomial with have an item
##' parameter vector of: \eqn{\omega, \xi, \alpha_1, \tau_1, \alpha_2, \tau_2}{
##' omega, xi, alpha1, tau1, alpha2, tau2}.
##'
##' See Falk & Cai (in press) for more details as to how the
##' polynomial is constructed. In general, the polynomial looks like the
##' following, but coefficients, b, are not directly estimated, but
##' are a function of the item parameters.
##'
##' \deqn{\mathrm m(\theta) = \xi + b_1\theta + b_2\theta^2 + \dots + b_{2k+1}\theta^{2k+1}
##' }{m(theta) = xi + b_1*theta + b_2*theta^2 + \dots + b_(2k+1)*theta^{2k+1}}
##'
##' At the lowest order polynomial (k=0) the model reduces to the
##' two-parameter logistic (2PL) model. However, parameterization of the
##' slope parameter, \eqn{\omega}{omega}, is currently different than
##' the 2PL (i.e., slope = exp(\eqn{\omega}{omega})). This parameterization
##' ensures that the response function is always monotonically increasing
##' without requiring constrained optimization.
##'
##' Please note that the functions implementing this item model
##' may eventually be replaced or subsumed by an alternative
##' item model. That is, backwards compatability will not necessarily
##' be guaranteed and this item model should be considered experimental
##' until further notice.
##'
##' For example, Falk & Cai present a polytomous item model derived from
##' the generalized partial credit model that also uses a monotonic
##' polynomial as the linear predictor, referred to as a GPC-MP item
##' model. Since the GPC-MP reduces to the LMP when the number of
##' categories is 2, this is a potential candidate for replacing the
##' LMP item model. An alternative may include the retention of a
##' dichotomous response model, but with a lower (and upper)
##' asymptote that further reduces to the three-parameter logistic
##' (or four-parameter logistic) item model when k=0. Finally,
##' future versions may reparameterize \eqn{\omega}{omega}, or allow
##' the option to release constraints on monotonicity. For instance,
##' releasing constraints on \eqn{\omega}{omega} may be desirable
##' in cases where the user wishes to have the option of a
##' monotonically decreasing response function. Further releasing
##' constraints on \eqn{\tau}{tau} would allow nonmontonicity and would
##' be equivalent to replacing the linear predictor with a polynomial.
##'
##'
##' @param k a non-negative integer that controls the order of the
##' polynomial (2k+1) with a default of k=0 (1st order polynomial = 2PL).
##' @param multidimensional whether to use a multidimensional model.
##' Defaults to \code{FALSE}. The multidimensional version is not yet
##' available.
##' @return an item model
##' @export
##' @references Falk, C. F., & Cai, L. (in press). Maximum marginal likelihood
##' estimation of a monotonic polynomial generalized partial credit model with
##' applications to multiple group analysis. \emph{Psychometrika}.
##' \url{http://dx.doi.org/10.1007/s11336-014-9428-7}
##'
##' Liang (2007). \emph{A semi-parametric approach to estimating item response
##' functions}. Unpublished doctoral dissertation, Department of Psychology,
##' The Ohio State University.
##' @examples
##' spec <- rpf.lmp(1) # 3rd order polynomial
##' theta<-seq(-3,3,.1)
##' p<-rpf.prob(spec, c(-.11,.37,.24,-.21),theta)
##'
##' spec <- rpf.lmp(2) # 5th order polynomial
##' p<-rpf.prob(spec, c(.69,.71,-.5,-8.48,.52,-3.32),theta)

rpf.lmp <- function(k=0, multidimensional=FALSE) {
  if(!(k%%1==0)){
    stop("k must be an integer >= 0")
  }
  if(multidimensional){
      stop("Multidimensional LMP model is not yet supported")
  }
  m <- NULL
  id <- -1
  id <- rpf.id_of("lmp")
  m <- new("rpf.1dim.lmp",
           outcomes=2,
           factors=1)
  m@spec <- c(id, 2, m@factors, k)
  m
}

setMethod("rpf.rparam", signature(m="rpf.1dim.lmp"),
          function(m, version) {
            n <- 1
            k<-m$spec[4] ## ok to hardcode this index?
            ret<-c(omega=rnorm(n, 0, .5),xi=rnorm(n, 0, .75))
            if(k>0){
                for(i in 1:k){
                    ret<-c(ret,runif(n,-1,1),log(runif(n,.0001,1)))
                    names(ret)[(3+(i-1)*2):(2+(i*2))]<-c(paste("alpha",i,sep=""),paste("tau",i,sep=""))
                }
            }
            ret
        })
