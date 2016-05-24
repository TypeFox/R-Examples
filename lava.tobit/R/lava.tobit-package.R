##' Estimation and simulation of probit and tobit latent variable models
##'  
##' Framwork for estimating parameters and simulate data from Latent Variable
##' Models with binary and censored observations. Plugin for the \code{lava}
##' package
##' 
##' \tabular{ll}{ Package: \tab lava.tobit \cr Type: \tab Package \cr Version:
##' \tab 0.4-5 \cr Date: \tab 2012-03-15 \cr License: \tab GPL-3 \cr LazyLoad:
##' \tab yes \cr }
##' 
##' @name lava.tobit
##' @aliases lava.tobit lava.tobit-package
##' @docType package
##' @author Klaus K. Holst Maintainer: <kkho@@biostat.ku.dk>
##' @keywords package
##' @examples
##' 
##' m <- lvm(list(c(y,z) ~ x, y~z))
##' ## Simulate 200 observation from path analysis model
##' ## with all slopes and residual variances set to 1 and intercepts 0:
##' d <- sim(m,200)
##' ## Dichotomize y and introduce censoring on z
##' d <- transform(d, y=as.factor(y>0), z=Surv(z,z<2))
##' \donttest{
##' e <- estimate(m,d,control=list(trace=1))
##' effects(e,y~x)
##' }
##' 
NULL

##' For internal use
##'
##' @title For internal use
##' @name lava.tobit.estimate.hook
##' @rdname internal
##' @aliases tobit_gradient.lvm tobit_hessian.lvm tobit_logLik.lvm
##' tobit_method.lvm tobit_objective.lvm tobitw_gradient.lvm tobitw_hessian.lvm
##' tobitw_method.lvm lava.tobit.color.hook lava.tobit.estimate.hook
##' lava.tobit.init.hook lava.tobit.sim.hook
##' @author Klaus K. Holst
##' @keywords utilities
##' @export
NULL
