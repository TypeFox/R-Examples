##' rpf - Response Probability Functions
##'
##' The purpose of this package is to factor out logic and math common
##' to Item Factor Analysis fitting, diagnostics, and analysis.  It is
##' envisioned as core support code suitable for more specialized IFA
##' packages to build upon.
##'
##' This package provides optimized, low-level functions to map
##' parameters to response probabilities for dichotomous (1PL, 2PL and
##' 3PL) \code{\link{rpf.drm}} and polytomous (graded response
##' \code{\link{rpf.grm}}, partial credit/generalized partial credit
##' (via the nominal model), and nominal \code{\link{rpf.nrm}} items.
##'
##' Item model parameters are passed around as a numeric vector. A 1D
##' matrix is also acceptable. Regardless of model, parameters are
##' always ordered as follows: discrimination/slope ("a"),
##' difficulty/intercept ("b"), and pseudo guessing/upper-bound ("g"/"u"). If
##' person ability ranges from negative to positive then
##' probabilities are output from incorrect to correct. That is, a low
##' ability person (e.g., ability = -2) will be more likely to get an
##' item incorrect than correct. For example, a dichotomous model that
##' returns [.25, .75] indicates a probability of .25 for incorrect
##' and .75 for correct.  A polytomous model will have the most
##' incorrect probability at index 1 and the most correct probability
##' at the maximum index.
##'
##' All models are always in the logistic metric. To obtain normal
##' ogive discrimination parameters, divide slope parameters by
##' \code{\link{rpf.ogive}}. Item models are estimated in
##' slope-intercept form. Input/output matrices arranged in the way
##' most convenient for low-level processing in C. The maximum
##' absolute logit is 35 because f(x) := 1-exp(x) loses accuracy around f(-35)
##' and equals 1 at f(-38) due to the limited accuracy of double
##' precision floating point.
##'
##' This package could also accrete functions to support plotting (but
##' not the actual plot functions).
##'
##' @docType package
##' @rdname rpf.introduction
##' @name An introduction
##' @useDynLib rpf
##' @references Thissen, D. and Steinberg, L. (1986). A taxonomy of
##' item response models. \emph{Psychometrika 51}(4), 567-577.
##' @seealso
##' See \code{\link{rpf.rparam}} to create item parameters.
NULL

##' The base class for response probability functions.
##'
##' Item specifications should not be modified after creation.
##'
##' @name Class rpf.base
##' @rdname rpf.base-class
##' @aliases rpf.base-class
##' $,rpf.base-method
##' $<-,rpf.base-method
##' @export
setClass("rpf.base",
         representation(spec="numeric",
                        outcomes="numeric",
                        factors="numeric",
                        "VIRTUAL"))

imxExtractSlot <- function(x, name) {
	if (!.hasSlot(x, name)) {
		return(NULL)
	} else {
		return(slot(x, name))
	}
}

setMethod("$", "rpf.base", imxExtractSlot)

setReplaceMethod("$", "rpf.base", function(x, name, value) {
    stop("Slots are read-only")
})

##' The base class for 1 dimensional response probability functions.
##' @name Class rpf.1dim
##' @rdname rpf.1dim-class
##' @aliases rpf.1dim-class
##' @export
setClass("rpf.1dim", contains='rpf.base',
         representation("VIRTUAL"))

##' The base class for multi-dimensional response probability functions.
##' @name Class rpf.mdim
##' @rdname rpf.mdim-class
##' @aliases rpf.mdim-class
##' @export
setClass("rpf.mdim", contains='rpf.base',
         representation("VIRTUAL"))

##' Length of the item model vector
##' @param m item model
##' @aliases
##' rpf.numSpec,rpf.base-method
##' rpf_numSpec_wrapper
##' @examples
##' rpf.numSpec(rpf.grm(outcomes=3))
##' rpf.numSpec(rpf.nrm(outcomes=3))
setGeneric("rpf.numSpec", function(m) standardGeneric("rpf.numSpec"))

setMethod("rpf.numSpec", signature(m="rpf.base"),
          function(m) {
            if (length(m@spec)==0) {
              stop("Not implemented")
            } else {
              .Call(rpf_numSpec_wrapper, m@spec)
            }
          })

##' Length of the item parameter vector
##'
##' @param m item model
##' @aliases
##' rpf.numParam,rpf.base-method
##' rpf_numParam_wrapper
##' @examples
##' rpf.numParam(rpf.grm(outcomes=3))
##' rpf.numParam(rpf.nrm(outcomes=3))
setGeneric("rpf.numParam", function(m) standardGeneric("rpf.numParam"))

setMethod("rpf.numParam", signature(m="rpf.base"),
          function(m) {
            .Call(rpf_numParam_wrapper, m@spec)
          })

##' Create a similar item specification with the given number of factors
##'
##' @param m item model
##' @param factors the number of factors/dimensions
##' @aliases
##' rpf.modify,rpf.mdim.drm,numeric-method
##' rpf.modify,rpf.mdim.graded,numeric-method
##' rpf.modify,rpf.mdim.nrm,numeric-method
##' @examples
##' s1 <- rpf.grm(factors=3)
##' rpf.rparam(s1)
##' s2 <- rpf.modify(s1, 1)
##' rpf.rparam(s2)
setGeneric("rpf.modify", function(m, factors) standardGeneric("rpf.modify"))

##' Retrieve a description of the given parameter
##' @param m item model
##' @param num vector of parameters (defaults to all)
##' @return a list containing the type, upper bound, and lower bound
##' @aliases
##' rpf.paramInfo,rpf.base-method
##' rpf_paramInfo_wrapper
##' @examples
##' rpf.paramInfo(rpf.drm())
setGeneric("rpf.paramInfo", function(m, num=NULL) standardGeneric("rpf.paramInfo"))

setMethod("rpf.paramInfo", signature(m="rpf.base"),
          function(m, num=NULL) {
            if (missing(num)) {
              num <- 1:rpf.numParam(m)
            }
            if (length(num) == 0) {
              stop("Which parameter?")
            } else if (length(num) == 1) {
              .Call(rpf_paramInfo_wrapper, m@spec, num-1)
            } else {
              sapply(num, function (px) .Call(rpf_paramInfo_wrapper, m@spec, px-1))
            }
          })

##' Item parameter derivatives
##'
##' Evaluate the partial derivatives of the log likelihood with
##' respect to each parameter at \code{where} with \code{weight}.
##'
##' It is not easy to write an example for this function. To evaluate
##' the derivative, you need to sum the derivatives across a
##' quadrature. You also need response outcome weights at each
##' quadrature point. It is not anticipated that this function will be
##' often used in R code. It's mainly to expose a C-level function for
##' occasional debugging.
##'
##' @param m item model
##' @param param item parameters
##' @param where location in the latent space
##' @param weight per outcome weights (typically derived by observation)
##' @return first and second order partial derivatives of the log
##' likelihood evaluated at \code{where}. For p parameters, the first
##' p values are the first derivative and the next p(p+1)/2 columns
##' are the lower triangle of the second derivative.
##' @seealso
##' The numDeriv package.
##' @aliases
##' rpf.dLL,rpf.base,numeric,numeric,numeric-method
##' rpf.dLL,rpf.base,numeric,NULL,numeric-method
##' rpf_dLL_wrapper
setGeneric("rpf.dLL", function(m, param, where, weight) standardGeneric("rpf.dLL"))

setMethod("rpf.dLL", signature(m="rpf.base", param="numeric",
                                 where="numeric", weight="numeric"),
          function(m, param, where, weight) {
            if (length(m@spec)==0) {
              stop("Not implemented")
            } else {
              .Call(rpf_dLL_wrapper, m@spec, param, where, weight)
            }
          })

setMethod("rpf.dLL", signature(m="rpf.base", param="numeric",
                                 where="NULL", weight="numeric"),
          function(m, param, where, weight) {
            if (length(m@spec)==0) {
              stop("Not implemented")
            } else {
              .Call(rpf_dLL_wrapper, m@spec, param, where, weight)
            }
          })

##' Item derivatives with respect to the location in the latent space
##'
##' Evaluate the partial derivatives of the response probability with
##' respect to ability. See \link{rpf.info} for an application.
##'
##' @param m item model
##' @param param item parameters
##' @param where location in the latent distribution
##' @param dir if more than 1 factor, a basis vector]
##' @aliases
##' rpf_dTheta_wrapper
##' rpf.dTheta,rpf.base,numeric,numeric,numeric-method
##' rpf.dTheta,rpf.base,numeric,matrix,numeric-method
setGeneric("rpf.dTheta", function(m, param, where, dir) standardGeneric("rpf.dTheta"))

setMethod("rpf.dTheta", signature(m="rpf.base", param="numeric",
                                  where="numeric", dir="numeric"),
          function(m, param, where, dir) {
            if (length(m@spec)==0) {
              stop("Not implemented")
            } else {
              .Call(rpf_dTheta_wrapper, m@spec, param, where, dir)
            }
          })

setMethod("rpf.dTheta", signature(m="rpf.base", param="numeric",
                                  where="matrix", dir="numeric"),
          function(m, param, where, dir) {
            dP.raw <- apply(where, 2, function(w) rpf.dTheta(m, param, w, dir))
            list(gradient=sapply(dP.raw, function(deriv) deriv$gradient),
                 hessian=sapply(dP.raw, function(deriv) deriv$hessian))
          })

##' Rescale item parameters
##'
##' Adjust item parameters for changes in mean and covariance of the
##' latent distribution.
##'
##' @param m item model
##' @param param item parameters
##' @param mean vector of means
##' @param cov covariance matrix
##' @aliases
##' rpf_rescale_wrapper
##' rpf.rescale,rpf.base,numeric,numeric,matrix-method
##' @examples
##' spec <- rpf.grm()
##' p1 <- rpf.rparam(spec)
##' testPoint <- rnorm(1)
##' move <- rnorm(1)
##' cov <- as.matrix(rlnorm(1))
##' Icov <- solve(cov)
##' padj <- rpf.rescale(spec, p1, move, cov)
##' pr1 <- rpf.prob(spec, padj, (testPoint-move) %*% Icov)
##' pr2 <- rpf.prob(spec, p1, testPoint)
##' abs(pr1 - pr2) < 1e9
setGeneric("rpf.rescale", function(m, param, mean, cov) standardGeneric("rpf.rescale"))

setMethod("rpf.rescale", signature(m="rpf.base", param="numeric",
                                   mean="numeric", cov="matrix"),
          function(m, param, mean, cov) {
            if (length(m@spec)==0) {
              stop("Not implemented")
            } else {
              .Call(rpf_rescale_wrapper, m@spec, param, mean, cov)
            }
          })

##' Map an item model, item parameters, and person trait score into a
##' probability vector
##'
##' @param m an item model
##' @param param item parameters
##' @param theta the trait score(s)
##' @return a vector of probabilities. For dichotomous items,
##' probabilities are returned in the order incorrect, correct.
##' Although redundent, both incorrect and correct probabilities are
##' returned in the dichotomous case for API consistency with
##' polytomous item models.
##' @docType methods
##' @aliases
##' rpf.prob,rpf.1dim,numeric,numeric-method
##' rpf.prob,rpf.mdim,numeric,NULL-method
##' rpf.prob,rpf.mdim,numeric,numeric-method
##' rpf.prob,rpf.mdim,numeric,matrix-method
##' rpf.prob,rpf.base,data.frame,numeric-method
##' rpf.prob,rpf.base,matrix,numeric-method
##' rpf.prob,rpf.base,matrix,matrix-method
##' rpf.prob,rpf.1dim,numeric,matrix-method
##' rpf.prob,rpf.1dim.grm,numeric,numeric-method
##' rpf.prob,rpf.mdim.grm,numeric,numeric-method
##' rpf.prob,rpf.mdim.nrm,numeric,matrix-method
##' rpf.prob,rpf.mdim.mcm,numeric,matrix-method
##' rpf.prob,rpf.mdim.grm,numeric,matrix-method
##' rpf_prob_wrapper
##' @export
##' @examples
##' i1 <- rpf.drm()
##' i1.p <- rpf.rparam(i1)
##' rpf.prob(i1, c(i1.p), -1)   # low trait score
##' rpf.prob(i1, c(i1.p), c(0,1))    # average and high trait score
setGeneric("rpf.prob", function(m, param, theta) standardGeneric("rpf.prob"))

##' Map an item model, item parameters, and person trait score into a
##' probability vector
##'
##' Note that in general, exp(rpf.logprob(..)) != rpf.prob(..) because
##' the range of logits is much wider than the range of probabilities
##' due to limitations of floating point numerical precision.
##'
##' @param m an item model
##' @param param item parameters
##' @param theta the trait score(s)
##' @return a vector of probabilities. For dichotomous items,
##' probabilities are returned in the order incorrect, correct.
##' Although redundent, both incorrect and correct probabilities are
##' returned in the dichotomous case for API consistency with
##' polytomous item models.
##' @docType methods
##' @aliases
##' rpf.logprob,rpf.1dim,numeric,numeric-method
##' rpf.logprob,rpf.1dim,numeric,matrix-method
##' rpf.logprob,rpf.mdim,numeric,matrix-method
##' rpf.logprob,rpf.mdim,numeric,numeric-method
##' rpf.logprob,rpf.mdim,numeric,NULL-method
##' rpf_logprob_wrapper
##' @export
##' @examples
##' i1 <- rpf.drm()
##' i1.p <- rpf.rparam(i1)
##' rpf.logprob(i1, c(i1.p), -1)   # low trait score
##' rpf.logprob(i1, c(i1.p), c(0,1))    # average and high trait score
setGeneric("rpf.logprob", function(m, param, theta) standardGeneric("rpf.logprob"))

setMethod("rpf.logprob", signature(m="rpf.1dim", param="numeric", theta="numeric"),
          function(m, param, theta) {
            if (length(m@spec)==0) {
              stop("Not implemented")
            } else {
              .Call(rpf_logprob_wrapper, m@spec, param, theta)
            }
          })

setMethod("rpf.logprob", signature(m="rpf.mdim", param="numeric", theta="matrix"),
          function(m, param, theta) {
            if (length(m@spec)==0) {
              stop("Not implemented")
            } else {
              .Call(rpf_logprob_wrapper, m@spec, param, theta)
            }
          })

setMethod("rpf.logprob", signature(m="rpf.mdim", param="numeric", theta="numeric"),
          function(m, param, theta) {
            rpf.logprob(m, param, as.matrix(theta))
          })

setMethod("rpf.logprob", signature(m="rpf.mdim", param="numeric", theta="NULL"),
          function(m, param, theta) {
            if (length(m@spec)==0) {
              stop("Not implemented")
            } else {
              .Call(rpf_logprob_wrapper, m@spec, param, theta)
            }
          })

setMethod("rpf.logprob", signature(m="rpf.1dim", param="numeric", theta="matrix"),
          function(m, param, theta) {
            rpf.logprob(m, param, as.numeric(theta))
          })

setMethod("rpf.prob", signature(m="rpf.1dim", param="numeric", theta="numeric"),
          function(m, param, theta) {
            if (length(m@spec)==0) {
              exp(rpf.logprob(m, param, theta))
            } else {
              .Call(rpf_prob_wrapper, m@spec, param, theta)
            }
          })

setMethod("rpf.prob", signature(m="rpf.mdim", param="numeric", theta="NULL"),
          function(m, param, theta) {
            if (length(m@spec)==0) {
              exp(rpf.logprob(m, param, theta))
            } else {
              .Call(rpf_prob_wrapper, m@spec, param, theta)
            }
          })

setMethod("rpf.prob", signature(m="rpf.mdim", param="numeric", theta="matrix"),
          function(m, param, theta) {
            if (length(m@spec)==0) {
              exp(rpf.logprob(m, param, theta))
            } else {
              .Call(rpf_prob_wrapper, m@spec, param, theta)
            }
          })

setMethod("rpf.prob", signature(m="rpf.base", param="data.frame", theta="numeric"),
          function(m, param, theta) {
            exp(rpf.logprob(m, as.numeric(param), theta))
          })

setMethod("rpf.prob", signature(m="rpf.base", param="matrix", theta="numeric"),
          function(m, param, theta) {
            rpf.prob(m, as.numeric(param), theta)
          })

setMethod("rpf.prob", signature(m="rpf.base", param="matrix", theta="matrix"),
          function(m, param, theta) {
            rpf.prob(m, as.numeric(param), theta)
          })

setMethod("rpf.prob", signature(m="rpf.1dim", param="numeric", theta="matrix"),
          function(m, param, theta) {
            rpf.prob(m, param, as.numeric(theta))
          })

setMethod("rpf.prob", signature(m="rpf.mdim", param="numeric", theta="numeric"),
          function(m, param, theta) {
            rpf.prob(m, param, as.matrix(theta))
          })

##' Map an item model, item parameters, and person trait score into a
##' information vector
##'
##' @param ii an item model
##' @param ii.p item parameters
##' @param where the location in the latent distribution
##' @param basis if more than 1 factor, a positive basis vector
##' @return Fisher information
##' @export
##' @examples
##' i1 <- rpf.drm()
##' i1.p <- c(.6,1,.1,.95)
##' theta <- seq(0,3,.05)
##' plot(theta, rpf.info(i1, i1.p, t(theta)), type="l")
##' @references Dodd, B. G., De Ayala, R. J. & Koch,
##' W. R. (1995). Computerized adaptive testing with polytomous items.
##' \emph{Applied psychological measurement 19}(1), 5-22.
rpf.info <- function(ii, ii.p, where, basis=1) {
  if (any(basis < 0)) warning("All components of the basis vector should be positive")
  if (!missing(basis)) {
    basis <- basis/sqrt(sum(basis^2))
  }
  P <- rpf.prob(ii, ii.p, where)
  dP <- rpf.dTheta(ii, ii.p, where, basis)
  colSums(dP$gradient^2 / P - dP$hessian)
}

##' Generates item parameters
##'
##' This function generates random item parameters. The version
##' argument is available if you are writing a test that depends on
##' reproducable random parameters (using \code{set.seed}).
##'
##' @param m an item model
##' @param version the version of random parameters
##' @return item parameters
##' @docType methods
##' @aliases
##' rpf.rparam,rpf.1dim.drm-method
##' rpf.rparam,rpf.mdim.drm-method
##' rpf.rparam,rpf.1dim.graded-method
##' rpf.rparam,rpf.mdim.graded-method
##' rpf.rparam,rpf.mdim.nrm-method
##' rpf.rparam,rpf.mdim.mcm-method
##' rpf.rparam,rpf.1dim.lmp-method
##' @export
##' @examples
##' i1 <- rpf.drm()
##' rpf.rparam(i1)
setGeneric("rpf.rparam", function(m, version=2L) standardGeneric("rpf.rparam"))

##' The ogive constant
##'
##' The ogive constant can be multiplied by the discrimination
##' parameter to obtain a response curve very similar to the Normal
##' cumulative distribution function (Haley, 1952; Molenaar, 1974).
##' Recently, Savalei (2006) proposed a new constant of 1.749 based on
##' Kullback-Leibler information.
##'
##' In recent years, the logistic has grown in favor, and therefore,
##' this package does not offer any special support for this
##' transformation (Baker & Kim, 2004, pp. 14-18).
##'
##' @export
##' @references Camilli, G. (1994). Teacher's corner: Origin of the
##' scaling constant d=1.7 in Item Response Theory. \emph{Journal of
##' Educational and Behavioral Statistics, 19}(3), 293-295.
##'
##' Baker & Kim (2004). \emph{Item Response Theory: Parameter
##' Estimation Techniques.} Marcel Dekker, Inc.
##'
##' Haley, D. C. (1952). \emph{Estimation of the dosage mortality
##' relationship when the dose is subject to error} (Technical Report
##' No. 15). Stanford University Applied Mathematics and Statistics
##' Laboratory, Stanford, CA.
##'
##' Molenaar, W. (1974). De logistische en de normale kromme [The
##' logistic and the normal curve]. \emph{Nederlands Tijdschrift voor de
##' Psychologie} 29, 415-420.
##'
##' Savalei, V. (2006). Logistic approximation to the normal: The KL
##' rationale. \emph{Psychometrika, 71}(4), 763--767.
rpf.ogive <- 1.702

##' The base class for 1 dimensional graded response probability functions.
##'
##' This class contains methods common to both the generalized partial
##' credit model and the graded response model.
##'
##' @name Class rpf.1dim.graded
##' @rdname rpf.1dim.graded-class
##' @aliases rpf.1dim.graded-class
##' @export
setClass("rpf.1dim.graded", contains='rpf.1dim',
         representation("VIRTUAL"))

##' The base class for multi-dimensional graded response probability
##' functions.
##'
##' This class contains methods common to both the generalized partial
##' credit model and the graded response model.
##'
##' @name Class rpf.mdim.graded
##' @rdname rpf.mdim.graded-class
##' @aliases rpf.mdim.graded-class
##' @export
setClass("rpf.mdim.graded", contains='rpf.mdim',
         representation("VIRTUAL"))

##' The unidimensional graded response item model.
##'
##' @export
##' @name Class rpf.1dim.grm
##' @rdname rpf.1dim.grm-class
##' @aliases rpf.1dim.grm-class
setClass("rpf.1dim.grm", contains='rpf.1dim.graded')

##' Unidimensional dichotomous item models (1PL, 2PL, and 3PL).
##'
##' @export
##' @name Class rpf.1dim.drm
##' @rdname rpf.1dim.drm-class
##' @aliases rpf.1dim.drm-class
setClass("rpf.1dim.drm", contains='rpf.1dim')

##' Multidimensional dichotomous item models (M1PL, M2PL, and M3PL).
##'
##' @export
##' @name Class rpf.mdim.drm
##' @rdname rpf.mdim.drm-class
##' @aliases rpf.mdim.drm-class
setClass("rpf.mdim.drm", contains='rpf.mdim')

##' The multidimensional graded response item model.
##'
##' @export
##' @name Class rpf.mdim.grm
##' @rdname rpf.mdim.grm-class
##' @aliases rpf.mdim.grm-class
setClass("rpf.mdim.grm", contains='rpf.mdim.graded')

##' The nominal response item model (both unidimensional and
##' multidimensional models have the same parameterization).
##'
##' @export
##' @name Class rpf.mdim.nrm
##' @rdname rpf.mdim.nrm-class
##' @aliases rpf.mdim.nrm-class
setClass("rpf.mdim.nrm", contains='rpf.mdim')

##' The multiple-choice response item model (both unidimensional and
##' multidimensional models have the same parameterization).
##'
##' @export
##' @name Class rpf.mdim.mcm
##' @rdname rpf.mdim.mcm-class
##' @aliases rpf.mdim.mcm-class
setClass("rpf.mdim.mcm", contains='rpf.mdim')

##' Unidimensional logistic function of a monotonic polynomial.
##'
##' @export
##' @name Class rpf.1dim.lmp
##' @rdname rpf.1dim.lmp-class
##' @aliases rpf.1dim.lmp-class
##'
setClass("rpf.1dim.lmp", contains='rpf.1dim')

##' Convert an rpf item model name to an ID
##'
##' This is an internal function and should not be used.
##'
##' @param name name of the item model (string)
##' @return the integer ID assigned to the given model
rpf.id_of <- function(name) {
   .Call(get_model_names, name)
}
