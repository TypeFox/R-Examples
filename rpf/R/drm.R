##' Create a dichotomous response model
##'
##' For slope vector a, intercept c, pseudo-guessing parameter g,
##' upper bound u, and latent ability vector theta, the response probability
##' function is
##' \deqn{\mathrm P(\mathrm{pick}=0|a,c,g,u,\theta) = 1- \mathrm P(\mathrm{pick}=1|a,c,g,u,\theta)
##' }{P(pick=0|a,c,g,u,th) = 1-P(pick=1|a,c,g,u,th)}
##' \deqn{\mathrm P(\mathrm{pick}=1|a,c,g,u,\theta) = g+(u-g)\frac{1}{1+\exp(-(a\theta + c))}
##' }{P(pick=1|a,c,g,u,th) = g+(u-g)/(1+exp(-(a th + c)))}
##'
##' The pseudo-guessing and upper bound parameter are specified in
##' logit units (see \code{\link{logit}}).
##' 
##' For discussion on the choice of priors see Cai, Yang, and
##' Hansen (2011, p. 246).
##'
##' @param factors the number of factors
##' @param multidimensional whether to use a multidimensional model.
##' Defaults to \code{TRUE}.
##' @param poor if TRUE, use the traditional parameterization of
##' the 1d model instead of the slope-intercept parameterization
##' @return an item model
##' @export
##' @references Cai, L., Yang, J. S., & Hansen, M. (2011). Generalized
##' Full-Information Item Bifactor Analysis.  \emph{Psychological
##' Methods, 16}(3), 221-248.
##' @examples
##' spec <- rpf.drm()
##' rpf.prob(spec, rpf.rparam(spec), 0)
rpf.drm <- function(factors=1, multidimensional=TRUE, poor=FALSE) {
  if (!multidimensional && factors > 1) {
    stop("More than 1 dimension must use a multidimensional model")
  }
  m <- NULL
  id <- -1
  if (!multidimensional) {
    if (!poor) stop("The old parameterization is no longer available")
    id <- rpf.id_of("drm1-")
    m <- new("rpf.1dim.drm",
             outcomes=2,
             factors=1)
  } else {
    id <- rpf.id_of("drm")
    m <- new("rpf.mdim.drm",
             outcomes=2,
             factors=factors)
  }
  m@spec <- c(id, 2, m@factors)
  m
}

### 1dim

setMethod("rpf.rparam", signature(m="rpf.1dim.drm"),
          function(m, version) {
            n <- 1
            c(a=rlnorm(n, meanlog=0, sdlog=.5),
              b=rnorm(n),
              g=rbeta(1, 5,17),
              u=rbeta(1, 17,5))
          })

### mdim

setMethod("rpf.modify", signature(m="rpf.mdim.drm", factors="numeric"),
          function(m, factors) {
              rpf.drm(factors)
          })

##' Transform from [0,1] to the reals
##'
##' The logit function is a standard transformation from [0,1] (such
##' as a probability) to the real number line. This function is
##' exactly the same as qlogis.
##'
##' @param p a number between 0 and 1
##' @param location see qlogis
##' @param scale see qlogis
##' @param lower.tail see qlogis
##' @param log.p see qlogis
##' @examples
##' logit(.5)  # 0
##' logit(.25) # -1.098
##' logit(0)   # -Inf
##' @seealso
##' qlogis, plogis
logit <- qlogis

setMethod("rpf.rparam", signature(m="rpf.mdim.drm"),
          function(m, version) {
		  if (m@factors == 0) {
			  return(c(b=rnorm(1)))
		  }
		  if (version == 1L) {
			  c(a=rlnorm(m@factors, meanlog=0, sdlog=.5),
			    b=rnorm(1),
			    g=logit(rbeta(1, 5,17)),
			    u=logit(rbeta(1, 17,5)))
		  } else {
			  a <- rlnorm(m@factors, meanlog=0, sdlog=.5)
			  b <- rnorm(1) * sqrt(sum(a^2))
			  g <- logit(rbeta(1, 5,17))
			  u <- logit(rbeta(1, 17,5))
			  c(a=a, b=b, g=g, u=u)
		  }
          })

# Not sure if this is correct because of rotation
## as.loadings <- function(m, param) {
##   loading <- vector(mode="numeric", m@factors)
##   for (d in 1:m@factors) {
##     loading[d] <- param[d] / sqrt(1+sum(param[d:m@factors]^2))
##   }
##   loading
## }
