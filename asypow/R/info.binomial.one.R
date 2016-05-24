info.binomial.one <- function(model="linear", link="logistic",
                                                        theta, covariate) {
#-----------------------------------------------------------------------
#  Returns the information matrix for the binomial model with a
#  single covariate value
#
# model: One of {"linear", "quadratic"} Only enough to ensure a unique
#           match need be supplied.
#
# link: One of {"logistic", "complementary log"} Only enough to ensure
#           a unique match need be supplied.
#
# theta: Vector of parameters of the linear part of the model.
#
# covariate: Scalar value of the covariate.
#
#
# Returns: The contribution to the information matrix of a single
#           observation at value covariate.
#
#-----------------------------------------------------------------------

      model <- pmatch(model, c("linear","quadratic"))
      if (is.na(model)) stop("model must be one of {'linear','quadratic'}")
      if (model == 1) u <- theta[1] + theta[2]*covariate
      else u <- theta[1] + theta[2]*covariate + theta[3]*covariate^2

      if(model != (length(theta)-1))
                stop("theta inconsistant with model")

      link <- pmatch(link, c("logistic", "complementary log"))
      if (is.na(link))
                stop("link must be one of {'logistic','complementary log'}")
      if (link == 1) p <- exp(u) / (1 + exp(u)) else p <- 1 - exp(-exp(u))

      d2lldp2 <- -1 / (p *(1-p))

      if (link == 1) dpdu <- p * (1-p) else dpdu <- exp(u) * (1-p)

      dpda <- dpdu
      dpdb <- dpdu * covariate
      if (model == 2) dpdc <- dpdu * covariate^2

      if (model == 1) hess <- matrix(c(dpda^2,rep(dpda*dpdb,2),dpdb^2),2,2) else 
      hess <- matrix(c(dpda^2,dpda*dpdb,dpda*dpdc,
                     dpdb*dpda,dpdb^2,dpdb*dpdc,
                     dpdc*dpda,dpdc*dpdb,dpdc^2),3,3)

      return(- d2lldp2 * hess)
}
