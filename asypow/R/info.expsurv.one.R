info.expsurv.one <- function(model="linear", L, theta, covariate) {
#-----------------------------------------------------------------------
#  Returns the information matrix for the exponential survival model
#  with a single covariate value
#
# model: One of {"linear", "quadratic"} Only enough to ensure a unique
#           match need be supplied.
#
# L : Length of the experiment
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

      if(model != (length(theta)-1))
                stop("theta inconsistant with model")

      x <- covariate
      x2 <- x*x
      x3 <- x2*x
      x4 <- x3*x

      if (model == 1) u <- theta[1] + theta[2]*x  else 
         u <- theta[1] + theta[2]*x + theta[3]*x2

      w <- exp(u)
      wL <- w*L
      enwL <- exp(-wL)
      const <- ((enwL*(wL+2) + w*L - 2) - (enwL*(wL+1) - 1))/wL
      if (model == 1) 
         info <- const*matrix(c(1,x,x,x2),nrow=2,ncol=2) else
         info <- const*matrix(c(1,x,x2,x,x2,x3,x2,x3,x4),nrow=3,ncol=3)

      return(info)

}
