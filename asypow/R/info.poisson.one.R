info.poisson.one <- function(model="linear", theta, covariate) {
#-----------------------------------------------------------------------
#  Returns the information matrix for the poisson model with a
#  single covariate value
#
# model: One of {"linear", "quadratic"} Only enough to ensure a unique
#           match need be supplied.
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
      if (model == 1) u <- theta[1] + theta[2]*covariate else 
           u <- theta[1] + theta[2]*covariate + theta[3]*covariate^2

      if(model != (length(theta)-1))
                stop("theta inconsistant with model")

      lambda <- exp(u)

      d2lldl2 <- -1 / lambda

      dldu <- exp(u)

      dlda <- dldu
      dldb <- dldu * covariate
      if (model == 2) dldc <- dldu * covariate**2

      if (model == 1) hess <- matrix(c(dlda^2,rep(dlda*dldb,2),dldb^2),2,2)
      else hess <- matrix(c(dlda^2,dlda*dldb,dlda*dldc,
                     dldb*dlda,dldb^2,dldb*dldc,
                     dldc*dlda,dldc*dldb,dldc^2),3,3)

      return(- d2lldl2 * hess)
}
