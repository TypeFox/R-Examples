info.ordinal.one <- function(model="linear", link="logistic",
                                                    theta, covariate) {
#-----------------------------------------------------------------------
#  Returns the information matrix for the ordinal model with covariate
#
# model: One of {"linear", "quadratic"} Only enough to ensure a unique
#           match need be supplied.
#
# link: One of {"logistic", "complementary log"} Only enough to ensure
#           a unique match need be supplied.
#
# theta: Vector of parameters of the linear part of the model.
#          If there are K categories then
#    model = "linear" :
#      U[i] = THETA[i] + THETA[k]*covariate    i = 1,(k-1)
#    model = "quadratic" :
#      U[i] = THETA[i] + THETA[k]*covariate + THETA[k+1]*covariate^2
#                                                              i = 1,(k-1)
#                    
# covariate: Scalar value of the covariate.
#
#
# Returns: The contribution to the information matrix of a single
#           observation at value covariate.
#
#-----------------------------------------------------------------------

      imodel <- pmatch(model, c("linear","quadratic"))
      if (is.na(imodel)) stop("model must be one of {'linear','quadratic'}")
	
      link <- pmatch(link, c("logistic","complementary log"))
      if (is.na(link))
                stop("link must be one of {'logistic','complementary log'}")

	
      ntheta <- length(theta)
      ans <- matrix(0,ntheta,ntheta)
	
      if (imodel==1) ncat <- ntheta else 
             ncat <- ntheta - 1
      
      for (icat in 1:ncat) ans <- ans + info.ordinal.cat(
      imodel, link, theta, covariate, icat )

      return( ans )
}
