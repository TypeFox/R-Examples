#----------------------------------------------------------------------#
# Choose the tuning parameter of a ridge regression using              #
# cross-validation.                                                    #
#----------------------------------------------------------------------#
# Inputs :                                                             #
#                                                                      #
#  y      Length n response vector.                                    #
#                                                                      #
#  x      n x p matrix for covariates with p variables and             #
#         n sample size.                                               #
#                                                                      #
#  lambda A numeric vector for candidate tuning parameters for a       #
#         ridge regression.                                            #
#                                                                      #
#  fold   fold-cross validation used to choose the tuning parameter.   #
#                                                                      #
# Outputs :                                                            #
#                                                                      #
#  lambda The selected tuning parameter, which minimizes the           #
#         prediction error.                                            #
#                                                                      #
#   spe   The prediction error for all of the candidate lambda values. #
#----------------------------------------------------------------------#
ne.lambda.cv <- function(y, x, lambda, fold) {
    
  x <- cbind(y,x)

  n <- nrow(x)

  cv <- {1L:n} %% fold + 1L
  cv <- cbind(cv, sample(1:n))

  k <- length(lambda)
  spe <- numeric(k)

  for( i in 1L:fold ) {

    sets <- splitSets(cv = cv, i = i, x = x)

    coef <- svdFunc(x = sets$train[,-1L,drop = FALSE],
                    y = sets$train[,1L],
                    lambda = lambda)

    spe <- spe + colSums( {sets$test[,1L] - 
                           sets$test[,-1L,drop = FALSE] %*% coef}^2) 
  }

  return( list("lambda" = lambda,
               "spe" = spe / fold) )
}
