#----------------------------------------------------------------------#
# Chooses the tuning parameter of the ridge inverse by minimizing      #
# cross validation estimates of the total prediction errors of the p   #
# separate ridge regressions.                                          #
#----------------------------------------------------------------------#
# Inputs :                                                             #
#                                                                      #
#   x      An n by p data matrix.                                      #
#                                                                      #
#   lambda A numeric vector of candidate tuning parameters.            #
#                                                                      #
#   fold   fold-cross validation is performed.                         #
#                                                                      #
# Outputs :                                                            #
#                                                                      #
#   lambda The selected tuning parameter, which minimizes the total    #
#          prediction errors.                                          #
#                                                                      #
#   spe    The total prediction error for all the candidate lambda     #
#          values.                                                     #
#----------------------------------------------------------------------#
lambda.cv <- function(x, lambda, fold) {

  n <- nrow(x)
  p <- ncol(x)

  cv <- {1L:n} %% fold + 1L
  cv <- cbind(cv,sample(1:n))

  k <- length(lambda)

  spe <- numeric(k)

  for( i in 1L:fold ) {

    sets <- splitSets(cv = cv, i = i, x = x)

    speMat <- matrix(data = 0.0, nrow = p, ncol = k)

    for( j in 1L:p ) {

      coef <- svdFunc(x = sets$train[,-j,drop = FALSE],
                      y = sets$train[,j],
                      lambda = lambda)

      speMat[j,] <- colSums( {sets$test[,j] - 
                              sets$test[,-j,drop=FALSE] %*% coef}^2 )
    }

    spe <- spe + colSums(speMat)

  }

  return(list("lambda" = lambda,
              "spe" = spe / {fold * p}))
}
