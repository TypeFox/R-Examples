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
#   pcut   A vector of candidate cutoffs of pvalues.                   #
#                                                                      #
#   fold   fold-cross validation is performed.                         #
#                                                                      #
# Outputs :                                                            #
#                                                                      #
#   The total prediction errors for all lambda (row-wise) and          #
#   pcut (column-wise)                                                 #
#----------------------------------------------------------------------#
lambda.pcut.cv <- function(x, lambda, pcut, fold = 10L) {

  n <- nrow(x)

  cv <- 1L:n %% fold + 1L
  cv <- cbind(cv, sample(1:n)) 

  PE <- matrix(0.0,
               nrow = length(lambda),
               ncol = length(pcut),
               dimnames = list(lambda, pcut))

  k <- length(lambda)

  cat("cv:")
  for( i in 1L:fold ) {
    cat(" ", i)

    sets <- splitSets(cv = cv, i = i, x = x)

    PE <- PE + lambda.pcut.cv1(train = sets$train,
                               test = sets$test,
                               lambda = lambda, 
                               pcut = pcut)
  }
  cat("\n")
  return(PE / fold)
}
