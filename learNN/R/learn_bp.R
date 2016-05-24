#' Learn Back Propagation
#'
#' @references \url{http://qua.st/handcoding-neural-network/}
#' @param X input data
#' @param y output data
#' @importFrom stats runif
#' @references \url{http://qua.st/handcoding-neural-network/}
#' \url{http://iamtrask.github.io/2015/07/12/basic-python-network/}
#' @export
#' @examples
#' # create data
#' X = matrix(c(0,0,1,
#'              0,1,1,
#'              1,0,1,
#'              1,1,1), nrow=4, byrow=TRUE)
#'
#' y = matrix(c(0,
#'              1,
#'              1,
#'              0),
#'              nrow=4)
#'
#' # run full function
#' learn_bp(X, y)


learn_bp <- function(X, y) {
  # no importing here

  nonlin = function(x,deriv=FALSE) {
    if(deriv==TRUE)
      return( x*(1-x) )

    return( 1/(1+exp(-x)) )
  }

  set.seed(1)

  # initialize weights randomly with mean 0
  syn0 = matrix(runif(n = 12, min=-1, max=1), nrow=3)
  syn1 = matrix(runif(n =  4, min=-1, max=1), nrow=4)

  for (j in 1:60000) {

    # Feed forward through layers 0, 1, and 2
    l0 = X
    l1 = nonlin(l0%*%syn0)
    l2 = nonlin(l1%*%syn1)

    # how much did we miss the target value?
    l2_error = y - l2

    if (j %% 10000 == 0)
      print(paste("Error:", mean(abs(l2_error))))

    # in what direction is the target value?
    # were we really sure? if so, don't change too much.
    l2_delta = l2_error*nonlin(l2,deriv=TRUE)

    # how much did each L1 value contribute to the error (according to the weights)?
    l1_error = l2_delta %*% t(syn1)

    # in what direction is the target l1?
    # were we really sure? if so, don't change too much.
    l1_delta = l1_error * nonlin(l1, deriv=TRUE)

    syn1 <<- syn1 + t(l1) %*% l2_delta
    syn0 <<- syn0 + t(l0) %*% l1_delta                     }
  print("Output After Training:")
  print(l1)
}


#' @name learn_bp11
#' @title Learn Back Propagation in 11 lines
#' @param X input data
#' @param y output data
#' @export
#' @references \url{http://qua.st/handcoding-neural-network/}
#' \url{http://iamtrask.github.io/2015/07/12/basic-python-network/}
#' @seealso \code{\link{learn_bp}}
#' @examples
#' # construct new data
#' X = matrix(c(0,0,1,0,1,1,1,0,1,1,1,1), nrow=4, byrow=TRUE)
#' y = matrix(c(0,1,1,0),nrow=4)
#'
#' # run 11 lines function
#' learn_bp11(X, y)
#'
#' # view output
#' syn0
#' syn1

learn_bp11 <- function(X, y) {
  syn0 = matrix(runif(n = 12, min=-1, max=1), nrow=3)
  syn1 = matrix(runif(n =  4, min=-1, max=1), nrow=4)
  for (j in 1:60000) {
    l1 = 1 / ( 1 + exp(-( X%*%syn0)) )
    l2 = 1 / ( 1 + exp(-(l1%*%syn1)) )
    l2_delta = (y-l2) * (l2*(1-l2))
    l1_delta = (l2_delta %*% t(syn1)) * (l1 * (1-l1))
    syn1 <<- syn1 + t(l1) %*% l2_delta
    syn0 <<- syn0 + t(X) %*% l1_delta                         }
}
