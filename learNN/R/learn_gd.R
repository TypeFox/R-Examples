#' Learn Gradient Descent
#'
#' @param X input data
#' @param y output data
#' @param alpha fraction of gradient descent
#' @param hiddenSize size of the hidden layer
#' @importFrom stats runif
#' @export
#' @references \url{http://qua.st/handcoding-gradient-descent/}
#' \url{http://iamtrask.github.io/2015/07/27/python-network-part2/}
#' @examples
#' # input dataset
#' X = matrix(c(0,0,1,
#'              0,1,1,
#'              1,0,1,
#'              1,1,1), nrow=4, byrow=TRUE)
#'
#' # output dataset
#' y = matrix(c(0,
#'              1,
#'              1,
#'              0), nrow=4)
#'
#' # set parameters
#' alpha = 0.1
#' hiddenSize = 32
#' # also try using:
#' # alphas = c(0.001,0.01,0.1,1,10,100,1000)
#' # for (alpha in alphas) {
#' #   print(paste("Training With Alpha", alpha))
#' #   learn_gd(X, y, alpha, hiddenSize)         }
#'
#' # run gradient descent function
#' learn_gd(X, y, alpha, hiddenSize)


learn_gd <- function(X, y, alpha, hiddenSize) {
  # no importing here

  # compute sigmoid nonlinearity
  sigmoid = function(x) {
    output = 1 / (1+exp(-x))
    return(output)            }

  signmoid_output_to_derivative = function(output)
    return( output*(1-output) )

  # randomly initialize our weights with mean 0
  set.seed(1)
  synapse_0 = matrix(runif(n = 3*hiddenSize, min=-1, max=1), nrow=3)
  synapse_1 = matrix(runif(n = hiddenSize,   min=-1, max=1), ncol=1)

  for (j in 1:60000) {

    # Feed forward through layers 0, 1, and 2
    layer_0 = X
    layer_1 = sigmoid(layer_0%*%synapse_0)
    layer_2 = sigmoid(layer_1%*%synapse_1)

    # how much did we miss the target value?
    layer_2_error = layer_2 - y

    if (j %% 10000 == 0)
      print(paste("Error after", j, "iterations:", mean(abs(layer_2_error))))

    # in what direction is the target value?
    # were we really sure? if so, don't change too much.
    layer_2_delta = layer_2_error * signmoid_output_to_derivative(layer_2)

    # how much did each l1 value contribute to the l2 error (according to the weights)?
    layer_1_error = layer_2_delta %*% t(synapse_1)

    # in what direction is the target l1?
    # were we really sure? if so, don't change too much.
    layer_1_delta = layer_1_error * signmoid_output_to_derivative(layer_1)

    syanpse_1 = synapse_1 - alpha * ( t(layer_1) %*% layer_2_delta )
    synapse_0 = synapse_0 - alpha * ( t(layer_0)%*%layer_1_delta   )                   }

  print("Output After Training (transposed):")
  print(t(layer_1))
}

#' @name learn_gd13
#' @title Learn Gradient Descent in 13 lines
#' @param X input data
#' @param y output data
#' @param alpha alpha to be used
#' @param hidden_dim dimension of the hidden layer
#' @importFrom stats runif
#' @export
#' @references \url{http://qua.st/handcoding-gradient-descent/}
#' \url{http://iamtrask.github.io/2015/07/27/python-network-part2/}
#' @seealso \code{\link{learn_gd}}
#' @examples
#' # create new data
#' alpha = 0.5
#' hidden_dim = 4
#' X = matrix(c(0,0,1,0,1,1,1,0,1,1,1,1), nrow=4, byrow=TRUE)
#' y = matrix(c(0,1,1,0),nrow=4)
#'
#' # run 13 lines function
#' learn_gd13(X, y, alpha, hidden_dim)

learn_gd13 <- function(X, y, alpha, hidden_dim) {
  # no importing here
  synapse_0 = matrix(runif(n = 3*hidden_dim, min=-1, max=1), nrow=3)
  synapse_1 = matrix(runif(n = hidden_dim, min=-1, max=1), ncol=1)
  for (j in 1:60000) {
    layer_1 = 1 / ( 1 + exp(-( X%*%synapse_0)) )
    layer_2 = 1 / ( 1 + exp(-(layer_1%*%synapse_1)) )
    layer_2_delta = (layer_2-y)*(layer_2*(1-layer_2))
    layer_1_delta = (layer_2_delta %*% t(synapse_1)) * (layer_1*(1-layer_1))
    synapse_1 <<- synapse_1 - alpha * ( t(layer_1) %*% layer_2_delta )
    synapse_0 <<- synapse_0 - alpha * ( t(X) %*% layer_1_delta )            }
}

#
