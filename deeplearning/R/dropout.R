#' Generates dropout masks for dnn
#'
#' This function generates dropout maks for dnn
#' @param darch, a DArch instance
#' @param dropout_input, the dropout rate for the input layer
#' @param dropout_hidden, the dropout rate for the hidden layer
#'
#' @importFrom darch getLayers getLayerWeights
#' @references Dropout: A Simple Way to Prevent Neural Networks from
#'  Overfitting, Nitish Srivastava
#' @seealso \url{https://www.cs.toronto.edu/~hinton/absps/JMLRdropout.pdf}



generateDropoutMasksForDarch <- function(darch, dropout_input, dropout_hidden)
{
  dropoutMasks <- list()
  numLayers <- length(getLayers(darch))
  # generate dropout masks
  darch@dropoutMasks[[1]]<-
    generateDropoutMask(nrow(getLayerWeights(darch, 1)[]) - 1,
                        dropout_input)

  for (i in 1:(numLayers - 1))
  {
    darch@dropoutMasks[[i + 1]] <-
      generateDropoutMask(nrow(getLayerWeights(darch, i+1)[])-1,
                          dropout_hidden)
  }

  return (darch)
}

#' Generates the dropout mask for the deep neural network
#'
#' This function generates the dropout mask for the deep neural network
#' @param length, the dimension of the layer
#' @param dropoutRate, the dropout rate
#'
#' @references Dropout: A Simple Way to Prevent Neural Networks from
#'  Overfitting, Nitish Srivastava
#' @seealso \url{https://www.cs.toronto.edu/~hinton/absps/JMLRdropout.pdf}



generateDropoutMask <- function(length, dropoutRate)
{
  if (dropoutRate == 0)
  {
    ret <- rep(1, length)
  }
  else
  {
    ret <- sample(c(0, 1/(1 - dropoutRate)), length, replace = T,
                  prob = c(dropoutRate, 1 - dropoutRate))
  }

  return (ret)
}



#' Applies the given dropout mask to the given data row-wise.
#'
#' This function multiplies each row with the dropout mask. To apply the dropout
#' mask by row, it can simply be multiplied with the data matrix. This does not
#' work of the mask is to be applied row-wise, hence this function.
#'
#' @param data Data to which the dropout mask should be applied
#' @param mask The dropout mask, a vector of 0 and 1.
#' @return Data with applied dropout mask
#'
#' @references Dropout: A Simple Way to Prevent Neural Networks from
#'  Overfitting, Nitish Srivastava
#' @seealso \url{https://www.cs.toronto.edu/~hinton/absps/JMLRdropout.pdf}


applyDropoutMask <- function(data, mask)
{
  return (data * matrix(rep(mask, nrow(data)), nrow=nrow(data), byrow=T))
}

