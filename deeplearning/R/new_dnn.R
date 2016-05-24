#' Creats a new instance of darch class
#'
#' This function creates a new instance of darch class
#'
#' @param layer_structure a int vector that specifies the number and width of layers
#' @param layer_functions a list of activation functions used by each layer
#' @param output_layer_default the activation function for the output layer
#' @param hidden_layer_default the activation function for the hidden layers
#' @param weight_initiliazaiton function that initialize a layer's weight matrix
#'
#' @importFrom darch linearUnitDerivative generateWeights  createDataSet
#' @importFrom methods new
#' @importClassesFrom darch DArch
#' @examples
#' # create a new deep neural network for classificaiton
#' dnn_regression <- new_dnn(
#'  c(2, 50, 50, 20, 1),
#'  # The layer structure of the deep neural network.
#'  # The first element is the number of input variables.
#'  # The last element is the number of output variables.
#'  hidden_layer_default = rectified_linear_unit_function,
#'  # for hidden layers, use rectified_linear_unit_function
#'  output_layer_default = sigmoidUnitDerivative
#'  # for classification, use sigmoidUnitDerivative function
#' )
#'
#' # create a new deep neural network for classificaiton
#'dnn_regression <- new_dnn(
#'  c(2, 50, 50, 20, 1),
#'  # The layer structure of the deep neural network.
#'  # The first element is the number of input variables.
#'  # The last element is the number of output variables.
#'  hidden_layer_default = rectified_linear_unit_function,
#'  # for hidden layers, use rectified_linear_unit_function
#'  output_layer_default = linearUnitDerivative
#'  # for regression, use linearUnitDerivative function
#')
#' @export

new_dnn <- function(layer_structure,
                    layer_functions = NULL,
                    output_layer_default = linearUnitDerivative,
                    hidden_layer_default = rectified_linear_unit_function,
                    weight_initiliazaiton = generateWeights) {
  if (!is.null(layer_structure)) {
    # new a darch instance
    darch <-new("DArch")

    # set up the darch stats veriable
    darch@stats <-
      list("dataErrors" = list("raw"=c(), "class" = c()),
           "validErrors" = list("raw"=c(), "class" = c()),
           "times" = c(), "preTrainTime" = 0, "fineTuneTime" = 0)

    # set up the layers
    numLayers <- length(layer_structure)
    for (i in 1:(numLayers -1)) # first layer is an input layer
    {
      layer <- list()
      # element 1: initialize the layer weights
      dim_1 <- layer_structure[[i]]
      dim_2 <- layer_structure[[i + 1]]
      layer[[1]]  <- weight_initiliazaiton(dim_1 + 1, dim_2)

      # element 2: set up the layer activation function
      if (is.null(layer_functions[[as.character(i)]]))
      {
        if (i < (numLayers - 1)) {
          layer[[2]] <- hidden_layer_default
        } else {
          layer[[2]] <- output_layer_default
        }
      }
      else
      {
        layer[[2]] <- layer_functions[[as.character(i)]]
      }

      # element 3: weight increase
      layer[[3]] <- matrix(0, dim_1, dim_2)

      # element 4: gamma coefficient in batch normalization
      layer[[4]] <- rep(1, dim_2)

      # element 5: mu coefficient in batch normalization
      layer[[5]] <- rep(0, dim_2)

      # element 6: sigma_2 coefficient in batch normalization
      layer[[6]] <- rep(1 - exp(-12), dim_2)

      # add layer to darch@layers
      darch@layers[[i]] <- layer
    }

    # set up the slots necessary for predict.DArch function
    darch@dataSet <- createDataSet(matrix(0, 1, layer_structure[[1]]), NULL)
    darch@ff <- F

    # set up the execution function
    darch@executeFunction <- run_dnn
  } else {
    darch <- NULL
    flog.fatal("Illegal layer structures!")
  }
  return (darch)
}
