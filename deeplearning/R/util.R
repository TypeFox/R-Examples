#' Prints out the weight of a deep neural network
#'
#' This function prints out the weight in a heat map, 3D surface, or histogram
#'
#' @param darch DArch instance
#' @param num_of_layer the number of the layer to print
#' @param show_derivative T to show the weight value. F to show the percentage
#' weight change in the finetuning stage. This helps spot the network saturation problem.
#' @param type type of the graph. It supports "heatmap", "surface", and "histogram"
#'
#' @importFrom darch getLayer
#'
#' @examples
#' # Example of Regression
#'
#' input <- matrix(runif(1000), 500, 2)
#' input_valid <- matrix(runif(100), 50, 2)
#' target <- rowSums(input + input^2)
#' target_valid <- rowSums(input_valid + input_valid^2)
#' # create a new deep neural network for classificaiton
#' dnn_regression <- new_dnn(
#'  c(2, 50, 50, 20, 1),  # The layer structure of the deep neural network.
#'  # The first element is the number of input variables.
#'  # The last element is the number of output variables.
#'  hidden_layer_default = rectified_linear_unit_function,
#'  # for hidden layers, use rectified_linear_unit_function
#'  output_layer_default = linearUnitDerivative
#'  # for regression, use linearUnitDerivative function
#')
#'
#' # print the layer weights
#' # this function can print heatmap, histogram, or a surface
#' print_weight(dnn_regression, 1, type = "heatmap")
#'
#' print_weight(dnn_regression, 2, type = "surface")
#'
#' print_weight(dnn_regression, 3, type = "histogram")
#'
#'
#' @export

print_weight <- function(darch, num_of_layer, show_derivative = F, type = "heatmap") {
  weight <- getLayer(darch, num_of_layer)[[1]]
  weight_change <- getLayer(darch, num_of_layer)[[3]] / weight[1:(dim(weight)[[1]]-1), ]

  if(type == "histogram") {
    if(!show_derivative) {
      plotly::plot_ly(x = c(weight), type = type)
    } else {
      plotly::plot_ly(x = c(weight_change), type = type)
    }
  } else {
    if(!show_derivative) {
      plotly::plot_ly(z = weight, type = type, colorscale = "hot")
    } else {
      plotly::plot_ly(z = weight_change, type = type, colorscale = "hot")
    }
  }
}



#' Calculates the outer product of two matricies
#'
#' Calcualtes the outer product of two matrices
#'
#' @param data the date matrix
#' @param weight the weight matrix
#'

matMult <- function(data, weight) {
  return(data %*% weight)
}




#' Data proprosess function that covnerts a categorical input to continuous input or
#' vectorize it
#'
#' Proprosess a data set. It converts categorical data into binary variables
#' if it is unordered or continuous variable from 0 to 1 if it is ordinal
#' @param x input variable
#' @param type ordinal or other
#' @param ordered_list the rank ordering of an ordinal variable. Users are expected to
#' provide a complete list of the rank ordering. Otherwise, a default rank ordering
#' will be used.
#' @param var_name the name of the input variable. This is used to to create vectorized
#' input variables
#' @param ... other inputs
#'
#' @export

convert_categorical <- function(x,
                                type = "ordinal",
                                ordered_list = list(),
                                var_name = "var",
                                ...) {

  if(type == "ordinal") {
    unique_x <- unique(x)

    if(is.null(ordered_list)) {
      ordered_list <- sort(unique_x) # list_x has all unique values in vector x
    }

    if(any(!(unique_x %in% ordered_list))) {
      ordered_list <- sort(unique_x) # list_x has all unique values in vector x
    }

    num_categories <- length(ordered_list)
    mapped_value <- c(0:(num_categories - 1))/(num_categories - 1)
    ret <- mapped_value[match(x, ordered_list)]
  } else {
    unique_x <- unique(x)
    ordered_list <- sort(unique_x)
    num_categories <- length(ordered_list)
    mapped_value <- c(1:num_categories)
    numeric_x <- mapped_value[match(x, ordered_list)]
    vectorized_x <- matrix(0, nrow = length(x), ncol = length(unique_x))
    for( i in 1:length(x)) {
      vectorized_x[i, numeric_x[i]] <- 1
    }

    ret <- data.frame(vectorized_x)
    colnames(ret) <- paste0(var_name, " = ", ordered_list)
  }

  return(ret)
}

#' Creates a matrix by repeating a row vector N times
#'
#' helper function that repeat a row vector N times
#'
#' @param vector the row vector
#' @param N number of rows in the output matirx
#' @return a matrix

verticalize <- function(vector, N) {
  return(matrix(rep(vector, N), N, byrow = T))
}


