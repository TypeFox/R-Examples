#' @title Predict FLCs given new PLCs
#'
#' @description 
#' This function predicts FLCs given new PLCs based on the estimated 
#' \eqn{\epsilon} mappings and estimated conditional distributions.
#' 
#' @param train a list of training examples with
#' LC observations (a list of \code{PLC} and \code{FLC}), \code{weight.matrix}, and \code{pdfs}
#' @param test a list of test examples with PLC observations and/or the
#' \code{weight.matrix} associated with the PLC observations.
#' @param method estimation method for estimating PLC and FLC distributions
#' @param type prediction: \code{'mean'}, \code{'median'}, \code{'weightedmean'},
#' or \code{'mode'}.
#' @return
#' \eqn{N \times K} matrix
#' @keywords methods
#' @export
#'

predict_FLC_given_PLC <- function(train = list(data = list(FLC = NULL, PLC = NULL), 
                                               weight.matrix = NULL,
                                               pdfs = list(FLC = NULL, PLC = NULL)),
                                  test = list(PLC = NULL, weight.matrix = NULL),
                                  type = c("weighted.mean", "mean",
                                           "median", "mode"),
                                  method = list(FLC = "nonparametric",
                                                PLC = "normal")) {
  # TODO: make it work for multivariate FLC input
  type <- match.arg(type)  
  if (is.null(train$weight.matrix)) {
    stop("The trained model weight matrix must be provided.")
  }
  num.states <- ncol(train$weight.matrix)
  if (is.null(test$weight.matrix)) {
    pdf.of.PLC.test <- estimate_LC_pdfs(LCs = train$data$PLC, 
                                        weight.matrix = train$weight.matrix,
                                        method = method[["PLC"]],
                                        eval.LCs = test$PLC)
    test$weight.matrix <- 
      estimate_state_probs(weight.matrix = train$weight.matrix,
                           states = NULL, 
                           pdfs = list(PLC = pdf.of.PLC.test, 
                                       FLC = NULL),
                           num.states = num.states)
  }

  if (is.null(dim(train$data$FLC))) {
    train$data$FLC <- matrix(train$data$FLC, ncol = 1)
  }
  if (type == "weighted.mean") {
    # sum up observed FLCs (weighted sum)
    pred.per.state <- t(train$weight.matrix) %*% train$data$FLC
    # divide by sample size of each column
    pred.per.state <- sweep(pred.per.state, 1, colSums(train$weight.matrix), "/")
  } else if (type %in% c("mean", "median")) {
    state.vector <- weight_matrix2states(train$weight.matrix)
    if (type == "mean") {
      pred.per.state <- sapply(seq_len(num.states),
                               function(x) {
                                 colMeans(as.matrix(train$data$FLC[state.vector == x, ]))
                               })
    } else if (type == "median") {
      pred.per.state <- sapply(seq_len(max(state.vector)),
                               function(x) {
                                 apply(as.matrix(train$data$FLC[state.vector == x, ]),
                                       2, median)
                               })
    }
  } else if (type == "mode") {
    if (is.null(train$pdfs$FLC)) {
      train$pdfs$FLC <- estimate_LC_pdfs(train$data$FLC, train$weight.matrix)
    }
    pred.per.state <- train$data$FLC[apply(train$pdfs$FLC, 2, which.max), ]
  }
  
  if (!is.null(dim(pred.per.state))) {
    if (nrow(pred.per.state) != ncol(test$weight.matrix)) {
      pred.per.state <- t(pred.per.state) #t(pred.per.state)
    } 
  }
  # make pred_tests for test data by weighting with weight matrix
  pred.test <- as.matrix(test$weight.matrix %*% pred.per.state)
  pred.test <- sweep(pred.test, 1, rowSums(test$weight.matrix), "/")
  return(pred.test)
}
