#' An S4 class that defines the look-up table and all other components required for prediction using this table.
#'
#' @import dplyr
#' @importFrom  data.table data.table
#' @import methods
#' @slot table the look-up table with entries to be retrieved as prediction results
#' @slot feature.con a vector of continuous feature names
#' @slot feature.cat a vector of categorical feature names
#' @slot feature.boundaries a list of boundaries for each input feature (inferred during construction from input data)
#' @slot response the name of the response variable for the look-up table
#' @slot default the default value for cells corresponding to a missing combination of input values
#' @slot response.categories sequence of all categories (order-dependent) for the response variable, if it's categorical
#' @export lookupTable
#' @exportClass lookupTable
lookupTable <- setClass('lookupTable',
                        slots = c(feature.boundaries = 'list',
                                  feature.con = 'character',
                                  feature.cat = 'character',
                                  table = 'data.table',
                                  response = 'character',
                                  default = 'numeric',
                                  response.categories = 'character')
)

#' Initialize and construct a lookupTable object
#'
#' @param .Object the prototype object
#' @param df.input training data set containing columns with names found in features.con and features.cat vectors
#' @param response name of the response variable
#' @param feature.boundaries a list of thresholds for each continuous feature (names contained in feature.con) to construct bins. Should use -Inf and Inf as the first and last values, respectively.
#' @param features.con a vector of continuous feature names
#' @param features.cat a vector of categorical feature names
#' @param fill.method the method to fill entries of the table ('mean' or 'median')
#' @return A lookupTable object with a table trained with df.input data
#'
setMethod('initialize', 'lookupTable', function(.Object, df.input, response, feature.boundaries,
                                                features.con = character(0), features.cat = character(0),
                                                fill.method = 'mean') {
  df.input <- df.input[,c(response, features.con, features.cat)]

  if (class(df.input[[response]]) == 'factor')
    .Object@response.categories <- levels(df.input[[response]])
  else
    .Object@response.categories <- character(0)

  if (dim(df.input)[2] > 2) {
    if (length(features.cat) > 0) {
      for (i in 1:length(features.cat)) {
        df.input[, features.cat[i]] <- as.factor(df.input[, features.cat[i]])
        if (length(levels(df.input[[features.cat[i]]])) > (nrow(df.input)^(1/3)))
          stop('At least one categorical feature has too many levels, which should not be more than the cubic root of the sample size.')
      }
    }
  }
  if (length(features.con) > 0) {
    for (i in 1:length(features.con)) {
      df.input <- data.frame(df.input, bb = rep(0, dim(df.input)[1]))
      names.orig <- names(df.input)[1:(length(names(df.input)) - 1)]
      df.input <-  dplyr::mutate(df.input, bb = cut(df.input[[features.con[i]]], feature.boundaries[[i]]))
      names(df.input) <- c(names.orig, paste(features.con[i], 'cut', sep='_'))
    }
  }

  dots <- list()
  dt.key <- c()
  if (length(features.con) > 0) {
    dots <- lapply(paste(features.con, '_cut', sep=''), as.symbol)
    dt.key <- paste(features.con, '_cut', sep='')
  }
  features.cat.symbol <- lapply(features.cat, as.symbol)
  dots <- c(dots, features.cat.symbol)
  dt.key <- c(dt.key, features.cat)

  df.input <- dplyr::summarise_(dplyr::group_by_(df.input, .dots=dots),
                                paste(fill.method, '(', response, ')', sep=''))

  colnames(df.input)[length(df.input)] <- 'lookup'
  dt.output <- data.table::data.table(df.input, key=dt.key)

  .Object@feature.boundaries <- feature.boundaries
  .Object@feature.con <- features.con
  .Object@feature.cat <- features.cat
  .Object@table <- dt.output
  .Object@response <- response
  .Object@default <- mean(dt.output[['lookup']][!is.na(dt.output[['lookup']])])
  return(.Object)
}
)

#' \code{\link{predict}} method for \code{\linkS4class{lookupTable}} objects
#'
#' @title Predictions from a look-up table
#' @param object a fitted lookupTable object
#' @param newdata data.frame from which to evaluate predictions
#' @return a numeric vector of predicted values
#' @param newparams new parameters to use in evaluating predictions
#' @param ... optional additional parameters.  None are used at present.
#' @examples
#' df.input <- cars
#' response <- 'dist'
#' feature.boundaries <- list(c(-Inf, 5, 10, 15, 20, 25, Inf))
#' features.con <- c('speed')
#' dist.table <- lookupTable(df.input, response, feature.boundaries, features.con)
#' df.test <- data.frame(speed = c(2, 23, 41, 5, 9, 8))
#' predict(dist.table, df.test)
#' @method predict lookupTable
#' @export
predict.lookupTable <- function(object, newdata, newparams = NULL, ...) {
  features <- object@feature.con
  if (length(object@feature.cat) > 0) {
    for (i in 1:length(object@feature.cat)) {
      newdata[, object@feature.cat[i]] <- as.factor(newdata[, object@feature.cat[i]])
    }
  }
  bins <- list()
  if (length(features) > 0) {
    bins <- list(levels(object@table[[paste(features[1], '_cut', sep = '')]])[findInterval(newdata[[features[1]]], object@feature.boundaries[[1]])])
    if (length(features) > 1) {
      for (i in 2:length(features)) {
        bins <- c(bins, list(levels(object@table[[paste(features[i], '_cut', sep = '')]])[findInterval(newdata[[features[i]]], object@feature.boundaries[[i]])]))
      }
    }
  }
  if (length(object@feature.cat) > 0) {
    bins <- c(bins, list(newdata[, object@feature.cat[1]]) )
    if (length(object@feature.cat) > 1) {
      for (i in 2:length(object@feature.cat)) {
        bins <- c(bins,  list(newdata[, object@feature.cat[i]]))
      }
    }
  }
  arr.output <- object@table[bins][['lookup']]
  arr.output[is.na(arr.output)] <- object@default
  if (length(object@response.categories) > 0)
    arr.output <- object@response.categories[round(arr.output)]
  return(arr.output)
}


