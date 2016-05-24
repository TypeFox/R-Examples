#' Data shaping tools
#'
#' Generate a new dataset in a format that is accepted by \code{\link{compare}}. Dummy variables are
#' created for categorical variables.
#'
#' This function can be used to prepare a dataset before applying the
#' \code{\link{compare}} function. The outcome column number must be
#' specified, and specific predictors and observation subsets may be specified. 2-level
#' categorical variables will be converted to binary, while dummy variables will be
#' created for categorical predictors with greater than two levels.
#'
#' The "datashaped" dataset should be saved to a new object.
#'
#' @import utils
#'
#' @param dataset   a dataframe or matrix containing user data.
#' @param y         the column number of the outcome variable in dataset d.
#' @param x         a vector containing the explanatory or predictor variables of interest.
#'                  The new data matrix will only contain these variables and the outcome
#'                  variable.
#' @param subs  a vector describing a subset of rows from dataset d to include in the
#'              returned data matrix.
#' @return This function returns a matrix conforming to the specifications supplied to
#' the datashape function.
#' @section Warning:
#'  \code{datashape} will not function if missing values are present.
#'
#'  @examples
#'## Preparing the iris dataset
#'    data(iris)
#'    iris.shaped <- datashape(dataset = iris, y = 4)
#'    head(iris.shaped)
#'## Creating a copy of iris with sepal-related predictors and a subset of observations.
#'    iris.sub <- datashape(dataset = iris, y = 4, x = c(1,2), subs = c(1:20, 50:70))
#'    head(iris.sub)

datashape<- function(dataset, y, x, subs) {

  if (missing(dataset)) stop("Dataset required for shaping")
  if (missing(y)) stop("Outcome variable must be specified")

  nc <- dim(dataset)[2]
  if (missing(x)) {
    if (y == 1) { d2 <- dataset[, c(2:nc, 1)]}
    else if (y ==  nc) {d2 <- dataset}
    else {d2 <- dataset[, c(1:(y - 1), (y + 1):nc, y)]}
  }
  else d2 <- dataset[, c(x, y)]

  zero_one <- function(dataset){

    dataset <- as.data.frame(dataset)
    if (!is.numeric(dataset)) {
      nncols <- sapply(dataset, is.numeric)
      nncols <- as.numeric(nncols)
      nncols <- which(nncols == 0)
    }
    for (k in nncols){dataset[, k] <- zero_one_vec(dataset[, k])}
    if (is.numeric(dataset)){dataset <- as.matrix(dataset)}
    return(dataset)
  }

  zero_one_vec <- function(v){

    if (length(unique(v)) < 3){
      zero <- sort(unique(v))[1]
      y <- as.numeric(v == zero)
      return(y)
    } else {return(v)}
  }

  dummify <- function(dataset){

    dataset <- as.data.frame(dataset)
    n <- dim(dataset)[2]
    D1 <- list("1")
    for (k in 1:n) {D1[[k]] <- dummify1(dataset[, k])}
    y <- do.call("cbind", D1)
    if (is.numeric(y)){y <- as.matrix(y)}
    colnames(y) <- varnames(dataset)
    return(y)
  }

  dummify1 <- function(v, g){

    if (missing(g)) g = 1
    if (!is.numeric(v)) {
      m <- length(v)
      vlevels <- sort(unique(v))
      nl <- length(vlevels)
      y <- rep(0, m * nl)
      y <- matrix(y, ncol = nl)
      for (k in 1:nl) y[, k] <- as.numeric(v == vlevels[k])
      y <- y[, -g]
    } else {y <- v}
    return(y)
  }

  varnames <- function(dataset){
    y <- colnames(dataset)
    y <- as.list(y)
    n <- length(y)
    for (k in 1:n){
      if (!is.numeric(dataset[, k])) {
        y[[k]] <- as.character(tail(sort(unique(dataset[, k])),
                                    length(unique(dataset[,  k])) - 1))
        y[[k]] <- paste0("ind_", y[[k]])
      }
    }
    z <- y[[1]]
    for (k in 2:n) {z <- c(z, y[[k]])}
    return(z)
  }

  d2 <- d2[subs,]
  d2 <- zero_one(d2)
  d2 <- dummify(d2)
  d2 <- as.matrix(d2)

  return(d2)
}
