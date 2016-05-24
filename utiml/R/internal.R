#' Conditional value selection
#'
#' @param test an object which can be coerced to logical mode.
#' @param yes object that will be returned when the test value is true.
#' @param no object that will be returned when the test value is false
#' @return The respective value yes or no based on test value. This is an
#' alternative way to use a single logical value for avoid the real if/else for
#' choice lists, matrices and other composed data.
#'
#' @examples
#' \dontrun{
#' utiml_ifelse(TRUE, dataframe1, dataframe2) ## dataframe1
#' utiml_ifelse(length(my.list) > 10, my.list[1:10], my.list)
#' }
utiml_ifelse <- function(test, yes, no) {
  list(yes, no)[c(test, !test)][[1]]
}

#' Select the suitable method lapply or mclaplly
#'
#' @param mylist a list to iterate.
#' @param myfnc The function to be applied to each element of the mylist.
#' @param utiml.cores The number of cores to use. If 1 use lapply oterwise use
#'    mclapply.
#' @param utiml.seed A numeric value to set a seed to execute in parallel mode.
#' @param ... Extra arguments to myfnc.
#' @return A list with the results of the specified method.
utiml_lapply <- function(mylist, myfnc, utiml.cores, utiml.seed = NA, ...) {
  mylist <- as.list(mylist)
  indexes <- seq_along(mylist)
  names(indexes) <- names(mylist)

  if (anyNA(utiml.seed)) {
    thefunc <- function (i, ...) {
      myfnc(mylist[[i]], ...)
    }
  } else {
    thefunc <- function (i, ...) {
      set.seed(utiml_ifelse(is.null(utiml.seed),
                            NULL, as.numeric(utiml.seed) + i))
      myfnc(mylist[[i]], ...)
    }
  }

  if (requireNamespace("parallel", quietly = TRUE)) {
    parallel::mclapply(indexes,
                       thefunc,
                       mc.cores = min(utiml.cores, length(mylist)),
                       ...)
  }
  else {
    lapply(indexes, thefunc, ...)
  }
}

#' Internal normalize data function
#'
#' @param data a set of numbers.
#' @param max.val The maximum value to normalize. If NULL use the max value
#'   present in the data. (default: \code{NULL})
#' @param min.val The minimum value to normalize. If NULL use the min value
#'   present in the data (default: \code{NULL})
#' @return The normalized data
#'
#' @examples
#' \dontrun{
#' utiml_normalize(c(1,2,3,4,5))
#' #--> 0 0.25 0.5 0.75 1
#'
#' utiml_normalize(c(1,2,3,4,5), 10, 0)
#' #--> 0.1 0.2 0.3 0.4 0.5
#' }
utiml_normalize <- function(data, max.val = NULL, min.val = NULL) {
  max.val <- ifelse(is.null(max.val), max(data, na.rm = TRUE), max.val)
  min.val <- ifelse(is.null(min.val), min(data, na.rm = TRUE), min.val)
  (data - min.val)/(max.val - min.val)
}

#' Return the newdata to a data.frame or matrix
#'
#' @param newdata The data.frame or mldr data
#' @return A dataframe or matrix containing only dataset
#'
#' @examples
#' \dontrun{
#' test <- emotions$dataset[,emotions$attributesIndexes]
#' all(test == utiml_newdata(emotions)) # TRUE
#' all(test == utiml_newdata(test)) # TRUE
#' }
utiml_newdata <- function(newdata) {
  UseMethod("utiml_newdata")
}

#' @describeIn utiml_newdata Return the data in the original format
utiml_newdata.default <- function(newdata) {
  newdata
}

#' @describeIn utiml_newdata Return the dataset from the mldr dataset
utiml_newdata.mldr <- function(newdata) {
  newdata$dataset[, newdata$attributesIndexes]
}

#' Preserve current seed
utiml_preserve_seed <- function () {
  if (exists('.Random.seed', envir = .GlobalEnv, inherits = FALSE)) {
    current.seed <- get('.Random.seed', envir = .GlobalEnv, inherits = FALSE)
    scope <- parent.frame()
    scope$utiml.current.seed <- current.seed
  }
}

#' Rename the list using the names values or its own content
#'
#' @param X A list
#' @param names The list names, If empty the content of X is used
#' @return A list with the new names
#' @export
#'
#' @examples
#' utiml_rename(c("a", "b", "c"))
#' ## c(a="a", b="b", c="c")
#'
#' utiml_rename(c(1, 2, 3), c("a", "b", "c"))
#' ## c(a=1, b=2, c=3)
utiml_rename <- function (X, names = NULL) {
  names(X) <- utiml_ifelse(is.null(names), X, names)
  X
}

#' Restore the current seed
utiml_restore_seed <- function () {
  scope <- parent.frame()
  if (!is.null(scope$utiml.current.seed)) {
    assign('.Random.seed', scope$utiml.current.seed, envir = .GlobalEnv)
  }
}

#' Define if two sets are equals independently of the order of the elements
#'
#' @param a A list
#' @param b Other list
#' @return Logical value where TRUE the sets are equals and FALSE otherwise.
#' @examples
#' \dontrun{
#' utiml_is_equal_sets(c(1, 2, 3), c(3, 2, 1))
#' ## TRUE
#'
#' utiml_is_equal_sets(c(1, 2, 3), c(1, 2, 3, 4))
#' ## FALSE
#' }
utiml_is_equal_sets <- function (a, b) {
  length(setdiff(union(a, b), intersect(a, b))) == 0
}
