#' Accessory Function for Sample Splitting
#'
#' Randomly partition a dataset into a training set and a test set.
#'
#' @param dataset     a dataset for splitting.
#' @param fract  the proportion of observations to be designated to the training set.
#' @return Returns a list containing the training set and test set.
#'
#' @note This function is not designed to be called directly,
#'       but acts within \code{splitval}

randpart <- function(dataset, fract){

  n <- dim(dataset)[1]
  y <- sample(1:n, n)
  ds <- dataset[y, ]
  nhalf <- floor(n * fract)
  d2 <- ds[1:nhalf, ]
  d1 <- ds[(nhalf + 1):n, ]
  dsplit <- list(d2, d1)
  return(dsplit)
}
