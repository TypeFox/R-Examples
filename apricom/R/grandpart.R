#' Accessory Function for Cross-validation
#'
#' Randomly partition a dataset into k approximately equally sized partitions.
#'
#' @param dataset     a dataset for partitioning.
#' @param k     the number of partitions; the cross-validation fold number.
#'
#' @return Returns a list containing the partitioned datasets.
#'
#' @note This function is not designed to be called directly,
#'       but acts within \code{kcrossval}

grandpart <- function(dataset, k) {

  n <- dim(dataset)[1]
  y <- sample(1:n, n)
  dc <- dataset[y, ]         ## permuted data
  q <- floor(n / k)
  r <- n %% k          ## remainder

  if (r > 0){
    d1 <- dc[1:(q + 1), ]
    Dats <- list(d1)
  }
  else {
    d1 <- dc[1:q, ]
    Dats <- list(d1)        ## creates list of k subsets
  }

  if (1 < r){
     for (l in 2:r) {
         Dats[[l]] <- dc[((l - 1) * (q + 1) + 1):(l * (q + 1)), ]
     }
  }

  for (l in (r + 1):k) {
    Dats[[l]] <- dc[((l - 1) * q + r + 1):(l * q + r), ]
  }
  return(Dats)
}

