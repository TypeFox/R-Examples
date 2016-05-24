
##' Get optimal number of neighbors for bnnSurvival by cross validation
##' 
##' @title Get optimal number of neighbors
##' @param formula Formula
##' @param data Data 
##' @param k Number of neighbors
##' @param ... Further arguments passed to bnnSurvival 
##' @return Optimal k
##' @importFrom prodlim Hist
##' @importFrom pec pec
get_best_k <- function(formula, data, k, ...) {

  ## Use 5-fold cross validation
  K_cross <- 5

  ## Split data
  n <- nrow(data)
  folds <- split(sample(n, n), cut(1:n, K_cross, labels = FALSE))

  ## Compute integrated Brier score for each fold
  ibs <- sapply(folds, function(fold) {
    dat_test <- data[fold, ]
    dat_train <- data[-fold, ]

    ## Models
    models <- lapply(k, bnnSurvival, formula = formula, data = dat_train, ...)

    ## Compute integrated Brier score
    fitpec <- pec::pec(models, formula, dat_test, times = sort(unique(dat_train$time)),
                  cens.model = 'marginal', reference = FALSE)
    return(pec::crps(fitpec)[,1])
  })

  ## Return best k
  k[which.min(rowMeans(ibs))]
}
