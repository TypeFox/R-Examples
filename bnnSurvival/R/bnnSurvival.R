


##' Bootstrap aggregated (bagged) version of the k-nearest neighbors survival probability prediction method (Lowsky et al. 2013). 
##' In addition to the bootstrapping of training samples, the features can be subsampled in each base learner. 
##'
##' For a description of the k-nearest neighbors survival probability prediction method see (Lowsky et al. 2013). 
##' Please note, that parallel processing, as currently implemented, does not work on Microsoft Windows platforms.
##' 
##' The weighting function needs to be defined for all distances >= 0. 
##' The default function is constant 1, a possible alternative is w(x) = 1/(1+x). 
##' 
##' To use the non-bagged version as in Lowsky et al. 2013, use \code{num_base_learners=1}, \code{replace=FALSE} and \code{sample_fraction=1}.
##' 
##' @title Bagged k-nearest neighbors survival prediction
##' @param formula Object of class formula or character describing the model to fit.
##' @param data Training data of class data.frame.
##' @param k Number nearest neighbors to use. If a vector is given, the optimal k of these values is found using 5-fold cross validation.
##' @param num_base_learners Number of base learners to use for bootstrapping.
##' @param num_features_per_base_learner Number of features randomly selected in each base learner. Default: all.
##' @param metric Metric d(x,y) used to measure the distance between observations. Currently only "mahalanobis".
##' @param weighting_function Weighting function w(d(,x,y)) used to weight the observations based on their distance.
##' @param replace Sample with or without replacement. 
##' @param sample_fraction Fraction of observations to sample in [0,1]. Default is 1 for \code{replace=TRUE}, and 0.6321 for \code{replace=FALSE}. 
##' @return bnnSurvivalEnsemble object. Use predict() with a new data set to predict survival probabilites.
##' @examples
##' require(bnnSurvival)
##' 
##' ## Use only 1 core
##' options(mc.cores = 1)
##' 
##' ## Load a dataset and split in training and test data
##' require(survival)
##' n <- nrow(veteran)
##' idx <- sample(n, 2/3*n)
##' train_data <- veteran[idx, ]
##' test_data <- veteran[-idx, ]
##' 
##' ## Create model with training data and predict for test data
##' model <- bnnSurvival(Surv(time, status) ~ trt + karno + diagtime + age + prior, train_data, 
##'                      k = 20, num_base_learners = 10, num_features_per_base_learner = 3)
##' result <- predict(model, test_data)
##' 
##' ## Plot survival curve for the first observations
##' plot(timepoints(result), predictions(result)[1, ])
##' 
##' @references
##'   Lowsky, D.J. et al. (2013). A K-nearest neighbors survival probability prediction method. Stat Med, 32(12), 2062-2069.
##' @seealso \code{\link[=predict,bnnSurvivalEnsemble-method]{predict}}
##' @author Marvin N. Wright
##' @import stats 
##' @importFrom Rcpp evalCpp
##' @useDynLib bnnSurvival
##' @export
bnnSurvival <- function(formula, data, k = max(1, nrow(data)/10), num_base_learners = 50,
                        num_features_per_base_learner = NULL, metric = "mahalanobis",
                        weighting_function = function(x){x*0+1}, 
                        replace = TRUE, sample_fraction = NULL) {

  ## Generate model and matrix for training data
  formula <- formula(formula)
  if (class(formula) != "formula") {
    stop("Error: Invalid formula.")
  }
  train_model <- model.frame(formula, data)
  train_matrix <- data.matrix(cbind(train_model[, 1][, c(1,2)], train_model[, -1]))

  ## Check arguments
  if (is.numeric(num_base_learners) & !is.na(num_base_learners) & num_base_learners > 0) {
    num_base_learners <- as.integer(num_base_learners)
  } else {
    stop("num_base_learners is no positive number.")
  }
  if (is.null(num_features_per_base_learner)) {
    num_features_per_base_learner <- as.integer(ncol(train_matrix) - 2)
  } else if (is.numeric(num_features_per_base_learner) & !is.na(num_features_per_base_learner)) {
    if (num_features_per_base_learner <= 0) {
      stop("num_features_per_base_learner is no positive number.")
    } else if (num_features_per_base_learner > ncol(train_matrix) - 2) {
      stop("num_features_per_base_learner is larger than number of features.")
    } else {
      num_features_per_base_learner <- as.integer(num_features_per_base_learner)
    }   
  } else {
    stop("num_features_per_base_learner is no number.")
  }
  if (is.null(sample_fraction)) {
    if (replace) {
      sample_fraction <- 1
    } else {
      sample_fraction <- 0.6321
    }
  } 
  if (sample_fraction > 1 | sample_fraction < 0) {
    stop("sample_fraction is not in [0,1] interval.")
  }
  if (is.numeric(k) && length(k) == 1 && !is.na(k) && k > 0) {
    k <- as.integer(k)
  } else if (is.numeric(k) && length(k) > 1 && all(!is.na(k)) && all(k > 0)) {
    k <- get_best_k(formula = formula, data = data, k = k, 
                    num_base_learners = num_base_learners,
                    num_features_per_base_learner = num_features_per_base_learner, 
                    metric = metric, weighting_function = weighting_function, 
                    replace = replace, sample_fraction = sample_fraction)
    k <- as.integer(k)
  } else {
    stop("k is no positive number.")
  }
  if (k > nrow(train_matrix)) {
    stop("k cannot be larger than the number of observations in data.")
  }

  ## Create ensemble of base learners
  ensemble <- bnnSurvivalEnsemble(train_data = train_matrix,
                  formula = formula,
                  num_base_learners = num_base_learners,
                  num_features_per_base_learner = num_features_per_base_learner,
                  k = k, metric = metric, weighting_function = weighting_function, 
                  replace = replace, sample_fraction = sample_fraction)

  return(ensemble)
}

