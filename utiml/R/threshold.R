# FIXED ------------------------------------------------------------------------
#' Apply a fixed threshold in the results
#'
#' Transfom a prediction matrix with scores/probabilities in a mlresult applying
#' a fixed threshold. A global fixed threshold can be used of all labels or
#' different fixed thresholds, one for each label.
#'
#' @family threshold
#' @param prediction A matrix with scores/probabilities where the columns
#'    are the labels and the rows are the instances.
#' @param threshold A single value between 0 and 1 or a list with threshold
#'    values contained one value per label.
#' @param probability A logical value. If \code{TRUE} the predicted values are
#'  the score between 0 and 1, otherwise the values are bipartition 0 or 1.
#'  (Default: \code{FALSE})
#' @return A mlresult object.
#' @references
#'  Al-Otaibi, R., Flach, P., & Kull, M. (2014). Multi-label Classification: A
#'  Comparative Study on Threshold Selection Methods. In First International
#'  Workshop on Learning over Multiple Contexts (LMCE) at ECML-PKDD 2014.
#' @export
#'
#' @examples
#' # Create a prediction matrix with scores
#' result <- matrix(
#'  data = rnorm(9, 0.5, 0.2),
#'  ncol = 3,
#'  dimnames = list(NULL, c('lbl1',  'lb2', 'lb3'))
#' )
#'
#' # Use 0.5 as threshold
#' fixed_threshold(result)
#'
#' # Use an threshold for each label
#' fixed_threshold(result, c(0.4, 0.6, 0.7))
fixed_threshold <- function (prediction, threshold = 0.5, probability = FALSE) {
  UseMethod("fixed_threshold")
}

#' @describeIn fixed_threshold Fixed Threshold for matrix or data.frame
#' @export
fixed_threshold.default <- function(prediction, threshold = 0.5,
                                    probability = FALSE) {
  if (length(threshold) == 1) {
    threshold <- rep(threshold, ncol(prediction))
  }
  else if (length(threshold) != ncol(prediction)) {
    stop(paste("The threshold values must be a single value or the same",
               "number of labels"))
  }

  bipartition <- do.call(cbind, lapply(seq(ncol(prediction)), function(col) {
    as.integer(prediction[, col] >= threshold[col])
  }))
  dimnames(bipartition) <- dimnames(prediction)

  multilabel_prediction(bipartition, prediction, probability)
}

#' @describeIn fixed_threshold Fixed Threshold for mlresult
#' @export
fixed_threshold.mlresult <- function (prediction, threshold = 0.5,
                                      probability = FALSE) {
  fixed_threshold.default(as.probability(prediction), threshold, probability)
}

# MCUT -------------------------------------------------------------------------
#' Maximum Cut Thresholding (MCut)
#'
#' The Maximum Cut (MCut) automatically determines a threshold for each instance
#' that selects a subset of labels with higher scores than others. This leads to
#' the selection of the middle of the interval defined by these two scores as
#' the threshold.
#'
#' @family threshold
#' @param prediction A matrix or mlresult.
#' @param probability A logical value. If \code{TRUE} the predicted values are
#'  the score between 0 and 1, otherwise the values are bipartition 0 or 1.
#'  (Default: \code{FALSE})
#' @return A mlresult object.
#' @references
#' Largeron, C., Moulin, C., & Gery, M. (2012). MCut: A Thresholding Strategy
#'  for Multi-label Classification. In 11th International Symposium, IDA 2012
#'  (pp. 172-183).
#' @export
#'
#' @examples
#' prediction <- matrix(runif(16), ncol = 4)
#' mcut_threshold(prediction)
mcut_threshold <- function (prediction, probability = FALSE) {
  UseMethod("mcut_threshold")
}

#' @describeIn mcut_threshold Maximum Cut Thresholding (MCut) method for matrix
#' @export
mcut_threshold.default <- function (prediction, probability = FALSE) {
  result <- apply(prediction, 1, function (row) {
    sorted.row <- sort(row, decreasing = TRUE)
    difs <- unlist(lapply(seq(length(row)-1), function (i) {
      sorted.row[i] - sorted.row[i+1]
    }))
    t <- which.max(difs)
    mcut <- (sorted.row[t] + sorted.row[t+1]) / 2
    row <- ifelse(row > mcut, 1, 0)
    row
  })

  multilabel_prediction(t(result), prediction, probability)
}

#' @describeIn mcut_threshold Maximum Cut Thresholding (MCut) for mlresult
#' @export
mcut_threshold.mlresult <- function (prediction, probability = FALSE) {
  mcut_threshold.default(as.probability(prediction), probability)
}

# PCUT -------------------------------------------------------------------------
#' Proportional Thresholding (PCut)
#'
#' Define the proportion of examples for each label will be positive.
#' The Proportion Cut (PCut) method can be a label-wise or global method that
#' calibrates the threshold(s) from the training data globally or per label.
#'
#' @family threshold
#' @param prediction A matrix or mlresult.
#' @param ratio A single value between 0 and 1 or a list with ratio values
#'  contained one value per label.
#' @param probability A logical value. If \code{TRUE} the predicted values are
#'  the score between 0 and 1, otherwise the values are bipartition 0 or 1.
#'  (Default: \code{FALSE})
#' @return A mlresult object.
#' @references
#' Al-Otaibi, R., Flach, P., & Kull, M. (2014). Multi-label Classification: A
#'  Comparative Study on Threshold Selection Methods. In First International
#'  Workshop on Learning over Multiple Contexts (LMCE) at ECML-PKDD 2014.
#'
#' Largeron, C., Moulin, C., & Gery, M. (2012). MCut: A Thresholding Strategy
#'  for Multi-label Classification. In 11th International Symposium, IDA 2012
#'  (pp. 172-183).
#' @export
#'
#' @examples
#' prediction <- matrix(runif(16), ncol = 4)
#' pcut_threshold(prediction, .45)
pcut_threshold <- function (prediction, ratio, probability = FALSE) {
  UseMethod("pcut_threshold")
}

#' @describeIn pcut_threshold Proportional Thresholding (PCut) method for matrix
#' @export
pcut_threshold.default <- function (prediction, ratio, probability = FALSE) {
  n <- nrow(prediction)
  num.elem <- ceiling(ratio * n)
  if (length(num.elem) == 1) {
    num.elem <- rep(num.elem, ncol(prediction))
    names(num.elem) <- colnames(prediction)
  }
  else if (length(num.elem) != ncol(prediction)) {
    stop(paste("The number of elements values must be a single value or the",
               "same number of labels"))
  }
  else if (is.null(names(num.elem))) {
    names(num.elem) <- colnames(prediction)
  }

  indexes <- utiml_rename(seq(ncol(prediction)), colnames(prediction))
  result <- do.call(cbind, lapply(indexes, function (ncol) {
    values <- c(rep(1, num.elem[ncol]), rep(0, n - num.elem[ncol]))
    prediction[order(prediction[, ncol], decreasing=TRUE), ncol] <- values
    prediction[, ncol]
  }))

  multilabel_prediction(result, prediction, probability)
}

#' @describeIn pcut_threshold Proportional Thresholding (PCut) for mlresult
#' @export
pcut_threshold.mlresult <- function (prediction, ratio, probability = FALSE) {
  pcut_threshold.default(as.probability(prediction), ratio, probability)
}

# RCUT -------------------------------------------------------------------------
#' Rank Cut (RCut) threshold method
#'
#' The Rank Cut (RCut) method is an instance-wise strategy, which outputs the k
#' labels with the highest scores for each instance at the deployment.
#'
#' @family threshold
#' @param prediction A matrix or mlresult.
#' @param k The number of elements that will be positive.
#' @param probability A logical value. If \code{TRUE} the predicted values are
#'  the score between 0 and 1, otherwise the values are bipartition 0 or 1.
#'  (Default: \code{FALSE})
#' @return A mlresult object.
#' @references
#'  Al-Otaibi, R., Flach, P., & Kull, M. (2014). Multi-label Classification: A
#'  Comparative Study on Threshold Selection Methods. In First International
#'  Workshop on Learning over Multiple Contexts (LMCE) at ECML-PKDD 2014.
#' @export
#'
#' @examples
#' prediction <- matrix(runif(16), ncol = 4)
#' rcut_threshold(prediction, 2)
rcut_threshold <- function (prediction, k, probability = FALSE) {
  UseMethod("rcut_threshold")
}

#' @describeIn rcut_threshold Rank Cut (RCut) threshold method for matrix
#' @export
rcut_threshold.default <- function (prediction, k, probability = FALSE) {
  values <- c(rep(1, k), rep(0, ncol(prediction) - k))
  result <- apply(prediction, 1, function (row) {
    row[order(row, decreasing = TRUE)] <- values
    row
  })
  multilabel_prediction(t(result), prediction, probability)
}

#' @describeIn rcut_threshold Rank Cut (RCut) threshold method for mlresult
#' @export
rcut_threshold.mlresult <- function (prediction, k, probability = FALSE) {
  rcut_threshold.default(as.probability(prediction), k, probability)
}

# SCORE DRIVEN -----------------------------------------------------------------
score_driven_threshold <- function () {
  #TODO
}
# #' Cost-based loss function for multi-label classification
# #'
# #' @param mdata A mldr dataset containing the test data.
# #' @param mlresult An object of mlresult that contain the scores and bipartition
# #'  values.
# #' @param cost The cost of classification each positive label. If a single value
# #'  is informed then the all labels have tha same cost.
# #' @references
# #'  Al-Otaibi, R., Flach, P., & Kull, M. (2014). Multi-label Classification: A
# #'  Comparative Study on Threshold Selection Methods. In First International
# #'  Workshop on Learning over Multiple Contexts (LMCE) at ECML-PKDD 2014.
# multilabel_loss_function <- function (mdata, mlresult, cost = 0.5) {
#   if (length(cost) == 1) {
#     cost <- rep(cost, mdata$measures$num.labels)
#     names(cost) <- rownames(mdata$labels)
#   }
#   else if (is.null(names(cost))) {
#     names(cost) <- rownames(mdata$label)
#   }
#
#   prediction <- as.bipartition(mlresult)
#   labels <- utiml_rename(rownames(mdata$labels))
#   partial.results <- lapply(labels, function (lname) {
#     FN <- sum(mdata$dataset[,lname] == 1 & prediction [,lname] == 0) /
#       mdata$measures$num.instances
#     FP <- sum(mdata$dataset[,lname] == 0 & prediction [,lname] == 1) /
#       mdata$measures$num.instances
#     freq <- mdata$labels[lname, "freq"]
#     2 * ((cost[lname] * freq * FN) + ((1 - cost[lname]) * (1 - freq) * FP))
#   })
#
#   mean(unlist(partial.results))
# }

# SCUT -------------------------------------------------------------------------
#' SCut Score-based method
#'
#' This is a label-wise method that adjusts the threshold for each label to
#' achieve a specific loss function using a validation set or cross validation.
#'
#' Different from the others threshold methods instead of return the bipartition
#' results it returs the threshold values for each label.
#'
#' @family threshold
#' @param prediction A matrix or mlresult.
#' @param expected The expected labels for the prediction. May be a matrix with
#'  the label values or a mldr object.
#' @param loss.function A loss function to be optmized. If you want to use your
#'  own error function see the notes and example. (Default: Mean Squared Error)
#' @param cores The number of cores to parallelize the computation Values higher
#'  than 1 require the \pkg{parallel} package. (Default:
#'  \code{options("utiml.cores", 1)})
#' @return A numeric vector with the threshold values for each label
#' @note The loss function is a R method that receive two vectors, the expected
#'  values of the label and the predicted values, respectively. Positive values
#'  are represented by the 1 and the negative by the 0.
#' @references
#'  Fan, R.-E., & Lin, C.-J. (2007). A study on threshold selection for
#'   multi-label classification. Department of Computer Science, National
#'   Taiwan University.
#'
#'  Al-Otaibi, R., Flach, P., & Kull, M. (2014). Multi-label Classification: A
#'   Comparative Study on Threshold Selection Methods. In First International
#'   Workshop on Learning over Multiple Contexts (LMCE) at ECML-PKDD 2014.
#' @export
#'
#' @examples
#' names <- list(1:10, c("a", "b", "c"))
#' prediction <- matrix(runif(30), ncol = 3, dimnames = names)
#' classes <- matrix(sample(0:1, 30, rep = TRUE), ncol = 3, dimnames = names)
#' thresholds <- scut_threshold(prediction, classes)
#' fixed_threshold(prediction, thresholds)
#'
#' \dontrun{
#' # Penalizes only FP predictions
#' mylossfunc <- function (real, predicted) {
#'    mean(predicted - real * predicted)
#' }
#' prediction <- predict(br(toyml, "RANDOM"), toyml)
#' scut_threshold(prediction, toyml, loss.function = mylossfunc, cores = 5)
#' }
scut_threshold <- function (prediction, expected, loss.function = NA,
                            cores = getOption("utiml.cores", 1)) {
  UseMethod("scut_threshold")
}

#' @describeIn scut_threshold Default scut_threshold
#' @export
scut_threshold.default <- function (prediction, expected, loss.function = NA,
                                    cores = getOption("utiml.cores", 1)) {
  if (cores < 1) {
    stop("Cores must be a positive value")
  }

  if (!is.function(loss.function)) {
    # Mean Squared Error
    loss.function <- function(real, predicted) {
      mean((real - predicted) ^ 2)
    }
  }

  if (class(expected) == "mldr") {
    expected <- expected$dataset[expected$labels$index]
  }

  labels <- utiml_rename(colnames(prediction))
  thresholds <- utiml_lapply(labels, function (col) {
    scores <- prediction[, col]
    index <- order(scores)
    ones <- which(expected[index, col] == 1)
    difs <- c(Inf)
    for (i in seq(length(ones)-1)) {
      difs <- c(difs, ones[i+1] - ones[i])
    }

    evaluated.thresholds <- c()
    result <- c()
    for (i in ones[which(difs > 1)]) {
      thr <- scores[index[i]]
      res <- loss.function(expected[, col], ifelse(scores < thr, 0, 1))
      evaluated.thresholds <- c(evaluated.thresholds, thr)
      result <- c(result, res)
    }

    ifelse(length(ones) > 0,
           as.numeric(evaluated.thresholds[which.min(result)]),
           max(scores) + 0.0001) # All expected values are in the negative class
  }, cores)

  unlist(thresholds)
}

#' @describeIn scut_threshold Mlresult scut_threshold
#' @export
scut_threshold.mlresult <- function (prediction, expected, loss.function = NA,
                                     cores = getOption("utiml.cores", 1)) {
  scut_threshold.default(as.probability(prediction), expected,
                         loss.function, cores)
}

# SUBSET CORRECTION ------------------------------------------------------------
#' Subset Correction of a predicted result
#'
#' This method restrict a multi-label learner to predict only label combinations
#' whose existence is present in the (training) data. To this all labelsets
#' that are predicted but are not found on training data is replaced by the most
#' similar labelset.
#'
#' If the most simillar is not unique, those label combinations with higher
#' frequency in the training data are preferred. The Hamming loss distance is
#' used to determine the difference between the labelsets.
#'
#' @family threshold
#' @param mlresult An object of mlresult that contain the scores and bipartition
#'  values.
#' @param train_y A matrix/data.frame with all labels values of the training
#'  dataset or a mldr train dataset.
#' @param base.threshold A numeric value between 0 and 1 to use as base to
#'  determine which values needs be reescaled to preserve the corrected
#'  labelsets. If \code{NULL} the score correction is ignored. (Default: NULL)
#' @param probability A logical value. If \code{TRUE} the predicted values are
#'  the score between 0 and 1, otherwise the values are bipartition 0 or 1.
#'  (Default: \code{FALSE})
#' @return A new mlresult where all results are present in the training
#'  labelsets.
#' @note The original paper describes a method to create only bipartitions
#'  result, but we adapeted the method to change the scores. Based on the
#'  base.threshold value the scores higher than the threshold value, but must be
#'  lower are changed to respect this restriction. If \code{NULL} this
#'  correction will be ignored.
#' @references
#'  Senge, R., Coz, J. J. del, & Hullermeier, E. (2013). Rectifying classifier
#'    chains for multi-label classification. In Workshop of Lernen, Wissen &
#'    Adaptivitat (LWA 2013) (pp. 162-169). Bamberg, Germany.
#' @export
#'
#' @examples
#' prediction <- predict(br(toyml, "RANDOM"), toyml)
#' subset_correction(prediction, toyml)
subset_correction <- function(mlresult, train_y, base.threshold = NULL,
                              probability = FALSE) {
  bip <- as.bipartition(mlresult)
  prob <- as.probability(mlresult)

  if (class(train_y) == "mldr") {
    train_y <- train_y$dataset[train_y$labels$index]
  }

  if (ncol(mlresult) != ncol(train_y)) {
    stop("The number of columns in the predicted result are different from the
         training data")
  }

  # Bipartition correction
  labelsets <- as.matrix(unique(train_y))
  rownames(labelsets) <- apply(labelsets, 1, paste, collapse = "")

  order <- names(sort(table(apply(train_y, 1, paste, collapse = "")),
                      decreasing = TRUE))
  labelsets <- labelsets[order, ]

  new.pred <- t(apply(bip, 1, function(y) {
    labelsets[names(which.min(apply(labelsets, 1, function(row) {
      sum(row != y)
    }))), ]
  }))

  # Probabilities correction
  new.prob <- prob
  if (!is.null(base.threshold)) {
    for (r in seq(nrow(prob))) {
      row <- prob[r, ]

      max_index <- new.pred[r, ] - row > base.threshold
      min_index <- new.pred[r, ] - row <= -base.threshold

      indexes <- min_index | max_index
      max_v <- min(c(row[row > base.threshold & !indexes],
                     base.threshold + 0.1))
      min_v <- max(c(row[row < base.threshold & !indexes],
                     base.threshold - 0.1))

      # Normalize values
      new.prob[r, max_index] =
        row[max_index] * (max_v - base.threshold) + base.threshold
      new.prob[r, min_index] =
        row[min_index] * (base.threshold - min_v) + min_v
    }
  }

  multilabel_prediction(new.pred, new.prob, probability)
}
