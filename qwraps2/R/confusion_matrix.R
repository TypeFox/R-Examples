#' @title Confusion Matrices (Contingency Tables)
#'
#' @description Construction of confusion matrices, accuracy, sensitivity,
#' specificity, confidence intervals (Wilson's method and (optional
#' bootstrapping)).
#'
#' @param x prediction condition vector, a two level factor variable or a
#' variable that can be converted to one.
#' @param ... not currently used
#'
#' @details
#' Sensitivity and Specificity:
#' For the sensitivity and specificity function we expect the 2-by-2 confusion
#' matrix (contingency table) to be of the form:
#'
#' \tabular{lccc}{
#'                     \tab      \tab True \tab Condition \cr
#'                     \tab      \tab +    \tab -         \cr
#' Predicted Condition \tab +    \tab TP   \tab FP        \cr
#' Predicted Condition \tab -    \tab FN   \tab TN        \cr
#' }
#' where
#' \itemize{
#'   \item FN: False Negative, and
#'   \item FP: False Positive, 
#'   \item TN: True Negative,
#'   \item TP: True Positive.
#' }
#'
#' Recall: 
#' \itemize{
#'   \item sensitivity = TP / (TP + FN)
#'   \item specificity = TN / (TN + FP) 
#'   \item positive predictive value (PPV) = TP / (TP + FP)
#'   \item negative predictive value (NPV) = TN / (TN + FN)
#' }
#'
#' @return The sensitivity and specificity functions return numeric values.
#' \code{confusion_matrix} returns a list with elements:
#' \itemize{
#'   \item tab the confusion matrix,
#'   \item stats a matrix of summary statistics and confidence intervals.
#' }
#'
#' @examples
#' ## Example taken from caret::confusionMatrix
#' \donttest{ 
#' lvs <- c("normal", "abnormal")
#' truth <- factor(rep(lvs, times = c(86, 258)),
#'                 levels = rev(lvs))
#' pred <- factor(c(rep(lvs, times = c(54, 32)),
#'                  rep(lvs, times = c(27, 231))),               
#'                levels = rev(lvs))
#' 
#' confusion_matrix(pred, truth)
#' confusion_matrix(pred, truth)$stats
#' 
#' confusion_matrix(pred, truth, boot = TRUE)
#' confusion_matrix(pred, truth, positive = "normal", boot = TRUE)
#'
#' # Using formulas
#' test_data <- data.frame(xyz = pred, yyy = truth)
#' confusion_matrix(yyy ~ xyz, test_data)
#' confusion_matrix(yyy ~ xyz, test_data, positive = "normal") 
#' }
#'
#' @export
#' @rdname confusion_matrix
confusion_matrix <- function(x, ...) { 
  UseMethod("confusion_matrix")
}

#' @param y True Condition vector with the same possible values as x.
#' @param positive the level of x and y which is the positive outcome.  If
#' missing the first level of factor(y) will be used as the positive level.
#' @param boot boolean, should bootstrapped confidence intervals for the
#' sensitivity and specificity be computed?  Defaults to FALSE.
#' @param boot_samples number of bootstrapping sample to generate, defaults to
#' 1000L.  Ignored if \code{boot == FALSE}.
#' @param alpha 100(1-alpha)% confidence intervals for specificity and
#' sensitivity.  Ignored if \code{boot == FALSE}.
#' @export
#' @rdname confusion_matrix
confusion_matrix.default <- function(x, y, positive, boot = FALSE, boot_samples = 1000L, alpha = 0.05, ...) { 
  confusion_matrix.formula(stats::as.formula(paste(deparse(substitute(y)), deparse(substitute(x)), sep = "~")),
                           data = stats::setNames(data.frame(x,  y), c(deparse(substitute(x)), deparse(substitute(y)))),
                           positive,
                           boot, boot_samples, alpha)
}

#' @param formula column (known) ~ row (test) for building the confusion matrix
#' @param data environment containing the variables listed in the formula
#' @export
#' @rdname confusion_matrix
confusion_matrix.formula <- function(formula, data = parent.frame(), positive, boot = FALSE, boot_samples = 1000L, alpha = 0.05, ...) { 

  .data <- stats::model.frame(formula, data)
  .data[[1]] <- factor(.data[[1]])
  .data[[2]] <- factor(.data[[2]]) 

  if (!missing(positive)) {
    # Add error handing here
    .data[[1]] <- stats::relevel(.data[[1]], positive)
    .data[[2]] <- stats::relevel(.data[[2]], positive) 
  }

  if (nlevels(.data[[1]]) != nlevels(.data[[2]]) | nlevels(.data[[1]]) != 2) { 
    stop("qwraps2::confusion_matrix only supports factors with two levels.")
  }

  if (!all(levels(.data[[1]]) %in% levels(.data[[2]]))) { 
    stop("qwraps2::confusion_matrix expectes the same levels for the factors.")
  } 

  tab <- table(.data[[2]], .data[[1]], dnn = c("Prediction", "Truth"))#rev(names(.data)))

  stats <- rbind(Accuracy = accuracy(tab), 
                 Sensitivity = sensitivity(tab),
                 Specificity = specificity(tab), 
                 PPV = ppv(tab), 
                 NPV = npv(tab))

  stats <- cbind(stats, t(apply(stats, 1, wilson_score_interval, n = nrow(.data), alpha = alpha))) 
  colnames(stats) <- c("Est", "LCL", "UCL")

  if (boot) { 
    rows <- replicate(boot_samples, 
                      sample(seq(1, nrow(data), by = 1), nrow(data), replace = TRUE), 
                      simplify = FALSE)
    boot_stats <- 
      lapply(rows, 
             function(x) { 
               tab <- table(data[x, 1], data[x, 2])

               rbind(Accuracy = accuracy(tab), 
                     Sensitivity = sensitivity(tab),
                     Specificity = specificity(tab),
                     PPV = ppv(tab), 
                     NPV = npv(tab))
             }) 
    boot_stats <- do.call(cbind, boot_stats)

    boot_stats <- apply(boot_stats, 1, 
                         function(x) {
                           c(mean(x), stats::quantile(x, probs = c(alpha / 2, 1 - alpha / 2)))
                         })
    boot_stats <- t(boot_stats)
    colnames(boot_stats) <- c("Boot Est", "Boot LCL", "Boot UCL")

    stats <- cbind(stats, boot_stats)
  }

  rtn <- list(tab = tab, stats = stats)

  class(rtn) <- c("confusion_matrix", class(rtn))
  attr(rtn, "boot") <- boot
  attr(rtn, "alpha") <- alpha
  attr(rtn, "var_names") <- stats::setNames(as.list(names(.data)), c("Truth", "Prediction"))

  rtn 
}

#' @rdname confusion_matrix
#' @export
is.confusion_matrix <- function(x) inherits(x, "confusion_matrix")

#' @rdname confusion_matrix
#' @export
print.confusion_matrix <- function(x, ...) { 
  cat("\nTruth:      ", attr(x, "var_names")[[1]], 
      "\nPrediction: ", attr(x, "var_names")[[2]], "\n\n")
  print.table(x$tab) 
  print(x$stats) 
  invisible(x) 
}

accuracy <- function(tab) { 
  if (any(dim(tab) != 2)) { stop("Incorrect dim(tab)") } 
  as.numeric(sum(diag(tab)) / sum(tab))
}

ppv <- function(tab) { 
  as.numeric(tab[1, 1] / sum(tab[1, ]))
}

npv <- function(tab) { 
  as.numeric(tab[2, 2] / sum(tab[2, ]))
}

sensitivity <- function(tab) { 
  if (any(dim(tab) != 2)) { stop("Incorrect dim(tab)") } 
  as.numeric(tab[1, 1] / sum(tab[, 1]))
}

specificity <- function(tab, ...) { 
  if (length(dim(tab)) != 2 | any(dim(tab) != 2)) { stop("Incorrect dim(tab)") } 
  as.numeric(tab[2, 2] / sum(tab[, 2]))
}

wilson_score_interval <- function(p, n, alpha = 0.05) { 
  z <- stats::qnorm(1 - alpha/2) 
  1 / (1 + 1/n * z^2) * (p + 1 / (2 * n) * z^2 + c(-z, z) * sqrt( 1 / n * p * (1 - p) + 1 / (4 * n^2) * z^2)) 
}
