#
# ROCR
#
.rocr_wrapper <- function(testset, calc_auc = FALSE, store_res = TRUE) {
  if (!requireNamespace("ROCR", quietly = TRUE)) {
    stop("ROCR needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # Prepare data
  scores <- testset$get_scores()
  labels <- testset$get_labels()

  # Calculate Precision-Recall curve
  pred <- ROCR::prediction(scores, labels)
  perf <- ROCR::performance(pred, "prec", "rec")

  # Get AUC
  if (calc_auc) {
    x <- methods::slot(perf, "x.values")[[1]]
    y <- methods::slot(perf, "y.values")[[1]]

    # Copied the logic from .performance.auc of ROCR
    aucscore <- 0
    for (i in 2:length(x)) {
      aucscore <- aucscore + 0.5 * (x[i] - x[i-1]) * (y[i] + y[i-1])
    }

  } else {
    aucscore <- NA
  }

  # Return x and y values if requested
  if (store_res) {
    x <- methods::slot(perf, "x.values")[[1]]
    y <- methods::slot(perf, "y.values")[[1]]

    list(x = x, y = y, auc = aucscore)
  } else {
    NULL
  }
}

#
# PerfMeas
#
.pm_wrapper <- function(testset, calc_auc = FALSE, store_res = TRUE) {
  if (!requireNamespace("PerfMeas", quietly = TRUE)) {
    stop("PerfMeas needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # Prepare data
  scores <- testset$get_scores()
  labels <- testset$get_labels()

  # Calculate Precision-Recall curve
  prc <- PerfMeas::precision.at.all.recall.levels(scores, labels)

  # Get AUC
  if (calc_auc) {
    aucscore <- PerfMeas::AUPRC(list(prc), comp.precision = TRUE)
  } else {
    aucscore <- NA
  }

  # Return x and y values if requested
  if (store_res) {
    x <- prc[["recall"]]
    y <- prc[["precision"]]

    list(x = x, y = y, auc = aucscore)
  } else {
    NULL
  }
}

#
# PRROC
#
.prroc_wrapper <- function(testset, calc_auc = FALSE, store_res = TRUE,
                           curve = TRUE, minStepSize = 0.01) {
  if (!requireNamespace("PRROC", quietly = TRUE)) {
    stop("PRROC needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # Prepare data
  fg <- testset$get_fg()
  bg <- testset$get_bg()

  # Calculate Precision-Recall curve
  prc <- PRROC::pr.curve(fg, bg, curve = curve, minStepSize = minStepSize)

  # Get AUC
  if (calc_auc) {
    aucscore <- prc$auc.integral
  } else {
    aucscore <- NA
  }

  # Return x and y values if requested
  if (store_res) {
    x <- rev(prc[["curve"]][, 1])
    y <- rev(prc[["curve"]][, 2])

    list(x = x, y = y, auc = aucscore)
  } else {
    NULL
  }
}

#
# precrec
#
.precrec_wrapper <- function(testset, calc_auc = FALSE, store_res = TRUE) {
  if (!requireNamespace("precrec", quietly = TRUE)) {
    stop("precrec needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # Prepare data
  scores <- testset$get_scores()
  labels <- testset$get_labels()

  # Calculate Precision-Recall curve
  curves <- precrec::evalmod(scores = scores, labels = labels)

  # Get AUC
  if (calc_auc) {
    aucs <- precrec::auc(curves)
    aucscore <-  aucs[aucs$curvetypes == "PRC", ]
  } else {
    aucscore <- NA
  }

  # Return x and y values if requested
  if (store_res) {
    x <- curves[["prcs"]][[1]][["x"]]
    y <- curves[["prcs"]][[1]][["y"]]

    list(x = x, y = y, auc = aucscore)
  } else {
    NULL
  }
}
