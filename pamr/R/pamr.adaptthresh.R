pamr.adaptthresh <- function(object, ntries = 10, reduction.factor = 0.9, full.out = FALSE) {
  errors <- error.nsc(object)
  threshold <- object$threshold   
### Remove all but the first leading zero errors
  ifirst <- match(TRUE, object$errors > 0, FALSE)
  if (!ifirst)
    stop("Zero training error throughout!")
  else {
    ifirst <- max(ifirst, 1)
    threshold <- threshold[seq(ifirst, length(threshold))]
  }
### initialization
  tscales <- object$threshold.scale
  all.errors <- matrix(0, ntries + 1, length(tscales),
                       dimnames = list(NULL, names(tscales)))
  all.scales <- all.errors
  all.objects <- as.list(seq(ntries + 1))
  rocs <- double(ntries + 1)
  all.scales[1,  ] <- tscales
  all.errors[1,  ] <- errors
  rocs[1] <- roc.nsc(object)      # integrated size^(1/4)*error
  cat("Initial errors:", format(round(errors, 5)), "Roc",
      format(round(rocs[1], 5)), "\n")
  for (i in seq(ntries)) {
    cat("Update", i, "\n")
    j <- rev(order(errors))[1]      # identify the largest error
    tscales[j] <- tscales[j] * reduction.factor     
                                        # and reduce its scale
    all.scales[i + 1,  ] <- tscales/min(tscales)    # and renormalize
    iobject <- update(object, threshold = threshold, 
                      threshold.scale = all.scales[i + 1,  ], remove.zeros = 
                      FALSE)
    all.errors[i + 1,  ] <- errors <- error.nsc(iobject)
    rocs[i + 1] <- roc.nsc(iobject)
    cat("\nErrors", format(round(errors, 5)), "Roc",
        format(round(rocs[i + 1], 5)), "\n")
  }
  j <- order(rocs)[1]     # identify the scales with the smallest "roc"
  opt.scale <- all.scales[j,  ]
  if (full.out)
    list(errors = all.errors, scales = all.scales, rocs = rocs, 
         opt.scale = opt.scale)
  else
    opt.scale
}

pamr.batchadjust <- function(data) {
  if (is.null(data$batchlabels)) {
    stop("batch labels are not in data object")
  }
  lab <- data$batchlabels
  dd <- model.matrix( ~ factor(lab) - 1)
  data$x <- data$x - misreg.simple(dd, data$x)
  data
}


misreg.simple <- function(Y, x) {
###Y is a indicator response matrix
  nax <- is.na(x)
  nsamples <- (!nax)%*%Y
  x[nax] <- 0
  xsum <- x%*%Y
  xbar <- xsum/nsamples
  xbar %*% t(Y)
}
