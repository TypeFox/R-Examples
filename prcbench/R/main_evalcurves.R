#' Evaluate Precision-Recall curves with specified tools and test sets
#'
#' The \code{run_evalcurve} function runs several tests to evaluate
#'    the accuracy of Precision-Recall curves.
#'
#' @param testset A character vector to specify a test set generated by
#'   \code{\link{create_testset}}.
#'
#' @param toolset A character vector to specify a tool set generated by
#'   \code{\link{create_toolset}}.
#'
#' @return A data frame with validation results.
#'
#' @seealso \code{\link{create_testset}} to generate a test dataset.
#'    \code{\link{create_toolset}} to generate a tool set.
#'
#' @examples
#' ## Evaluate curves for c1, c2, c3 test sets and crv5 tool set
#' testset <- create_testset("curve", c("c1", "c2", "c3"))
#' toolset <- create_toolset(set_names = "crv5")
#' res1 <- run_evalcurve(testset, toolset)
#'
#' @export
run_evalcurve <- function(testset, toolset) {
  # Validate arguments
  new_args <- .validate_run_evalcurve_args(testset, toolset)

  # Prepare tool sets and test data sets
  new_testset <- rep(new_args$testset, length(new_args$toolset))
  new_toolset <- rep(new_args$toolset, length(new_args$testset))

  # Evaluate curves
  testres <- .run_curve_tests(new_testset, new_toolset)
  summres <- .summarize_scores(testres, new_args$testset)
  bptsres <- .get_base_points(new_args$testset)
  predres <- .predict_curves(new_testset, new_toolset)

  reslst <- list(testscores = testres, testsum = summres, basepoints = bptsres,
                 predictions = predres)

  # Create an S3 object
  s3obj <- structure(reslst, class = "evalcurve")
}

#
# Validate a Precision-Recall curve
#
.run_curve_tests <- function(testset, toolset) {
  tfunc <- function(i) {
    tool <- toolset[[i]]
    tset <- testset[[i]]
    tool$call(tset)

    vres <- .eval_x_range(tool)
    vres <- rbind(vres, .eval_y_range(tool))
    vres <- rbind(vres, .eval_fpoint(tset, tool))
    vres <- rbind(vres, .eval_intpts(tset, tool))
    vres <- rbind(vres, .eval_epoint(tset, tool))
    rownames(vres) <- NULL

    resdf <- data.frame(testitem = c("x_range", "y_range", "fpoint", "intpts",
                                     "epoint"))
    scoredf <- cbind(resdf, vres)
    basedf <- data.frame(testset = tset$get_tsname(),
                         toolset = tool$get_setname(),
                         toolname = tool$get_toolname())
    cbind(basedf, scoredf)
  }
  res <- do.call(rbind, lapply(seq_along(testset), tfunc))
  sres <- res[order(res$testset, res$toolset, res$toolname), ]
  rownames(sres) <- NULL
  sres
}

#
# Check the x value range of a Precision-Recall curve
#
.eval_x_range <- function(tool) {

  # Test 1
  if (all(tool$get_x() >= 0, na.rm = T)) {
    test1 <- TRUE
  } else {
    test1 <- FALSE
  }

  # Test 2
  if (all(tool$get_x() <= 1, na.rm = T)) {
    test2 <- TRUE
  } else {
    test2 <- FALSE
  }

  if (test1 && test2) {
    success <- 1
  } else {
    success <- 0
  }

  scores <- c(success, 1)
  names(scores) <-  c("success", "total")
  scores
}

#
# Check the y value range of a Precision-Recall curve
#
.eval_y_range <- function(tool) {

  # Test 1
  if (all(tool$get_y() >= 0, na.rm = T)) {
    test1 <- TRUE
  } else {
    test1 <- FALSE
  }

  # Test 2
  if (all(tool$get_y() <= 1, na.rm = T)) {
    test2 <- TRUE
  } else {
    test2 <- FALSE
  }

  if (test1 && test2) {
    success <- 1
  } else {
    success <- 0
  }

  scores <- c(success, 1)
  names(scores) <-  c("success", "total")
  scores
}

#
# Check the first point of a Precision-Recall curve
#
.eval_fpoint <- function(tset, tool, tolerance = 1e-4) {
  bx <- tset$get_basepoints_x()
  by <- tset$get_basepoints_y()

  if (!is.na(tool$get_x()[1]) && !is.na(tool$get_y()[1])
        && (abs(tool$get_x()[1] - bx[1]) < tolerance)
        && (abs(tool$get_y()[1] - by[1]) < tolerance)) {
    success <- 1
  } else {
    success <- 0
  }

  scores <- c(success, 1)
  names(scores) <-  c("success", "total")
  scores
}

#
# Check intermediate points of a Precision-Recall curve
#
.eval_intpts <- function(tset, tool, tolerance = 1e-4) {
  bx <- tset$get_basepoints_x()
  by <- tset$get_basepoints_y()

  if (length(bx) > 2) {
    fpfunc <- function(i) {
      xidx <- abs(tool$get_x() - bx[i]) < tolerance
      ys <- tool$get_y()[xidx]

      if (any(abs(ys - by[i]) < tolerance, na.rm = T)) {
        return(1)
      } else {
        return(0)
      }
    }

    fcounts <- lapply(2:(length(bx) - 1), fpfunc)
    success <- do.call(sum, fcounts)
    total <- length(bx) - 2
    success <- min(success, total)
  } else {
    success <- 0
    total <- 0
  }

  scores <- c(success, total)
  names(scores) <-  c("success", "total")
  scores
}

#
# Check the end point of a Precision-Recall curve
#
.eval_epoint <- function(tset, tool, tolerance = 1e-4) {
  bx <- tset$get_basepoints_x()
  by <- tset$get_basepoints_y()

  epos1 <- length(tool$get_x())
  epos2 <- length(bx)

  if (!is.null(epos1) && !is.na(epos1) && (epos1 > 0)
      && !is.na(tool$get_x()[epos1])
      && !is.na(tool$get_y()[epos1])
      && (abs(tool$get_x()[epos1] - bx[epos2]) < tolerance)
      && (abs(tool$get_y()[epos1] - by[epos2]) < tolerance)) {
    success <- 1
  } else {
    success <- 0
  }

  scores <- c(success, 1)
  names(scores) <-  c("success", "total")
  scores
}

#
# Summarize curve evaluation results
#
.summarize_scores <- function(testres, testset) {
  sumdf <- stats::aggregate(testres[,c('success', 'total')],
                            by = list(testres$testset, testres$toolset,
                                      testres$toolname),
                            FUN = sum, na.rm = TRUE)
  colnames(sumdf)[1:3] <- c("testset", "toolset", "toolname")
  sumdf$label <- factor(paste0(sumdf$success, "/", sumdf$total))
  sumdf$lbl_pos_x <- 0
  sumdf$lbl_pos_y <- 0
  for (tset in testset) {
    tsname <- tset$get_tsname()
    sumdf[sumdf$testset == tsname, "lbl_pos_x"] <- tset$get_textpos_x()
    sumdf[sumdf$testset == tsname, "lbl_pos_y"] <- tset$get_textpos_y()
  }

  sres <- sumdf[order(sumdf$testset, sumdf$toolset, sumdf$toolname), ]
  rownames(sres) <- NULL
  sres
}

#
# Get base points from test datasets
#
.get_base_points <- function(testset) {
  bfunc <- function(tset) {
    tsname <- tset$get_tsname()
    bpx <- tset$get_basepoints_x()
    bpy <- tset$get_basepoints_y()
    data.frame(testset = rep(tsname, length(bpx)), x = bpx, y = bpy)
  }
  bpres <- do.call(rbind, lapply(testset, bfunc))
  rownames(bpres) <- NULL
  bpres
}

#
# Predict curves by tools with test datasets
#
.predict_curves <- function(testset, toolset) {
  pfunc <- function(i) {
    tool <- toolset[[i]]
    tset <- testset[[i]]
    tool$call(tset)

    tsname <- tset$get_tsname()
    setname <- tool$get_setname()
    toolname <- tool$get_toolname()
    x <- tool$get_x()
    y <- tool$get_y()

    data.frame(testset = rep(tsname, length(x)),
               toolset = rep(setname, length(x)),
               toolname = rep(toolname, length(x)), x = x, y = y)
  }
  predres <- do.call(rbind, lapply(seq_along(testset), pfunc))
  rownames(predres) <- NULL
  predres
}

#
# Validate arguments and return updated arguments
#
.validate_run_evalcurve_args <- function(testset, toolset) {

  assertthat::assert_that(is.list(testset))
  assertthat::assert_that(length(testset) > 0)
  for (tset in testset) {
    if (!methods::is(tset, "TestDataC")) {
      stop("Invalid testset", call. = FALSE)
    }
  }

  assertthat::assert_that(is.list(toolset))
  assertthat::assert_that(length(toolset) > 0)
  for (tool in toolset) {
    if (!methods::is(tool, "ToolIFBase")) {
      stop("Invalid toolset", call. = FALSE)
    }
    if(tool$get_setname() %in% c("auc5", "auc4")) {
      stop(paste0("Invalid predifend tool set: ", tool$get_setname()),
           call. = FALSE)
    }
    if(!environment(tool$clone)$private$def_store_res) {
      stop(paste0("Incorrect store_res value in ", tool$get_toolname()),
           call. = FALSE)
    }
  }

  list(testset = testset, toolset = toolset)
}