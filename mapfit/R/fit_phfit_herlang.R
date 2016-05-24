
phfit.herlang.options <- function() {
  list(maxiter = 2000,
       reltol = sqrt(.Machine$double.eps),
       abstol = +Inf)
}

phfit.herlang.verbose <- function() {
  list(emstep = FALSE,
    shape = TRUE,
    emprogress = 1)
}

phfit.herlang <- function(phsize, data, herlang, method = c("all", "increment"),
  lbound = 1, ubound = phsize, control = list(), verbose = list(), ...) {
  method <- match.arg(method)
  call <- match.call()

  con <- phfit.herlang.options()
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("unknown names in control: ", paste(noNms, collapse = ", "))

  ver <- phfit.herlang.verbose()
  nmsC <- names(ver)
  ver[(namc <- names(verbose))] <- verbose
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("unknown names in verbose: ", paste(noNms, collapse = ", "))

  ## init parameters
  if (method == "all") {
    tres <- system.time(result <- phfit.herlang.all(phsize, data, lbound=lbound, ubound=ubound,
      control=con, verbose=ver, ...))
  }
  else if (method == "increment") {
    tres <- system.time(result <- phfit.herlang.increment(phsize, data, lbound=lbound, ubound=ubound,
      control=con, verbose=ver, ...))
  }
  else {
    stop("Invalid method")
  }

  result <- c(result, list(lbound=lbound, ubound=ubound,
    data=data@data, ctime=tres[1], call=call))
  class(result) <- "phfit.result"
  result
}

## algorithms for shape parameter vectors

phfit.herlang.all <- function(phsize, data, lbound, ubound, control, verbose, ...) {
  maxllf <- -Inf
  allshape <- herlang.shape.all(phsize, lbound, ubound)
  for (s in allshape) {
    result <- emfit(herlang(shape=s), data, initialize=TRUE, control=control, verbose=verbose, ...)
    if (verbose$shape)
      cat("shape: ", s, gettextf(" llf=%.2f\n", result$llf))
    if (is.finite(result$llf)) {
      if (result$llf > maxllf) {
        maxllf <- result$llf
        maxresult <- result
      }
    }
  }
  return(maxresult)
}

phfit.herlang.increment <- function(phsize, data, lbound, ubound, control, verbose, ...) {
  maxllf <- -Inf
  # maxllfv <- rep(-Inf, phsize)
  shape <- rep(1, lbound)
  repeat {
    shapelist1 <- herlang.shape.increment(shape, ubound, phsize)
    shapelist2 <- herlang.shape.decrement(shape, lbound)
    shapelist <- c(shapelist1, shapelist2)
    shape <- NULL
    for (s in shapelist) {
      result <- emfit(herlang(shape=s), data, initialize=TRUE, control=control, verbose=verbose, ...)
      if (verbose$shape)
        cat("shape: ", s, gettextf(" llf=%.2f\n", result$llf))
      if (is.finite(result$llf)) {
        # if (result$llf > maxllfv[k]) {
        #   maxllfv[k] <- result$llf
        #   # shape <- s
        # }
        if (result$llf > maxllf) {
          maxllf <- result$llf
          maxresult <- result
          shape <- s
        }
      }
    }
    if (is.null(shape)) {
##      warning(message="break")
      break
    }
  }
  return(maxresult)
}

