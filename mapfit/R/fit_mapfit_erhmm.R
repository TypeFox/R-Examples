
mapfit.erhmm.options <- function() {
  list(maxiter = 2000,
       reltol = sqrt(.Machine$double.eps),
       abstol = +Inf)
}

mapfit.erhmm.verbose <- function() {
  list(emstep = FALSE,
    shape = TRUE,
    emprogress = 1)
}

mapfit.erhmm.structure.full <- function(size) {
  list(alpha=rep(1, size),
    P=matrix(1, size, size))
}

mapfit.erhmm.time <- function(phsize, difftime, erhmm,
  structure = mapfit.erhmm.structure.full,
  method = c("all", "increment"), lbound = 1, ubound = phsize,
  stationary = TRUE,
  control = list(), verbose = list(), ...) {
  data <- mapfit.time.data.frame(difftime)
  mapfit.erhmm(phsize, data, erhmm, structure, method, lbound, ubound, stationary, control, verbose)  
}

mapfit.erhmm <- function(phsize, data, erhmm,
  structure = mapfit.erhmm.structure.full,
  method = c("all", "increment"), lbound = 1, ubound = phsize,
  stationary = TRUE, 
  control = list(), verbose = list(), ...) {
  method <- match.arg(method)
  call <- match.call()

  con <- mapfit.erhmm.options()
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("unknown names in control: ", paste(noNms, collapse = ", "))

  ver <- mapfit.erhmm.verbose()
  nmsC <- names(ver)
  ver[(namc <- names(verbose))] <- verbose
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("unknown names in verbose: ", paste(noNms, collapse = ", "))

  ## init parameters
  if (method == "all") {
    tres <- system.time(result <- mapfit.erhmm.all(phsize, data, lbound, ubound,
      stationary, structure, control=con, verbose=ver, ...))
  }
  else if (method == "increment") {
    tres <- system.time(result <- mapfit.erhmm.increment(phsize, data, lbound, ubound,
      stationary, structure, control=con, verbose=ver, ...))
  }
  else {
    stop("Invalid method")
  }

  result <- c(result, list(data=data@data, ctime=tres[1], call=call))
  class(result) <- "mapfit.result"
  result
}

## algorithms for shape parameter vectors

mapfit.erhmm.all <- function(phsize, data, lbound, ubound,
  stationary, structure, control, verbose, ...) {
  maxllf <- -Inf
  allshape <- herlang.shape.all(phsize, lbound, ubound)
  for (s in allshape) {
    struct <- structure(length(s))
    m <- erhmm.param.kmeans(s, data@data, struct$alpha, struct$P, verbose, ...)
    result <- emfit(m, data, initialize=FALSE, stationary=stationary,
      control=control, verbose=verbose, ...)
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

mapfit.erhmm.increment <- function(phsize, data, lbound, ubound,
  stationary, structure, control, verbose, ...) {
  maxllf <- -Inf
  # maxllfv <- rep(-Inf, phsize)
  shape <- rep(1, lbound)
  repeat {
    shapelist1 <- herlang.shape.increment(shape, ubound, phsize)
    shapelist2 <- herlang.shape.decrement(shape, lbound)
    shapelist <- c(shapelist1, shapelist2)
    shape <- NULL
    for (s in shapelist) {
      struct <- structure(length(s))
      m <- erhmm.param.kmeans(s, data@data, struct$alpha, struct$P, verbose, ...)
      result <- emfit(m, data, initialize=FALSE, stationary=stationary,
        control=control, verbose=verbose, ...)
      if (verbose$shape)
        cat("shape: ", s, sprintf(" llf=%.2f\n", result$llf))
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

