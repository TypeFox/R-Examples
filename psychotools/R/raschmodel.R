## workhorse fitting function
raschmodel <- function(y, weights = NULL, start = NULL, reltol = 1e-10, 
  deriv = c("sum", "diff", "numeric"), hessian = TRUE, maxit = 100L,
  full = TRUE, gradtol = reltol, iterlim = maxit, ...)
{
  ## argument matching
  if(missing(reltol) && !missing(gradtol) && !is.null(gradtol)) reltol <- gradtol
  if(missing(maxit) && !missing(iterlim) && !is.null(iterlim)) maxit <- iterlim
  deriv <- match.arg(deriv)

  ## original data
  y <- as.matrix(y)
  k <- k_orig <- ncol(y)
  n <- nrow(y)
  if(is.null(colnames(y))) colnames(y) <- paste("Item", gsub(" ", "0", format(1:k)), sep = "")

  ## weights processing
  if(is.null(weights)) weights <- rep.int(1L, n)
  ## data and weights need to match
  stopifnot(length(weights) == n)

  ## omit zero weights
  weights_orig <- weights
  y_orig <- y
  y <- y[weights > 0, , drop = FALSE]
  weights <- weights[weights > 0]
  n <- nrow(y)

  ## all parameters identified?
  if(n < 2) stop("not enough observations")
  cm <- colMeans(y, na.rm = TRUE)
  status <- as.character(cut(cm, c(-Inf, 1/(2 * n), 1 - 1/(2 * n), Inf), labels = c("0", "0/1", "1")))
  status[is.na(status)] <- "NA"
  status <- factor(status, levels = c("0/1", "0", "1", "NA"))
  ident <- status == "0/1"
  names(status) <- colnames(y)

  ## just estimate identified parameters
  y_orig <- y_orig[,ident, drop = FALSE]
  y <- y[,ident, drop = FALSE]
  k <- ncol(y)
  y_na <- is.na(y)
  any_y_na <- any(y_na)

  if(!any_y_na) {
    ## compute likelihood/gradient/hessian on aggregated data
  
    ## data statistics
    cs <- colSums(y * weights)
    rs <- rowSums(y)
    rf <- as.vector(tapply(weights, factor(rs, levels = 0:k), sum))
    rf[is.na(rf)] <- 0

    ## starting values
    ## contrast: set parameter 1 to zero
    if(is.null(start)) {
      start <- log(sum(weights) - cs) - log(cs) #previously:# -qlogis(cs/sum(weights))
      start <- start[-1] - start[1]
    }
    rf <- rf[-1]
    cs <- cs[-1]

    ## objective function: conditional log-likelihood
    cloglik <- function(par) {
      ## obtain esf and apply contrast
      esf <- elementary_symmetric_functions(c(0, par), order = 0, diff = deriv == "diff")
      g <- esf[[1]][-1]

      ## conditional log-likelihood
      cll <- sum(-cs * par) - sum(rf * log(g))

      ## catch degenerated cases (typically cause by non-finite gamma)
      if(is.na(cll) | !is.finite(cll)) cll <- -.Machine$double.xmax

      return(-cll)
    }

    ## analytical gradient
    agrad <- function(par) {

        ## calculate esf
        esf <- elementary_symmetric_functions(c(0, par), order = 1, diff = deriv == "diff")

        ## calculate gradient
        - colSums(weights * (- y + esf[[2]][rs + 1, , drop = FALSE] / esf[[1]][rs + 1])[,-1, drop = FALSE])
    }

    ## analytical hessian
    ahessian <- function(par, esf) {
      ## obtain esf and apply contrast
      g <- esf[[1]][-1]
      g1 <- esf[[2]][-1, -1, drop = FALSE]
      g2 <- esf[[3]][-1, -1, -1, drop = FALSE]

      ## hessian
      hess <- matrix(0, ncol = k-1, nrow = k-1)
      g1s <- g1/g
      for (q in 1:(k-1)) hess[q,] <- colSums(rf * (g2[,q,]/g - (g1[,q]/g) * g1s))
    
      return(hess)
    }

  } else {
    ## compute likelihood/gradient/hessian on individual data

    ## process NA patterns and calculate static things once
    na_patterns <- factor(apply(y_na, 1, function(z) paste(which(z), collapse = "\r")))
    na_i <- wi_i <- wi2_i <- cs_i <- rs_i <- rf_i <- k_i <- vector("list", nlevels(na_patterns))
    na_pattern_levels <- levels(na_patterns)

    for (i in seq_along(na_pattern_levels)) {

       ## parse NA pattern
       na_level_i <- na_pattern_levels[i]
       wi_i[[i]] <- as.integer(strsplit(na_level_i, "\r")[[1]])
       wi2_i[[i]]<- if (length(wi_i[[i]]) < 1) 1:k else (1:k)[-wi_i[[i]]]
       k_i[[i]] <- length(wi2_i[[i]])

       ## select subset
       na_i[[i]] <- which(na_patterns == na_level_i)
       if(length(wi_i[[i]]) < 1) y_i <- y[na_i[[i]], , drop = FALSE]
       else y_i <- y[na_i[[i]], -wi_i[[i]], drop = FALSE]
       weights_i <- weights[na_i[[i]]]
       cs_i[[i]] <- colSums(y_i * weights_i)
       rs_i[[i]] <- rowSums(y_i)
       rf_i[[i]] <- as.vector(tapply(weights_i, factor(rs_i[[i]], levels = 0:ncol(y_i)), sum))
       rf_i[[i]][is.na(rf_i[[i]])] <- 0
    }

    ## starting values
    if(is.null(start)) {
      cs <- colSums(y * weights, na.rm = TRUE)
      ws <- colSums(!y_na * weights)
      start <- log(ws - cs) - log(cs) #previously:# -qlogis(cs/ws)
      start <- start[-1] - start[1]
    }

    ## convenience function
    zero_fill <- function(obj, at) {
      if(length(at) < 1) return(obj)
      if(is.null(dim(obj))) {
        rval <- rep.int(0, length(obj) + length(at))
	rval[-at] <- obj
      } else {
        rval <- matrix(0, ncol = ncol(obj) + length(at), nrow = nrow(obj) + length(at))
	rval[-at,-at] <- obj      
      }
      return(rval)
    }

    ## conditional log-likelihood function for NA pattern i
    cll_i <- function (cs_i, par_i, rf_i) {
        sum(-cs_i * par_i) - sum(rf_i * log(elementary_symmetric_functions(par_i, order = 0, diff = deriv == "diff")[[1]]))
    }

    ## objective function: conditional log-likelihood
    cloglik <- function(par) {

      ## initialize return values and extract esf parameters
      cll <- 0
      par_i <- lapply(wi_i, function(x) if (length(x) < 1) c(0, par) else c(0, par)[-x])
      
      ## conditional log-likelihood
      cll <- sum(mapply(cll_i, cs_i, par_i, rf_i, SIMPLIFY = TRUE, USE.NAMES = FALSE))

      ## catch degenerated cases (typically cause by non-finite gamma)
      if(is.na(cll) | !is.finite(cll)) cll <- -.Machine$double.xmax

      ## collect and return
      return(-cll)
    }
 
    ## analytical gradient
    agrad <- function(par) {

      ## initialize return value and esf parameters
      rval <- matrix(0, nrow = n, ncol = k)
      par_i <- lapply(wi_i, function(x) if (length(x) < 1) c(0, par) else c(0, par)[-x])
      esf_i <- mapply(elementary_symmetric_functions, par = par_i, MoreArgs = list(order = 1, diff = deriv == "diff"), SIMPLIFY = FALSE)

      ## loop over observed NA patterns	 
      for(i in seq_along(levels(na_patterns))) {
          rval[na_i[[i]], wi2_i[[i]]] <- weights[na_i[[i]]] * (- y[na_i[[i]], wi2_i[[i]], drop = FALSE] +
          esf_i[[i]][[2]][rs_i[[i]] + 1, , drop = FALSE] / esf_i[[i]][[1]][rs_i[[i]] + 1])
      }
    
      return(- colSums(rval[, -1, drop = FALSE]))
    }

    ## analytical hessian
    ahessian <- function(par, esf) {

      ## set up return value    
      rval <- matrix(0, ncol = k-1, nrow = k-1)

      ## loop over observed NA patterns      
      for(i in seq_along(levels(na_patterns))) {
	
        ## obtain esf
        g_i <- esf[[i]][[1]]
        g1_i <- esf[[i]][[2]]
        g2_i <- esf[[i]][[3]]

        ## hessian
	hess <- matrix(0, nrow = k_i[[i]], ncol = k_i[[i]])
        for (q in 1:k_i[[i]]) hess[q,] <- colSums(rf_i[[i]] * (g2_i[,q,]/g_i - (g1_i[,q]/g_i) * g1_i/g_i))
        rval <- rval + zero_fill(hess, wi_i[[i]])[-1, -1]
      }

      return(rval)
    }

  }
  
  ## optimization
  if(maxit > 0L) {
  opt <- optim(par = start, fn = cloglik, gr = agrad, method = "BFGS",
               hessian = (deriv == "numeric") & hessian, control = list(reltol = reltol, maxit = maxit, ...))
  } else {
    opt <- list(
      estimate = start,
      minimum = cloglik(start),
      hessian = if(deriv != "numeric") ahessian(start, esf) else NULL, ## no numeric Hessian available here
      iterations = 0,
      code = 4)
  }
  
  ## collect and annotate results
  cf <- opt$par
  names(cf) <- colnames(y)[-1]
  if(full) {
    esf <- if(any_y_na) {
      lapply(levels(na_patterns), function(z) {
        wi <- as.integer(strsplit(z, "\r")[[1]])
        cfi <- if(length(wi) < 1) c(0, cf) else c(0, cf)[-wi]
        elementary_symmetric_functions(cfi,
          order = 2 - (deriv == "numeric" | !hessian), 
	  diff = deriv == "diff")
      })
    } else {
      elementary_symmetric_functions(c(0, cf),
        order = 2 - (deriv == "numeric" | !hessian),
        diff = deriv == "diff")
    }
    if(any_y_na) names(esf) <- levels(na_patterns)
  
    if(hessian) {
      vc <- if(deriv == "numeric") opt$hessian else ahessian(cf, esf)
      vc <- solve(vc)
    } else {
      vc <- matrix(NA, nrow = length(cf), ncol = length(cf))
    }
    rownames(vc) <- colnames(vc) <- names(cf)
  } else {
    esf <- NULL
    vc <- NULL
  }

  ## collect, class, and return
  rval <- list(
    coefficients = cf,
    vcov = vc,
    loglik = -opt$value,
    df = k-1,
    data = y_orig,
    weights = if(identical(as.vector(weights_orig), rep(1L, nrow(y_orig)))) NULL else weights_orig,
    n = sum(weights_orig > 0),
    items = status,
    na = any_y_na,
    elementary_symmetric_functions = esf,
    code = opt$convergence,
    iterations = tail(na.omit(opt$counts), 1L),
    reltol = reltol,
    deriv = deriv        
  )
  class(rval) <- "raschmodel"
  return(rval)
}

## methods
coef.raschmodel <- function(object, ...) object$coefficients

vcov.raschmodel <- function(object, ...) object$vcov

logLik.raschmodel <- function(object, ...) structure(object$loglik, df = object$df, class = "logLik")

weights.raschmodel <- function(object, ...) if(!is.null(object$weights)) object$weights else rep(1, nrow(object$data))

print.raschmodel <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Rasch model difficulty parameters:\n")
  print(coef(x), digits = digits)
  invisible(x)
}

worth.raschmodel <- function(object, ...) {
    
    ## check if difficulty argument is used, if yes, return warning.
    addargs <- list(...)
    if ("difficulty" %in% names(addargs)) warning("The argument 'difficulty' is deprecated and not longer used.")

    ## call itempar
    return(itempar(object, ...))
}

summary.raschmodel <- function(object, vcov. = NULL, ...)
{
  ## coefficients
  cf <- coef(object)

  ## covariance matrix
  if(is.null(vcov.)) 
      vc <- vcov(object)
  else {
      if(is.function(vcov.)) vc <- vcov.(object)
        else vc <- vcov.
  }
  
  ## Wald test of each coefficient
  cf <- cbind(cf, sqrt(diag(vc)), cf/sqrt(diag(vc)), 2 * pnorm(-abs(cf/sqrt(diag(vc)))))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

  object$coefficients <- cf      
  class(object) <- "summary.raschmodel"
  return(object)
}

print.summary.raschmodel <- function(x, digits = max(3, getOption("digits") - 3), 
    signif.stars = getOption("show.signif.stars"), ...)
{
  if(is.null(x$call)) {
    cat("\nRasch model\n\n")  
  } else {
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  }

  if(any(x$items != "0/1")) cat("Excluded items:",
    paste(names(x$items)[x$items != "0/1"], collapse = ", "), "\n\n")

  cat("Difficulty parameters:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)

  cat("\nLog-likelihood:", format(signif(x$loglik, digits)),
    "(df =", paste(x$df, ")", sep = ""), "\n")
  cat("Number of iterations in BFGS optimization:", x$iterations, "\n\n")
  invisible(x)
}

plot.raschmodel <- function (x, type = c("profile", "curves", "regions", "information", "piplot"), ...)
{
   ## check input
  type <- match.arg(type)

  ## just call requested plotting function and pass arguments on
  switch(type,
         "curves" = curveplot(x, ...),
         "regions" = regionplot(x, ...),
         "profile" = profileplot(x, ...),
         "information" = infoplot(x, ...),
         "piplot" = piplot(x, ...))
}

predict.raschmodel <- function (object, newdata = NULL, type = c("probability",
  "cumprobability", "mode", "median", "mean", "category-information",
  "item-information", "test-information"), ref = NULL, ...)
{
  ## check type, process newdata, if NULL, use person parameters of given model object
  type <- match.arg(type)
  if (is.null(newdata)) {
    rs <- rowSums(object$data, na.rm = TRUE)
    rs <- rs[0 < rs & rs < ncol(object$data)]
    newdata <- personpar(object, vcov = FALSE)[rs]
    names(newdata) <- NULL
  }

  ## compute probabilities with internal function
  probs <- prm(theta = newdata, beta = itempar(object, ref = ref, vcov = FALSE))

  ## add category zero probabilities for consistency with other predict functions
  if (type %in% c("probability", "cumprobability", "category-information",
                  "item-information", "test-information")) {
    clnms <- colnames(probs)
    rwnms <- rownames(probs)
    nc <- ncol(probs)
    probs0 <- matrix(0, ncol = 2 * nc, nrow = nrow(probs))
    probs0[, seq(from = 2, by = 2, length.out = nc)] <- probs
    probs0[, seq(from = 1, by = 2, length.out = nc)] <- 1 - probs
    if (type == "cumprobability") {
      probs0[, seq(from = 1, by = 2, length.out = nc)] <- probs0[, seq(from = 1, by = 2, length.out = nc)] + probs
    }
    probs <- probs0
    rownames(probs) <- rwnms
    colnames(probs) <- as.vector(t(outer(clnms, c("C0", "C1"), paste, sep = if (type == "probability") "-" else ">=")))
  }

  ## if requested: compute test/item/category information (see Muraki, 1993, for details and further references)
  if (grepl("information", type)) {
    m <- length(clnms)
    if (type == "category-information") {
      info <- matrix(NA, nrow = nrow(probs), ncol = ncol(probs))
      colnames(info) <- as.vector(t(outer(clnms, c("C0", "C1"), paste, sep = "-")))
    } else {
      info <- matrix(NA, nrow = nrow(probs), ncol = m)
      colnames(info) <- clnms
    }

    for (j in 1:m) {
      idx <- grepl(clnms[j], colnames(probs))
      iteminfo <- apply(probs[, idx, drop = FALSE], 1, prod)
      if (type == "category-information") {
        info[, idx] <- probs[, idx, drop = FALSE] * iteminfo
      } else {
        info[, j] <- iteminfo
      }
    }
  }

  ## return as requested in type, for RM mode, median, mean is the same
  switch(type,
         "probability" = probs,
         "cumprobability" = probs,
         "mode" = round(probs),
         "median" = round(probs),
         "mean" = round(probs),
         "category-information" = info,
         "item-information" = info,
         "test-information" = matrix(rowSums(info), ncol = 1))
}

itempar.raschmodel <- function (object, ref = NULL, alias = TRUE, vcov = TRUE, ...)
{
  ## extract cf and labels, include restricted parameter
  cf <- c(0.00, coef(object))
  m <- length(cf)
  lbs <- names(cf)
  lbs[1] <- if(lbs[1] == "Item02") "Item01" else colnames(object$data)[1]

  ## process ref
  if (is.null(ref)) {
    ref <- 1:m
  } else if (is.vector(ref) && is.character(ref)) {
    stopifnot(all(ref %in% lbs))
    ref <- which(lbs %in% ref)
  } else if (is.vector(ref) && is.numeric(ref)) {
    ref <- as.integer(ref)
    stopifnot(all(ref %in% 1:m))
  } else if (is.matrix(ref) && is.numeric(ref)) {
    stopifnot(nrow(ref) == m && ncol(ref) == m)
  } else stop("Argument 'ref' is misspecified (see ?itempar for possible values).")

  ## if not given, specify contrast matrix
  if (is.matrix(ref)) {
    D <- ref
  } else {
    D <- diag(m)
    D[, ref] <- D[, ref] - 1/length(ref)
  }

  ## apply ref
  ## if vcov requested: adjust existing vcov
  cf <- as.vector(D %*% cf)
  if (vcov) {
    vc <- D %*% rbind(0, cbind(0, vcov(object))) %*% t(D)
  } else { # else return NA matrix
    vc <- matrix(NA, nrow = m, ncol = m)
  }

  ## set labels
  names(cf) <- rownames(vc) <- colnames(vc) <- lbs

  ## items solved by no or all subjects
  cf[object$items == "0"] <- -Inf
  cf[object$items == "1"] <- Inf

  ## process argument alias
  if (!alias) {
      if (is.matrix(ref)) {
          ## FIXME: Implement alias when ref is a specific constrast matrices -> detect linear dependent columns?
          stop("Processing of argument 'alias' not implemented with a contrast matrix given in argument 'ref'.")
      } else {
          cf <- cf[-ref[1]]
          vc <- vc[-ref[1], -ref[1]]
          alias <- paste0("I", ref[1])
          names(alias) <- lbs[ref[1]]
      }
  }
  
  ## setup and return result object
  rv <- structure(cf, class = "itempar", model = "RM", ref = ref, alias = alias, vcov = vc)
  return(rv)
}

threshpar.raschmodel <- function (object, type = c("mode", "median", "mean"), ref = NULL,
                                  alias = TRUE, relative = FALSE, cumulative = FALSE, vcov = TRUE, ...)
{
  ## check input
  type <- match.arg(type)

  ## extract relevant informations
  cf <- c(0.00, coef(object))
  m <- length(cf)
  ilbs <- names(cf)
  ilbs[1] <- if(ilbs[1] == "Item02") "Item01" else colnames(object$data)[1]
  clbs <- "C1"
  lbs <- paste(ilbs, clbs, sep = "-")

  ## process argument relative, type is not relevant for RM because mode = median = mean
  if (relative) {
    
    ## just check values given in ref
    if (is.null(ref)) {
      ref <- 1
    } else if (is.vector(ref) && is.character(ref)) {
      stopifnot(all(ref %in% clbs))
    } else if (is.vector(ref) && is.numeric(ref)) {
      stopifnot(all(as.integer(ref) %in% 1))
    } else if (is.matrix(ref) && is.numeric(ref)) {
      stopifnot(nrow(ref) == m && ncol(ref) == m)
    } else stop("Argument 'ref' is misspecified (see ?threshpar for possible values).")

    ## setup threshold parameters, if ref is a matrix, apply it.
    tp <- rep.int(0L, m)
    if (is.matrix(ref)) tp <- as.vector(ref %*% tp)
    tp <- as.list(tp)

    ## create vcov if requested
    if (vcov) {
      vc <- matrix(0L, nrow = m, ncol = m)
    } else {
      vc <- matrix(NA, nrow = m, ncol = m)
    }

  } else {

    ## process ref
    if (is.null(ref)) {
      ref <- 1:m
    } else if (is.vector(ref) && is.character(ref)) {
      stopifnot(all(ref %in% lbs))
      ref <- which(lbs %in% ref)
    } else if (is.vector(ref) && is.numeric(ref)) {
      ref <- as.integer(ref)
      stopifnot(all(ref %in% 1:m))
    } else if (is.matrix(ref)) {
      stopifnot(nrow(ref) == m && ncol(ref) == m)
    } else stop("Argument 'ref' is misspecified (see ?threshpar for possible values).")

    ## if not given, specify contrast matrix
    if (is.matrix(ref)) {
      D <- ref
    } else {
      D <- diag(m)
      D[, ref] <- D[, ref] - 1/length(ref)
    }
    
    ## impose contrast
    tp <- as.list(D %*% cf)

    ## create vcov if requested
    if (vcov) {
      vc <- D %*% rbind(0, cbind(0, vcov(object))) %*% t(D)
    } else {
      vc <- matrix(NA, nrow = m, ncol = m, dimnames = list(lbs, lbs))
    }

  }

  ## set labels
  names(tp) <- ilbs
  for (i in 1:m) names(tp[[i]]) <- "C1"
  rownames(vc) <- colnames(vc) <- lbs

  ## process argument alias
  if (!alias) {
      if (is.matrix(ref)) {
          ## FIXME: Implement alias when ref is a specific constrast matrices -> detect linear dependent columns?
          stop("Processing of argument 'alias' not implemented with a contrast matrix given in argument 'ref'.")
      } else {
          if (relative) { ## just return empty list/matrix as there are not unaliased parameters in the RM
              tp <- vector(mode = "list", length = m)
              for (i in 1:m) tp[[i]] <- numeric()
              vc <- matrix(0, ncol = 0, nrow = 0)
              alias <- as.list(rep(1, m))
              names(alias) <- ilbs
          } else {
              tp <- tp[-ref[1]]
              vc <- vc[-ref[1], -ref[1]]
              alias <- paste0("I", ref[1], "-C", ref[1])
              names(alias) <- ilbs[ref[1]]
          }
      }
  }

  ## return requested threshold parameters
  rv <- structure(tp, class = "threshpar", model = "RM", type = type, ref = ref, relative = relative, cumulative = cumulative, alias = alias, vcov = vc)
  return(rv)

}

discrpar.raschmodel <- function (object, ref = NULL, alias = TRUE, vcov = TRUE, ...)
{
  ## check input
  if (!is.null(ref)) warning("Argument 'ref' is currently not processed.") ## FIXME: Implement argument ref for discrimination parameters

  ## extract labels and number of items
  lbs <- c("", names(coef(object)))
  lbs[1] <- if(lbs[2] == "Item02") "Item01" else colnames(object$data)[1]
  m <- length(lbs)

  ## process argument alias
  if (alias) {
      dp <- rep.int(1, m)
      vc <- if (vcov) matrix(0, nrow = m, ncol = m, dimnames = list(lbs, lbs)) else matrix(NA, nrow = m, ncol = m, dimnames = list(lbs, lbs))
  } else {
      dp <- numeric()
      vc <- matrix(0, nrow = 0, ncol = 0)
      alias <- rep.int(1, m)
      names(alias) <- lbs
  }

  ## setup and return result object
  rv <- structure(dp, .Names = if (is.logical(alias)) lbs else NULL, class = "discrpar", model = "RM", ref = ref, alias = alias, vcov = vc)
  return(rv)
}

personpar.raschmodel <- function (object, ref = NULL, vcov = TRUE, interval = NULL, tol = 1e-8, ...)
{
  ## extract item parameters and relevant informations
  ip <- itempar(object, ref = ref, vcov = FALSE)
  m <- length(ip)
  rng <- 1:(m-1)

  ## calculate person parameters
  ## iterate over raw scores (from 1 until rmax-1), see Fischer & Molenaar, 1995, p.55
  if(is.null(interval)) interval <- c(-1, 1) * qlogis(1/m * 1e-3) #FIXME: 1e3 enough?
  pp <- sapply(rng, function(rawscore) uniroot(function(pp) rawscore - sum(plogis(pp - ip)),
    interval = interval, tol = tol)$root)

  if (vcov) {    
    ## relevant data for joint loglik approach
    y <- object$data[weights(object) > 0, , drop = FALSE]
    cs <- colSums(y)
    rs <- rowSums(y)
    rf <- tabulate(rs, nbins = m - 1)

    ## remove unidentified parameters
    rs <- rs[rf != 0]
    rng <- rng[rf != 0]
    pp <- pp[rf != 0]
    rf <- rf[rf != 0]

    ## obtain Hessian from objective function: joint log likelihood
    cloglik <- function (pp)  - sum(rf * rng * pp) + sum(cs * ip) + sum(rf * rowSums(log(1 + exp(outer(pp, ip, "-")))))
    vc <- solve(optim(pp, fn = cloglik, hessian = TRUE, method = "BFGS", control = list(reltol = tol, maxit = 0, ...))$hessian)
  } else {
    vc <- matrix(NA_real_, nrow = length(rng), ncol = length(rng))
  }
  colnames(vc) <- rownames(vc) <- rng

  ## setup and return result object
  rv <- structure(pp, .Names = rng, class = "personpar", model = "RM", vcov = vc)
  return(rv)

}

nobs.raschmodel <- function (object, ...)
{
  return(object$n)
}

bread.raschmodel <- function(x, ...) x$vcov * x$n

estfun.raschmodel <- function(x, ...) {
  ## extract data and parameters of interest
  par <- x$coefficients
  esf <- x$elementary_symmetric_functions
  y <- x$data
  weights_orig <- weights(x)
  y <- y[weights_orig > 0, , drop = FALSE]
  weights <- weights_orig[weights_orig > 0]
  rs <- rowSums(y)
  
  ## analytical gradient
  if(!x$na) {
    agrad <- weights * (- y + esf[[2]][rs + 1, , drop = FALSE] / esf[[1]][rs + 1])[,-1, drop = FALSE]
  } else {
    ## set up return value
    n <- nrow(y)
    k <- ncol(y)
    agrad <- matrix(0, nrow = n, ncol = k)

    ## observed NA patterns
    na_patterns <- factor(apply(is.na(y), 1, function(z) paste(which(z), collapse = "\r")))

    ## loop over observed NA patterns	   
    for(i in seq_along(levels(na_patterns))) {
      ## parse NA pattern
      lev_i <- levels(na_patterns)[i]
      na_i <- which(na_patterns == lev_i)
      wi_i <- as.integer(strsplit(lev_i, "\r")[[1]])
      wi_i <- if(length(wi_i) < 1) 1:k else (1:k)[-wi_i]

      ## compute gradient per pattern
      esf_i <- esf[[i]]
      rs_i <- rowSums(y[na_i, wi_i, drop = FALSE])
      agrad[na_i, wi_i] <- weights[na_i] * (- y[na_i, wi_i, drop = FALSE] +
    	esf_i[[2]][rs_i + 1, , drop = FALSE] / esf_i[[1]][rs_i + 1])
    }

    agrad <- agrad[, -1, drop = FALSE]
  }

  ## collect and return
  grad <- matrix(0, ncol = length(par), nrow = length(weights_orig))
  grad[weights_orig > 0,] <- agrad
  return(grad)
}


### misc. internal functions

## prm: calculate response probabilities for given thetas and betas under the RM.
prm <- function(theta = NULL, beta = NULL)
{
  ## check input
  stopifnot(!is.null(theta) && !is.null(beta))
  
  ## if list input, recurse...
  if (is.list(theta)) return(lapply(theta, prm, beta = beta))
  if (is.list(beta)) return(lapply(beta, prm, theta = theta))

  ## calculate probabilities
  return(plogis(outer(theta, beta, "-")))
}

## rrm: calculate response matrices for given thetas and betas under the RM.
rrm <- function(theta = NULL, beta = NULL, return_setting = TRUE)
{
  ## check input
  stopifnot(!is.null(theta) && !is.null(beta))
  
  ## if list input, recurse...
  if (is.list(theta)) return(lapply(theta, rrm, beta = beta, return_setting = return_setting))
  if (is.list(beta)) return(lapply(beta, rrm, theta = theta, return_setting = return_setting))

  ## calculate response probabilities and responses (randomized cutpoint like in eRm:::sim.rasch)
  n <- length(theta)
  m <- length(beta)
  probs <- prm(theta = theta, beta = beta)
  resp <- (matrix(runif(n * m), nrow = n, ncol = m) < probs) + 0

  ## return
  if (return_setting) list(theta = theta, beta = beta, data = resp) else resp
}
