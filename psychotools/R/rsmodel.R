## Function to fit a rating scale modell, see e.g. Andrich (1978), Andersen (1995)
rsmodel <- function (y, weights = NULL, start = NULL, reltol = 1e-10,
  deriv = c("sum", "diff"), hessian = TRUE, maxit = 100L, full = TRUE, ...)
{
  ## argument matching
  deriv <- match.arg(deriv)

  ## process input and data.
  diff <- deriv == "diff"
  y <- as.matrix(y)
  n_org <- nrow(y)
  m <- ncol(y)
  ident_items <- rep.int(TRUE, m)

  ## process weights
  if (is.null(weights)) weights <- rep.int(1L, n_org) else stopifnot(length(weights) == n_org)
  y <- y[weights > 0, , drop = FALSE]
  weights_org <- weights
  weights <- weights[weights > 0]
  n <- nrow(y)  

  ## move categories to zero if necessary, get number of categories per item
  mincat <- apply(y, 2, min, na.rm = TRUE)
  if(all(mincat > 0)) {
    warning("The minimum score is above zero for all items. Items are scaled down to zero.")
    y[, mincat > 0] <- scale.default(y[, mincat > 0], center = mincat[mincat > 0], scale = FALSE)
  }
  oj <- apply(y, 2, max, na.rm = TRUE)

  ## remove unidentified items (only one category used)
  unid <- oj == 0
  nunid <- sum(unid)
  if (nunid) {
    y <- y[, !unid, drop = FALSE]
    oj <- oj[!unid]
    m <- m - nunid
    ident_items[unid] <- FALSE
    warning("There were unidentified items (only one category used, ", paste("I", which(unid), sep = "", collapse = ", "), ").")
  }

  ## stop if not all oj equal or oj = 1 (dichotomous items).
  o <- max(oj)
  if (any(oj != o)) warning("Not all items have the same number of categories. Consider using pcmodel() for data with a variable number of categories per item.")
  oj <- rep.int(o, m)
  if (o == 1) stop("Rating scale model requires polytomous items. Use raschmodel() for dichotomous items.")

  ## check for 'total' null categories (null for all items), if present, stop, because tau is assessable
  cat_freq <- apply(y + 1, 2, tabulate, nbins = o + 1)
  if (any(rowSums(cat_freq) == 0)) stop("There are categories which are null for all items. ",
                   "Please use pcmodel() instead (see argument nullcats in ?pcmodel).")

  ## check for missing
  y_na <- is.na(y)
  any_na <- any(y_na)

  ## set starting values
  if (is.null(start)) start <- rep.int(0, m + o - 2)

  ## small helper function: conversion of a given rating scale parameter vector (and additional information)
  ## to a pcm parameter vector (item-category-parameters) (expects that beta_0 = 0 and tau_0 = 0.)
  rsm2pcm <- function(rsm_par, item_par_pos, cat_par_pos, oj, cat_scores) rep.int(c(0, rsm_par[item_par_pos]), oj) * cat_scores + c(0, rsm_par[cat_par_pos])

  ## calculate positions and helper variables: ipar/cpar without and mv/ov with restricted parameters
  ipar <- 1:(m-1)
  cpar <- m:(m + o - 2)
  mv <- 1:m
  ov <- 1:o

  ## calculate item, category and person totals
  itot <- colSums(y * weights)[-1]                                   # without item 1
  cs <- t(apply(y, 1, tabulate, nbins = o))
  ctot <- colSums(cs * weights)[-1] # without category 0 and 1
  rs <- rowSums(y, na.rm = TRUE)
  ptot <- as.vector(tapply(weights, factor(rs, levels = 0:sum(oj)), sum))
  ptot[is.na(ptot)] <- 0
  
  ## build log-likelihood
  if(!any_na) {

    ## objective function: conditional log-likelihood
    cloglik <- function (par) {

      ## transform rsm par to pcm par list for esf function
      esf_par <- split(rsm2pcm(par, ipar, cpar, oj, ov), rep.int(mv, oj))

      ## finally: the cll
      cll <- - sum(itot * par[ipar]) - sum(ctot * par[cpar]) - sum(ptot * log(elementary_symmetric_functions(par = esf_par, order = 0, diff = diff)[[1]]))

      ## catch degenerated cases (typically caused by non-finite gamma)
      if (is.na(cll) | !is.finite(cll)) cll <- -.Machine$double.xmax
      
      return(-cll)
    }
    
    ## static helper variables for analytical gradient
    betaindex <- rep.int(mv,  rep.int(o, m))
    tauindex <- rep.int(ov, m)

    ## analytical gradient
    agrad <- function (par) {
        
        ## transform RSM parameters to PCM parameters and calculate actual ESF
        parx <- rsm2pcm(par, ipar, cpar, oj, ov)
        esf <- elementary_symmetric_functions(par = split(parx, rep.int(mv, oj)), order = 1, diff = diff)
        gamma0 <- esf[[1]][rs + 1]
        gamma1_pcm <- esf[[2]]

        ## calculate transformed derivatives, select relevant derivatives with ptot, drop unindentified parameters.
        gamma1 <- matrix(0, nrow = nrow(gamma1_pcm), ncol = m + o)
        for (j in mv) gamma1[, j] <- gamma1_pcm[, betaindex == j, drop = FALSE] %*% ov
        for (k in ov) gamma1[, k + m] <- rowSums(gamma1_pcm[, tauindex == k, drop = FALSE])
        gamma1 <- apply(gamma1, 2, "[", rs + 1)

        ## finally: the gradient
        agrad <- matrix(0, nrow = n, ncol = m + o)
        agrad[, mv] <- weights * (- y + (gamma1 / gamma0)[, mv, drop = FALSE])
        agrad[, m + ov] <- weights * (- cs + (gamma1 / gamma0)[, m + ov, drop = FALSE])

        ## return aggegrated gradient
        return(- colSums(agrad[, -c(1, m + 1), drop = FALSE]))
    }

  } else {

    ## fetch different NA-patterns like Achim does
    na_patterns <- factor(apply(y_na, 1, function(z) paste(which(z), collapse = "\r")))
    lev_na_patterns <- levels(na_patterns)
    par_pos <- rep.int(mv, oj) ## calculate positions of pcm parameters (needed for ESF)

    ## setup na_pattern lists for various elements of loglik
    ptot_i <- itot_i <- ctot_i <- ipar_i <- oj_i <- pcm_par_i <- m_i <- mv_i <- vector("list", length(lev_na_patterns))
    na_obs_i <- rs_i <- weights_i <- vector("list", length(lev_na_patterns))

    ## loop over observed NA patterns, calculate fix things once and store in list
    for(i in seq_along(lev_na_patterns)) {
      
      ## from pattern i: get NA item(s), observations of and weights with this pattern
      na_items_i <- as.integer(strsplit(lev_na_patterns[i], "\r")[[1]])
      n_na_items_i <- length(na_items_i)
      na_obs_i[[i]] <- which(na_patterns == lev_na_patterns[i])
      weights_i[[i]] <- weights[na_obs_i[[i]]]
      
      ## select subset
      if(n_na_items_i < 1) {            # no missings
        y_i <- y[na_obs_i[[i]], , drop = FALSE]
        m_i[[i]] <- m
        mv_i[[i]] <- mv
        pcm_par_i[[i]] <- rep.int(TRUE, m*o)
        ipar_i[[i]] <- ipar
        oj_i[[i]] <- oj
      } else {                          # specific NA pattern
        y_i <- y[na_obs_i[[i]], -na_items_i, drop = FALSE]
        rs_i[[i]] <- rowSums(y_i)
        m_i[[i]] <- m - n_na_items_i
        mv_i[[i]] <- mv[-na_items_i]
        pcm_par_i[[i]] <- !(par_pos%in% na_items_i)
        ipar_i[[i]] <- if ((n_na_items_i == 1) && na_items_i == 1) ipar else ipar[-(setdiff(na_items_i, 1) - 1)] # beta_j is in col j-1
        oj_i[[i]] <- oj[-na_items_i]
    }

      ## calculate category and person totals for NA-group i
      itot_i[[i]] <- colSums(y_i * weights_i[[i]])
      if (!(1 %in% na_items_i)) itot_i[[i]] <- itot_i[[i]][-1]                           # remove item 1 if not already removed
      ctot_i[[i]] <- colSums(t(apply(y_i, 1, tabulate, nbins = o)) * weights_i[[i]])[-1] # without category 0 and 1
      rs_i[[i]] <- rowSums(y_i)
      ptot_i[[i]] <- as.vector(tapply(weights_i[[i]], factor(rs_i[[i]], levels = 0:sum(oj_i[[i]])), sum))
      ptot_i[[i]][is.na(ptot_i[[i]])] <- 0
    }

    ## fun for mapply (calculates cll contributions per na_group with variable par's)
    cll_i <- function (itot_i, item_par_i, ctot_i, ptot_i, esf_par_i, cat_par) {
      - sum(itot_i * item_par_i) - sum(ctot_i * cat_par) - sum(ptot_i * log(elementary_symmetric_functions(par = esf_par_i, order = 0, diff = diff)[[1]]))
    }

    ## objective function: conditional log-likelihood (build incrementally by looping over the NA-patterns
    cloglik <- function (par) {

      ## fetch parameter vectors
      esf_par <- rsm2pcm(par, ipar, cpar, oj, ov)
      esf_par_i <- lapply(pcm_par_i, function(x) esf_par[x])
      esf_par_i <- mapply(split, esf_par_i, mapply(function (x, y) rep.int(1:x, y), m_i, oj_i))
      item_par_i <- lapply(ipar_i, function(x) par[x])
      cat_par <- par[cpar]

      ## conditional log-likelihood
      cll <- sum(mapply(cll_i, itot_i, item_par_i, ctot_i, ptot_i, esf_par_i,
                        MoreArgs = list(cat_par = cat_par), SIMPLIFY = TRUE, USE.NAMES = FALSE))

      ## catch degenerated cases (typically cause by non-finite gamma)
      if(is.na(cll) | !is.finite(cll)) cll <- -.Machine$double.xmax

      return(-cll)
    }

    ## static helper variables for analytical gradient
    betaindex_i <- mapply(function (x, y, ...) rep.int(x, rep.int(o, y)), mv_i, m_i, MoreArgs = list(o = o), SIMPLIFY = FALSE)
    tauindex_i <- lapply(m_i, function (x) rep.int(ov, x))
    
    ## analytical gradient
    agrad <- function (par) {

        ## transform RSM parameters to PCM parameters
        parx <- rsm2pcm(par, ipar, cpar, oj, ov)
        esf_par_i <- lapply(pcm_par_i, function (x) parx[x])
        esf_par_i <- mapply(split, esf_par_i, mapply(function (x, y) rep.int(1:x, y), m_i, oj_i))
        esf_i <- mapply(elementary_symmetric_functions, par = esf_par_i, MoreArgs = list(order = 1, diff = diff), SIMPLIFY = FALSE)

        ## loop over obseved NA patterns and gradually built up analytical gradient
        grad <- matrix(0, nrow = n, ncol = m + o)
        for (i in seq_along(lev_na_patterns)) {

            ## fetch ESF with PCM parametrization
            gamma0_i <- esf_i[[i]][[1]][rs_i[[i]] + 1]
            gamma1_pcm_i <- esf_i[[i]][[2]]

            ## transform derivatives via delta rule to RSM parametrization
            gamma1_i <- matrix(0, nrow = nrow(gamma1_pcm_i), ncol = m_i[[i]] + o)
            for (j in 1:m_i[[i]]) gamma1_i[, j] <- gamma1_pcm_i[, betaindex_i[[i]] == mv_i[[i]][j], drop = FALSE] %*% ov
            for (k in ov) gamma1_i[, k + m_i[[i]]] <- rowSums(gamma1_pcm_i[, tauindex_i[[i]] == k, drop = FALSE])
            gamma1_i<- apply(gamma1_i, 2, "[", rs_i[[i]] + 1)
            if (!is.matrix(gamma1_i)) gamma1_i<- matrix(gamma1_i, nrow = 1)

            ## finally: the gradient for NA group i
            grad[na_obs_i[[i]], mv_i[[i]]] <- weights_i[[i]] * (- y[na_obs_i[[i]], mv_i[[i]], drop = FALSE] + (gamma1_i/gamma0_i)[, 1:m_i[[i]], drop = FALSE])
            grad[na_obs_i[[i]], m + ov] <- weights_i[[i]] * (- cs[na_obs_i[[i]], , drop = FALSE] + (gamma1_i/gamma0_i)[, m_i[[i]] + ov, drop = FALSE])

        }

        return(- colSums(grad[, -c(1, m + 1), drop = FALSE]))
    }
  }

  ## optimization
  opt <- optim(par = start, fn = cloglik, gr = agrad, method = "BFGS",
               hessian = hessian, control = list(reltol = reltol, maxit = maxit, ...))
 
  ## final estimates ...
  est <- opt$par
  names(est) <- if (is.null(colnames(y))) {
    c(paste("I", setdiff(which(ident_items), 1), sep = ""), paste("C", 2:o, sep = ""))
  } else c(colnames(y)[setdiff(which(ident_items), 1)], paste("C", 2:o, sep = ""))
    

  ## ... and (if requested) esf of these ...
  if (full) {
    parx <- rsm2pcm(est, ipar, cpar, oj, ov)
    esf <- if (any_na) {
      esf_par_i <- lapply(pcm_par_i, function (x) parx[x])
      esf_par_i <- mapply(split, esf_par_i, mapply(function (x, y) rep.int(1:x, y), m_i, oj_i))
      mapply(elementary_symmetric_functions, par = esf_par_i, MoreArgs = list(order = 1, diff = diff), SIMPLIFY = FALSE)
    } else {
      elementary_symmetric_functions(par = split(parx, rep.int(mv, oj)), order = 1, diff = diff)
    }
    
    ## ... as well as variance-bcovariance matrix
    if (hessian) {
      vc <- opt$hessian
      vc <- solve(vc)
    } else {
      vc <- matrix(NA, nrow = length(est), ncol = length(est))
    }
    rownames(vc) <- colnames(vc) <- names(est)
  } else {
    esf <- NULL
    vc <- NULL
  }


  ## collect results, set class, and return
  res <- list(coefficients = est,
              vcov = vc,
              data = y,
              items = ident_items,
              categories = oj,
              n = sum(weights_org > 0),
              n_org = n_org,
              weights = if (identical(as.vector(weights_org), rep.int(1L, n_org))) NULL else weights_org,
              na = any_na,
              esf = esf,
              loglik = -opt$value,
              df = length(est),
              code = opt$convergence,
              iterations = tail(na.omit(opt$counts), 1L),
              reltol = reltol)
  class(res) <- c("rsmodel")
  return(res )
}

print.rsmodel <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  cat("RSM item and threshold parameters:\n")
  print(coef(x), digits = digits)
  invisible(x)
}

coef.rsmodel <- function (object, ...) object$coefficients

vcov.rsmodel <- function (object, ...) object$vcov

logLik.rsmodel <- function (object, ...) structure(object$loglik, df = object$df, class = "logLik")

weights.rsmodel <- function (object, ...) if (is.null(object$weights)) rep.int(1L, object$n_org) else object$weights

summary.rsmodel <- function (object, vcov. = NULL, ...) {
  ## coefficients
  cf <- coef(object)

  ## covariance matrix
  if (is.null(vcov.)) 
    vc <- vcov(object)
  else {
    if (is.function(vcov.)) vc <- vcov.(object)
    else vc <- vcov.
  }
  
  ## Wald test of each coefficient
  cf <- cbind(cf, sqrt(diag(vc)), cf/sqrt(diag(vc)), 2 * pnorm(-abs(cf/sqrt(diag(vc)))))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

  object$coefficients <- cf      
  class(object) <- "summary.rsmodel"
  return(object)
}


print.summary.rsmodel <- function (x, digits = max(3, getOption("digits") - 3),
                                   signif.stars = getOption("show.signif.stars"), ...)
{
  if (is.null(x$call)) {
    cat("\nRating scale model\n\n")  
  } else {
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  }

  if (any(!x$items)) cat("Excluded items:",
                         paste(names(x$items)[!x$items], collapse = ", "), "\n\n")

  cat("Item location and threshold parameters:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)

  cat("\nLog-likelihood:", format(signif(x$loglik, digits)),
      "(df =", paste(x$df, ")", sep = ""), "\n")
  cat("Number of iterations in BFGS optimization:", x$iterations, "\n\n")
  invisible(x)
}

plot.rsmodel <- function (x, type = c("regions", "profile", "curves", "information", "piplot"), ...)
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

predict.rsmodel <- function (object, newdata = NULL, type = c("probability",
  "cumprobability", "mode", "median", "mean", "category-information",
  "item-information", "test-information"), ref = NULL, ...)
{
  ## check type, process newdata, if NULL, use person parameters of given model object
  type <- match.arg(type)
  if (is.null(newdata)) {
    rs <- rowSums(object$data, na.rm = TRUE)
    rs <- rs[0 < rs & rs < (max(object$data) * ncol(object$data))]
    newdata <- personpar(object, vcov = FALSE)[rs]
    names(newdata) <- NULL
  }
  nms <- names(newdata)

  ## itempars and raw probabilities
  ip <- itempar(object, ref = ref, vcov = FALSE)
  tp <- threshpar(object, type = "mode", ref = ref, relative = TRUE, vcov = FALSE)[[1]]
  ilbs <- names(ip)
  os <- 0:length(tp)
  o <- max(os) + 1
  probs <- prsm(theta = newdata, beta = ip, tau = tp)

  ## if requested: compute test/item/category information (see Muraki, 1993, for details and further references)
  if (grepl("information", type)) {
    m <- length(ilbs)
    if (type == "category-information") {
      info <- matrix(NA, nrow = length(newdata), ncol = m * o)
      colnames(info) <- as.vector(t(outer(ilbs, os, paste, sep = "-C")))
    } else {
      info <- matrix(NA, nrow = length(newdata), ncol = m)
      colnames(info) <- ilbs
    }

    for (j in 1:m) {
      ojm <- matrix(rep(os, each = length(newdata)), nrow = length(newdata), ncol = ncol(probs[[j]]))
      iteminfo <- rowSums((ojm - rowSums(ojm * probs[[j]]))^2 * probs[[j]])
      if (type == "category-information") {
        info[, (j - 1) * o + os + 1] <- probs[[j]] * iteminfo
      } else {
        info[, j] <- iteminfo
      }
    }
  }
  

  ## return as requested in type, for RM mode, median, mean is the same
  switch(type,
         "probability" = {
           ## construct labels
           lbs <- as.vector(t(outer(ilbs, os, paste, sep = "-C")))
           ## create (named) matrix with probabilities
           probs <- do.call("cbind", probs)
           dimnames(probs) <- list(nms, lbs)
           probs
           },
         "cumprobability" = {
           ## construct labels
           lbs <- as.vector(t(outer(ilbs, os, paste, sep = ">=C")))
           ## create (named) matrix with probabilities
           probs <- lapply(probs, function (probsj) t(apply(probsj, 1, function (x) rev(cumsum(rev(x))))))
           probs <- do.call("cbind", probs)
           dimnames(probs) <- list(nms, lbs)
           probs
           },
         "mode" = sapply(probs, apply, 1, which.max) - 1,
         "median" = sapply(probs, function (probsj) apply(probsj, 1, which.max) - 1),
         "mean" = round(sapply(probs, function (probsj) apply(probsj, 1, function (j) sum(j * 0:(length(j) - 1))))),round(probs),
         "category-information" = info,
         "item-information" = info,
         "test-information" = matrix(rowSums(info), ncol = 1))
}

itempar.rsmodel <- function (object, ref = NULL, alias = TRUE, vcov = TRUE, ...)
{
  ## extract cf and labels, include restricted parameters
  m <- sum(object$items)
  ip <- c(0.00, coef(object)[1:(m - 1)])
  ctp <- c(0.00, coef(object)[-(1:(m - 1))])
  cf <- c(ip, ctp)
  o <- length(ctp)  
  lbs <- names(ip)
  lbs[1] <- if (lbs[2] == "I2") "I1" else colnames(object$data)[1]

  ## process ref
  if (is.null(ref)) {
    ref <- 1:m
  } else if (is.vector(ref) && is.character(ref)) {
    stopifnot(all(ref %in% lbs))
    ref <- which(lbs %in% ref)
  } else if (is.vector(ref) && is.numeric(ref)) {
    ref <- as.integer(ref)
    stopifnot(all(ref %in% 1:m))
  } else if (is.matrix(ref) && (is.numeric(ref))) {
    stopifnot(nrow(ref) == m & ncol(ref) == m)
  } else stop("Argument 'ref' is misspecified (see ?itempar for possible values).")

  ## if not given, specify contrast matrix from ref
  if (is.matrix(ref)) {
    D <- ref
  } else {
    D <- diag(m)
    D[, ref] <- D[, ref] - 1/length(ref)
  }
  
  ## transform estimated item-specific parameters and estimated
  ## cumulative threshold parameters to mean absolute threshold parameters ...
  C <- matrix(0, nrow = m, ncol = m + o)
  C[1:m, 1:m] <- diag(m)
  C[, m + o] <- 1/o
  cf <- as.vector(C %*% cf)

  ## if requested: create adjusted vcov
  if (vcov) {
    vc <- matrix(0, nrow = m + o, ncol = m + o)
    vc[2:m, 2:m] <- vcov(object)[1:(m-1), 1:(m-1)]
    vc[2:m, (m + 2):(m + o)] <- vcov(object)[1:(m - 1), m:(m + o - 2)]
    vc[(m + 2):(m + o), 2:m] <- vcov(object)[m:(m + o - 2), 1:(m - 1)]
    vc[(m + 2):(m + o), (m + 2):(m + o)] <- vcov(object)[m:(m + o - 2), m:(m + o - 2)]
    vc <- C %*% vc %*% t(C)
  }
  
  ## finally apply ref
  ## if vcov requested: adjust existing vcov
  cf <- as.vector(D %*% cf)
  if (vcov) {
    vc <- D %*% vc %*% t(D)
  } else { # else return NA matrix
    vc <- matrix(NA, nrow = m, ncol = m)
  }
  
  ## set labels
  names(cf) <- rownames(vc) <- colnames(vc) <- lbs

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
  rv <- structure(cf, class = "itempar", model = "RSM", ref = ref, vcov = vc)
  return(rv)
}

threshpar.rsmodel <- function (object, type = c("mode", "median", "mean"), ref = NULL,
                               alias = TRUE, relative = FALSE, cumulative = FALSE, vcov = TRUE, ...)
{
  ## check input
  type <- match.arg(type)

  ## extract relevant informations
  m <- sum(object$items)
  ip <- c(0.00, coef(object)[1:(m - 1)])     # item specific parameters
  tp <- c(0.00, coef(object)[-c(1:(m - 1))]) # cumulative relative threshold parameters
  o <- length(tp)
  b <- m * o
  ilbs <- names(ip)
  ilbs[1] <- if (ilbs[2] == "I2") "I1" else colnames(object$data)[1]
  clbs <- paste0("C", 1:o)
  lbs <- as.vector(outer(clbs, ilbs, function (i, j) paste(j, i, sep = "-")))

  ## process argument relative
  if (relative) {

    if (type == "mode") {

      ## process ref
      if (is.null(ref)) {
        ref <- 1:o
      } else if (is.vector(ref) && is.character(ref)) {
        stopifnot(all(ref %in% clbs))
        ref <- which(clbs %in% ref)
      } else if (is.vector(ref) && is.numeric(ref)) {
        ref <- as.integer(ref)
        stopifnot(all(ref %in% 1:o))
      } else if (is.matrix(ref) && is.numeric(ref)) {
        stopifnot(nrow(ref) == o && ncol(ref) == o)
      } else stop("Argument 'ref' is misspecified (see ?threshpar for possible values).")

      ## transform estimated cumulative relative item threshold parameters to relative item threshold parameters
      C <- diag(o)
      for (i in 2:o) C[i, i - 1] <- -1
      tp <- as.vector(C %*% tp)
      ## if vcov requested: create adjusted vcov
      if (vcov) {
        vc0 <- rbind(0, cbind(0, vcov(object)))[-c(1:(m - 1)), -c(1:(m - 1))]
        vc0 <- C %*% vc0 %*% t(C)
      }

      ## if not given, specify contrast matrix
      if (is.matrix(ref)) {
        D <- ref
      } else {
        D <- diag(o)
        D[, ref] <- D[, ref] - 1/length(ref)
      }
      
      ## apply ref
      ## if vcov requested: adjust vcov too
      tp <- as.vector(D %*% tp)
      tp <- split(rep(tp, m), rep(1:m, each = o))
      if (vcov) {
        vc0 <- D %*% vc0 %*% t(D)
        vc <- matrix(0, nrow = m * o, ncol = m * o)
        for (i in 1:m) vc[1:o + (i - 1) * o, 1:o + (i - 1) * o] <- vc0
      } else { # else return NA matrix
        vc <- matrix(NA, nrow = m * o, ncol = m * o)
      }

      ## if cumulative relative threshold parameters are requested: transform relative item
      ## threshold parameters back to cumulative relative item threshold parameters
      if (cumulative) {
        tp <- lapply(tp, cumsum)
        if (vcov) {
          C <- matrix(0, nrow = o, ncol = o)
          for (k in 1:o) C[k, 1:k] <- 1
          for (i in 1:m) vc[1:o + (i - 1) * o, 1:o + (i - 1) * o] <- C %*% vc[1:o + (i - 1) * o, 1:o + (i - 1) * o] %*% t(C)
        }
      }
      
    } else stop("Relative threshold parameters not implemented for types other than mode.")
    
  } else {

    if (type == "mode") {

      ## process ref
      if (is.null(ref)) {
        ref <- 1:b
      } else if (is.vector(ref) && is.character(ref)) {
        stopifnot(all(ref %in% clbs))
        ref <- which(lbs %in% ref)
      } else if (is.vector(ref) && is.numeric(ref)) {
        ref <- as.integer(ref)
        stopifnot(all(ref %in% 1:b))
      } else if (is.matrix(ref) && is.numeric(ref)) {
        stopifnot(ncol(ref) == b && nrow(ref) == b)
      } else stop("Argument 'ref' is misspecified (see ?threshpar for possible values).")
      
      ## transform estimated item-specific parameters and estimated cumulative
      ## relative item threshold parameters to absolute item threshold parameters
      C <- matrix(0, nrow = b, ncol = m + o)
      for (j in 1:m) {
        C[(j - 1) * o + 1:o, j] <- 1
        C[(j - 1) * o + 1:o, m + 1:o] <- diag(o)
        for (k in 2:o) C[(j - 1) * o + k, m  + k - 1] <- -1
      }
      tp <- C %*% c(ip, tp)
      if (vcov) {
        vc <- matrix(0, nrow = m + o, ncol = m + o)
        vc[2:m, 2:m] <- vcov(object)[1:(m-1), 1:(m-1)]
        vc[2:m, (m + 2):(m + o)] <- vcov(object)[1:(m - 1), m:(m + o - 2)]
        vc[(m + 2):(m + o), 2:m] <- vcov(object)[m:(m + o - 2), 1:(m - 1)]
        vc[(m + 2):(m + o), (m + 2):(m + o)] <- vcov(object)[m:(m + o - 2), m:(m + o - 2)]
        vc <- C %*% vc %*% t(C)
      }
      
      ## if not given, specify contrast matrix
      if (is.matrix(ref)) {
        D <- ref
      } else {
        D <- diag(b)
        D[, ref] <- D[, ref] - 1/length(ref)
      }
      
      ## apply ref
      ## if vcov requested: adjust existing vcov
      tp <- as.vector(D %*% tp)
      tp <- split(tp, rep(1:m, each = o))
      if (vcov) {
        vc <- D %*% vc %*% t(D)
      } else { # else return NA matrix
        vc <- matrix(NA, nrow = b, ncol = b)
      }

      ## if cumulative absolute item threshold parameters are requested:
      ## cumulate calculated absolute item threshold parameters and adjust vcov if requested ...
      if (cumulative) {
        tp <- lapply(tp, cumsum)
        if (vcov) {
          C <- matrix(0, nrow = m * o, ncol = m * o)
          for (j in 1:m) {
            for (k in 1:o) C[(j - 1) * o + k, (j - 1) * o + 1:k] <- 1
          }
          vc <- C %*% vc %*% t(C)
        }
      }

    } else {

      ## process ref and setup NA vcov
      if (!is.null(ref)) warning("Argument 'ref' is not processed for types other than mode.")
      if (!alias) warning("Argument 'alias' is not processed for types other than mode.")
      tp <- lapply(ip, "+", tp)
      vc <- matrix(NA, nrow = b, ncol = b)

      if (type == "mean") {

        ## setup threshold parameters, expected scores, function to find locations on theta axis
        xpct <- 1:o - 0.5
        zexpct <- function (theta = NULL, delta = NULL, expct = NULL) ppcm(theta = theta, delta = delta) %*% 0:length(delta) - expct
        
        ## loop though items and find locations by means of zexpct() and uniroot()
        for (j in 1:m) tp[[j]] <- sapply(xpct, function (xp) uniroot(f = zexpct, interval = c(-10, 10), delta = tp[[j]], expct = xp)$root)

      }

      if (type == "median") {

        ## setup threshold parameters, expected scores, function to find locations on theta axis
        zmedian <- function (theta = NULL, delta = NULL, geq = NULL, ncat = NULL) rowSums(ppcm(theta = theta, delta = delta)[, (geq + 1):ncat, drop = FALSE]) - 0.5

        ## loop though items and find locations by means of zmedian() and uniroot()
        for (j in 1:m) tp[[j]] <- sapply(1:o, function (geq) uniroot(f = zmedian, interval = c(-10, 10), delta = tp[[j]], geq = geq, ncat = o + 1)$root)
        
      }

      ## if cumulative absolute item threshold parameters are requsted, just cumulate ...
      if (cumulative) tp <- lapply(tp, cumsum)
    }
    
  }

  ## set labels
  names(tp) <- ilbs
  rownames(vc) <- colnames(vc) <- lbs
  for (i in 1:m) names(tp[[i]]) <- paste0("C", 1:o)

  ## process argument alias
  if (!alias && type == "mode") {
    if (is.matrix(ref)) {
      ## FIXME: Implement alias when ref is a specific constrast matrices -> detect linear dependent columns?
      stop("Processing of argument 'alias' not implemented with a contrast matrix given in argument 'ref'.")
    } else {
      if (relative) {
        tp <- lapply(tp, function (j) j[-ref[1]])
        s <- -seq(from = ref[1], by = o, length.out = m)
        vc <- vc[-s, -s]
        alias <- as.list(rep(ref[1], m))
        names(alias) <- ilbs
      } else {
        ref1 <- ref[1]
        i <- split(1:b, rep(1:m, each = o))
        item <- which(sapply(i, function (j) ref1 %in% j))
        tp[[item]] <- tp[[item]][-which(ref1 == i[[item]])]
        vc <- vc[-ref1, -ref1]
        alias <- paste0("I", item, "-C", which(ref1 == i[[item]]))
        names(alias) <- ilbs[item]
      }
    }
  }

  ## setup and return result object
  rv <- structure(tp, class = "threshpar", model = "RSM", type = type, ref = ref, relative = relative, alias = alias, vcov = vc)
  return(rv)  
}

discrpar.rsmodel <- function (object, ref = NULL, alias = TRUE, vcov = TRUE, ...)
{
  ## check input
  if (!is.null(ref)) warning("Argument 'ref' is currently not processed.")  ## FIXME: Implement argument ref for discrimination parameters

  ## extract labels and number of items
  m <- sum(object$items)
  lbs <- c("", names(coef(object)[1:(m - 1)]))
  lbs[1] <- if (lbs[2] == "I2") "I1" else colnames(object$data)[1]

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
  rv <- structure(dp, .Names = if (is.logical(alias)) lbs, class = "discrpar", model = "RSM", ref = ref, alias = alias, vcov = vc)
  return(rv)
}

personpar.rsmodel <- function (object, ref = NULL, vcov = TRUE, interval = NULL, tol = 1e-8, ...)
{
  ## extract item parameters in pcm formulation
  tp <- threshpar(object, ref = ref, type = "mode", alias = TRUE, vcov = FALSE)
  tp <- lapply(tp, cumsum)
  m <- length(tp)
  o <- length(tp[[1]])
  os <- 1:o
  osv <- rep.int(os, m)
  rng <- 1:(m * o - 1)
  
  ## iterate over raw scores (from 1 until rmax-1), see Fischer & Molenaar, 1995, p.286, Eq. (15.43)
  if(is.null(interval)) interval <- c(-1, 1) * qlogis(1/m * 1e-3) #FIXME: 1e3 enough?
  pp <- sapply(rng, function(rawscore) uniroot(function(pp) rawscore - sum(osv * exp(osv * pp - unlist(tp)) / unlist(lapply(tp, function (beta) rep.int(1 + sum(exp(os * pp - beta)), o)))),
    interval = interval, tol = tol)$root)

  ## calculate person parameters
  if (vcov) {
    ## relevant data for joint loglik approach
    y <- object$data[weights(object) > 0, , drop = FALSE]
    rs <- rowSums(y)
    rf <- tabulate(rs, nbins = m * o - 1)
    cs <- apply(y, 2, tabulate, nbins = o)
    
    ## remove unidentified parameters
    rs <- rs[rf != 0]
    rng <- rng[rf != 0]
    pp <- pp[rf != 0]
    rf <- rf[rf != 0]

    ## transform ip to matrix (o x m)
    tp <- do.call("cbind", tp)

    ## obtain Hessian from objective function: joint log likelihood
    cloglik <- function (pp) - sum(rf * rng * pp) + sum(cs * tp) + sum(rf * colSums(log(1 + apply(outer(os, pp), 2, function (x) colSums(exp(x - tp))))))
    vc <- solve(optim(pp, fn = cloglik, hessian = TRUE, method = "BFGS", control = list(reltol = tol, maxit = 0, ...))$hessian)
  } else {
    vc <- matrix(NA, nrow = length(rng), ncol = length(rng))
  }
  colnames(vc) <- rownames(vc) <- rng
  
  ## setup and return result object
  rv <- structure(pp, .Names = rng, class = "personpar", model = "RSM", vcov = vc)
  return(rv)
}

nobs.rsmodel <- function (object, ...)
{
  return(object$n)
}

bread.rsmodel <- function(x, ...) x$vcov * x$n

estfun.rsmodel <- function(x, ...) {
  ## get relevant informations
  dat <- x$data                    # completely cleaned (downcoded, null cats treatment, weights) data.
  weights_org <- weights(x)
  weights <- weights_org[weights_org > 0]
  n <- nrow(dat)
  m <- ncol(dat)
  o <- max(dat, na.rm = TRUE)
  npar <- x$df
  ptot <- rowSums(dat, na.rm = TRUE) + 1 # +1 because gamma of score 0 is in row 1.
  ctot <- t(apply(dat, 1, tabulate, nbins = o))
  
  ## helper variables for first derivative calculation
  mv <- 1:m
  ov <- 1:o
  beta_index <- rep.int(mv,  rep.int(o, m))
  tau_index <- rep.int(ov, m)

  ## calculate gradient (matrix with 1:(m-1) cols for beta_j, j = 2, m, and o - 1 cols for tau^*_k, k = 2, ..o)
  if (!x$na) {
    
    ## select zero and first derivatives of gamma functions
    ## transform derivatives with item-category parameters to derivatives of RSM-parameters
    gamma0 <- x$esf[[1]][ptot]
    gamma1_pcm <- x$esf[[2]]
  
    ## calculate transformed derivatives, select relevant derivatives with ptot, drop unindentified parameters.
    gamma1 <- matrix(0, nrow = nrow(gamma1_pcm), ncol = m + o)
    for (j in mv) gamma1[, j] <- gamma1_pcm[, beta_index == j, drop = FALSE] %*% ov
    for (k in ov) gamma1[, k + m] <- rowSums(gamma1_pcm[, tau_index == k, drop = FALSE])
    gamma1 <- apply(gamma1, 2, "[", ptot)

    ## finally: the gradient
    agrad <- matrix(0, nrow = n, ncol = m + o)
    agrad[, mv] <- weights * (- dat + (gamma1 / gamma0)[, mv, drop = FALSE])
    agrad[, m + ov] <- weights * (- ctot + (gamma1 / gamma0)[, m + ov, drop = FALSE])

  } else {

    ## return value
    agrad <- matrix(0, nrow = n, ncol = m + o)

    ## observed NA patterns
    na_patterns <- factor(apply(is.na(dat), 1, function(z) paste(which(z), collapse = "\r")))

    ## loop through na patterns, select derivatives and calculate gradient
    for(i in seq_len(nlevels(na_patterns))) {

      ## parse NA patterns
      lev_i <- levels(na_patterns)[i]
      na_i <- which(na_patterns == lev_i)
      mv_i <- as.integer(strsplit(lev_i, "\r")[[1]])
      mv_i <- if(length(mv_i) < 1) mv else mv[-mv_i]
      m_i <- length(mv_i)
      mv_i_1 <- 1:m_i

      ## calculate/fetch necessary stuff for gradient
      weights_i <- weights[na_i]
      ptot_i <- ptot[na_i]
      gamma0 <- x$esf[[i]][[1]][ptot_i]
      gamma1_pcm <- x$esf[[i]][[2]]
      beta_index_i <- rep.int(mv_i, rep.int(o, m_i))
      tau_index_i <- rep.int(ov, m_i)

      ## calculate transformed derivatives, select relevant derivatives with ptot, drop unindentified parameters.
      gamma1 <- matrix(0, nrow = nrow(gamma1_pcm), ncol = m_i + o)
      for (j in mv_i_1) gamma1[, j] <- gamma1_pcm[, beta_index_i == mv_i[j], drop = FALSE] %*% ov
      for (k in ov) gamma1[, k + m_i] <- rowSums(gamma1_pcm[, tau_index_i == k, drop = FALSE])
      gamma1 <- apply(gamma1, 2, "[", ptot_i)
      if (!is.matrix(gamma1)) gamma1 <- matrix(gamma1, nrow = 1)

      ## finally: the gradient for NA group i
      agrad[na_i, mv_i] <- weights_i * (- dat[na_i, mv_i, drop = FALSE] + (gamma1 / gamma0)[, mv_i_1, drop = FALSE])
      agrad[na_i, m_i + ov] <- weights_i * (- ctot[na_i, , drop = FALSE] + (gamma1 / gamma0)[, m_i + ov, drop = FALSE])
    }

  }

  ## collect and return matrix of initial size with gradients plugged in.
  grad <- matrix(0, ncol = npar, nrow = length(weights_org))
  grad[weights_org > 0, ] <- agrad[, -c(1, m + 1)]
  return(grad)
}


### misc. internal functions

## prsm: calculate response probabilities for given thetas, betas and taus under the RSM.
prsm <- function(theta = NULL, beta = NULL, tau = NULL)
{
  ## check input
  stopifnot(!is.null(theta) && !is.null(beta) && !is.null(tau))

  ## if list input, recurse...
  if (is.list(theta)) return(lapply(theta, prsm, beta = beta, tau = tau))
  if (length(beta) > 1) return(lapply(as.list(beta), prsm, theta = theta, tau = tau))

  ## calculate probabilities
  num <- cbind(0, outer(theta - beta, 1:length(tau), "*"))         # k * (theta_i - beta_j ) for all i, j and k = 1:o, add 0 for cat 0
  num <- exp(t(apply(num, 1, function (x) x - cumsum(c(0, tau))))) # - sum(tau_k)
  denom <- rowSums(num)                                            # denominator: sum over all numerators
  return(num/denom)
}

## rpcm: calculate response matrices for given thetas, betas and taus under the RSM.
rrsm <- function(theta = NULL, beta = NULL, tau = NULL, nullcats = FALSE, return_setting = TRUE)
{
  ## check input
  stopifnot(!is.null(theta) && !is.null(beta) && !is.null(tau))
  if (is.list(beta)) stopifnot(is.list(tau) && (length(beta) == length(tau)))

  ## if list input, recurse...
  if (is.list(theta)) return(lapply(theta, rrsm, beta = beta, tau = tau))
  if (is.list(beta)) return(mapply(rrsm, beta, tau, MoreArgs =
                                   list(theta = theta, nullcats = nullcats, return_setting = return_setting), SIMPLIFY = FALSE))
  
  ## calculate response probabilities
  probs <- prsm(theta = theta, beta = beta, tau = tau)
  if(!is.list(probs)) probs <- list(probs)

  ## calculate response matrices for given set of item
  rsp_item_i <- function (probmat_i) {  #inline function to calculate responses per item
    oj <- ncol(probmat_i)               #calculate rsp vector once
    rsp <- apply(probmat_i, 1, function (p) sample.int(n = oj, size = 1, prob = p))
    if(!nullcats) {                     #recalculate if null categories not allowed
      while (length(unique.default(rsp)) != oj) {
        rsp <- apply(probmat_i, 1, function (p) sample.int(n = oj, size = 1, prob = p))
      }
    }
    return(rsp-1)
  }
  res <- lapply(probs, rsp_item_i)
  res <- do.call(cbind, res)
  
  if (return_setting)
    return(list(theta = theta, beta = beta, tau = tau, data = res))
  else
    return(res)
}
