## Function to fit a partial credit model, see e.g. Masters & Wright (1984)
pcmodel <- function (y, weights = NULL, nullcats = c("keep", "downcode", "ignore"), start = NULL,
  reltol = 1e-10, deriv = c("sum", "diff"), hessian = TRUE, maxit = 100L, full = TRUE, ...)
{
  ## argument matching
  deriv <- match.arg(deriv)

  ## process input and data
  nullcats <- match.arg(nullcats)
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
  if(any(mincat > 0)) {
    warning("Minimum score is not zero for all items (", paste("I", which(mincat > 0), sep = "", collapse = ", "), "). These items are scaled to zero.")
    y[, mincat > 0] <- scale.default(y[, mincat > 0], center = mincat[mincat > 0], scale = FALSE)
  }
  oj <- apply(y, 2, max, na.rm = TRUE)
  oj_vec <- lapply(oj, seq, from = 0)

  ## check for null categories and warn...
  nl_cats <- lapply(1:m, function (j) !(oj_vec[[j]] %in% y[, j]))
  nl_cats_items <- which(sapply(nl_cats, any))
  any_nl_cat_item <- length(nl_cats_items) > 0
  if (any_nl_cat_item) {
    warning("There are items with null categories (", paste("I", nl_cats_items, sep = "", collapse = ", "), ").")
    ## .. then treat according to nullcats (ignore = do nothing)
    if (nullcats == "downcode") {
      oj <- oj - sapply(nl_cats, sum)
      oj_vec <- lapply(1:m, function (j) oj_vec[[j]][1:(oj[j]+1)])
      for (j in nl_cats_items) {
        missing_scores <- sort(which(nl_cats[[j]]), decreasing = TRUE) - 1
        for (score in missing_scores) y[!is.na(y[, j]) & y[, j] > score, j] <-  y[!is.na(y[, j]) & y[, j] > score, j] - 1
      }
    }
    if (nullcats == "keep") {
      all_par <- rep.int(NA, sum(oj))
      est_par <- !unlist(lapply(nl_cats, "[", -1))
      all_par[1] <- 0
      est_par[1] <- FALSE
      oj <- oj - sapply(nl_cats, sum)
      oj_vec <- mapply(function(x, y) x[!y], oj_vec, nl_cats, SIMPLIFY = FALSE)
    }
  }
  oj_max <- sapply(oj_vec, max)         # maximum category number (different from number of categories if nullcats == "keep")

  ## remove unidentified items (only one category used)
  unid <- oj == 0
  nunid <- sum(unid)
  if (nunid) {
    y <- y[, !unid, drop = FALSE]
    oj <- oj[!unid]
    oj_vec <- oj_vec[!unid]
    oj_max <- oj_max[!unid]
    m <- m - nunid
    ident_items[unid] <- FALSE
    warning("There were unidentified items (only one category used, ", paste("I", which(unid), sep = "", collapse = ", "), ").")
  }
  mv <- 1:m

  ## check for missings
  y_na <- is.na(y)
  any_na <- any(y_na)

  ## calculate category totals and set starting values from these (calculation accounts for NAs)
  ## (generalisation of first proposal of Fischer & Molenaar, p. 50. If oj = 1: equal to starting values in raschmodel())
  ctot <- vector("list", length = m)
  for (j in seq_len(m)) ctot[[j]] <- as.vector(tapply(weights, factor(y[, j], levels = oj_vec[[j]]), sum))
  if (is.null(start)) {
    start <- lapply(ctot, function (x) - cumsum(diff.default(log(x)))) # delta_jk = log(x_j(k-1) - log(x_jk), beta_jk = cumsum_k=1^k(delta_jk)
    start <- unlist(start)
    start <- start[-1] - start[1]
    start[is.na(start)] <- 0
  }

  ## calculate person totals, catch missing person scores
  rs <- rowSums(y, na.rm = TRUE)
  ptot <- as.vector(tapply(weights, factor(rs, levels = 0:sum(oj_max)), sum))
  ptot[is.na(ptot)] <- 0

  ## unlist ctot, remove first category total (since beta_j0 = 0 for all 0), catch null categories
  ctot <- unlist(lapply(ctot, "[", -1))
  ctot[is.na(ctot)] <- 0
    
  ## build log-likelihood
  if(!any_na) {
    
    ## objective function: conditional log-likelihood
    cloglik <- function (par) {

      ## include epsilons for esf when there are nullcats & nullcats == "keep" (see Wilson, 1993)
      if (any_nl_cat_item && nullcats == "keep") {
        all_par[est_par] <- par
        esf_par <- all_par
      } else {
        esf_par <- c(0, par)
      }
      esf_par <- split(esf_par, rep.int(mv, oj_max))

      ## finally: the cll
      cll <- - sum(ctot * c(0, par)) - sum(ptot * log(elementary_symmetric_functions(par = esf_par, order = 0, diff = diff)[[1]]))

      ## catch degenerated cases (typically caused by non-finite gamma)
      if (is.na(cll) | !is.finite(cll)) cll <- -.Machine$double.xmax
      
      return(-cll)
    }

    ## static helper variables for analytical gradient
    parindex <- unlist(lapply(oj_vec, "[", -1))
    itemindex <- rep.int(mv, oj)
    xmat <- matrix(FALSE, nrow = n, ncol = sum(oj))
    for (i in 1:n) xmat[i, ] <- y[i, itemindex] == parindex 

    ## analytical gradient
    agrad <- function (par) {

        ## elementary symmetric functions
        if (any_nl_cat_item && nullcats == "keep") {
            all_par[est_par] <- par
            esf_par <- all_par
        } else {
            esf_par <- c(0, par)
        }
        esf <- elementary_symmetric_functions(par = split(esf_par, rep.int(mv, oj_max)), order = 1, diff = diff)

        gamma0 <- esf[[1]][rs + 1]
        gamma1 <- apply(esf[[2]], 2, "[", rs + 1)
        if (any_nl_cat_item && nullcats == "keep") gamma1 <- gamma1[, c(TRUE, est_par[-1]), drop = FALSE]

        ## calculate and return aggregated gradient
        return(- colSums((weights * (- xmat + (gamma1 / gamma0)))[, -1, drop = FALSE]))
    }

  } else {
    
    ## fetch different NA-patterns like Achim does, setup position vector of parameters
    na_patterns <- factor(apply(y_na, 1, function(z) paste(which(z), collapse = "\r")))
    lev_na_patterns <- levels(na_patterns)
    par_pos <- rep.int(mv, oj_max)

    ## setup na_pattern lists for various elements of loglik
    m_i <- mv_i <- rs_i <- ptot_i <- ctot_i <- par_i <- oj_max_i <- vector("list", length(lev_na_patterns))
    na_obs_i <- weights_i <- est_par_i <- xmat_i <- vector("list", length(lev_na_patterns))

    ## loop over observed NA patterns, calculate constant things once
    for(i in seq_along(lev_na_patterns)) {

      ## from pattern i: get NA item(s), observations of and weights with this pattern
      na_items_i <- as.integer(strsplit(lev_na_patterns[i], "\r")[[1]])
      n_na_items_i <- length(na_items_i)
      na_obs_i[[i]] <- which(na_patterns == lev_na_patterns[i])
      weights_i[[i]] <- weights[na_obs_i[[i]]]
        
      ## select subset
      if(n_na_items_i < 1) {            # no missings
        y_i <- y[na_obs_i[[i]], , drop = FALSE]
        par_i[[i]] <- rep.int(TRUE, sum(oj_max))
        m_i[[i]] <- m
        mv_i[[i]] <- mv
        oj_vec_i <- oj_vec
        oj_max_i[[i]] <- oj_max
      } else {                          # specific NA pattern
        y_i <- y[na_obs_i[[i]], -na_items_i, drop = FALSE]
        par_i[[i]] <- !(par_pos %in% na_items_i)
        m_i[[i]] <- m - n_na_items_i
        mv_i[[i]] <- mv[-na_items_i]
        oj_vec_i <- oj_vec[-na_items_i]
        oj_max_i[[i]] <- oj_max[-na_items_i]
      }

      ## calculate category totals and person totals for NA-group i
      ctot_i[[i]] <- vector("list", length = m_i[[i]])
      for (j in seq_len(m_i[[i]])) ctot_i[[i]][[j]] <- as.vector(tapply(weights_i[[i]], factor(y_i[, j], levels = oj_vec_i[[j]]), sum))
      ctot_i[[i]] <- unlist(lapply(ctot_i[[i]], "[", -1))
      ctot_i[[i]][is.na(ctot_i[[i]])] <- 0
      rs_i[[i]] <- rowSums(y_i)
      ptot_i[[i]] <- as.vector(tapply(weights_i[[i]], factor(rs_i[[i]], levels = 0:sum(oj_max_i[[i]])), sum))
      ptot_i[[i]][is.na(ptot_i[[i]])] <- 0

      ## calculate helper variables for gradient
      parindex_i <- unlist(lapply(oj_vec_i, "[", -1))
      itemindex_i <- unlist(rep.int(mv_i[[i]], oj[mv_i[[i]]]))
      est_par_i[[i]] <- !unlist(lapply(nl_cats, "[", -1)[mv_i[[i]]])
      xmat_i[[i]] <- matrix(FALSE, nrow = length(na_obs_i[[i]]), ncol = sum(oj[mv_i[[i]]]))
      for (j in 1:length(na_obs_i[[i]])) xmat_i[[i]][j, ] <- (y[na_obs_i[[i]], , drop = FALSE])[j, itemindex_i] == parindex_i
    }

    ## fun for mapply (calculates cll contributions per na_group with variable par's)
    cll_i <- function (ctot_i, par_i, ptot_i, esf_par_i) {
      - sum(ctot_i * par_i) - sum(ptot_i * log(elementary_symmetric_functions(par = esf_par_i, order = 0, diff = diff)[[1]]))
    }
    
    ## objective function: conditional log-likelihood (build incrementally by looping over the NA-patterns
    cloglik <- function (par) {

      ## add first zero par, then select current parameters
      ## include epsilons for esf when there are nullcats & nullcats == "keep" (see Wilson, 1993)
      if (any_nl_cat_item && nullcats == "keep") {
        all_par[est_par] <- par
        parx <- all_par
      } else {
        parx <- c(0, par)
      }
      esf_par <- lapply(par_i, function (x) parx[x])
      par <- lapply(esf_par, function (x) x[!is.na(x)])
      esf_par <- mapply(split, esf_par, mapply(function (x, y) rep(1:x, y), m_i, oj_max_i))

      ## conditional log-likelihood
      cll <- sum(mapply(cll_i, ctot_i, par, ptot_i, esf_par, SIMPLIFY = TRUE, USE.NAMES = FALSE))

      ## catch degenerated cases (typically cause by non-finite gamma)
      if(is.na(cll) | !is.finite(cll)) cll <- -.Machine$double.xmax

      return(-cll)
    }

    ## static helper variables for analytical gradient
    grad <- matrix(0, nrow = n, ncol = sum(oj))
    itemindex <- rep.int(mv, oj)

    ## analytical gradient
    agrad <- function (par) {
        
        ## elementary symmetric functions
        if (any_nl_cat_item && nullcats == "keep") {
            all_par[est_par] <- par
            esf_par <- all_par
        } else {
            esf_par <- c(0, par)
        }
        esf_par <- lapply(par_i, function (x) esf_par[x])
        esf_par <- mapply(split, esf_par, mapply(function (x, y) rep(1:x, y), m_i, oj_max_i))
        esf <- mapply(elementary_symmetric_functions, par = esf_par, MoreArgs = list(order = 1, diff = diff), SIMPLIFY = FALSE)

        ## loop over observed NA patterns and calculate gradients
        for(i in seq_along(lev_na_patterns)) {

            ## select gamma zero and first derivatives with rs_i
            gamma0_i <- esf[[i]][[1]][rs_i[[i]] + 1]
            gamma1_i <- apply(esf[[i]][[2]], 2, "[", rs_i[[i]] + 1)
            if (!is.matrix(gamma1_i)) gamma1_i <- matrix(gamma1_i, nrow = 1)

            ## finally: the gradient for NA group i
            if (any_nl_cat_item && nullcats == "keep") gamma1_i <- gamma1_i[, est_par_i[[i]], drop = FALSE]
            grad[na_obs_i[[i]], itemindex %in% mv_i[[i]]] <- weights_i[[i]] * (- xmat_i[[i]] + (gamma1_i/ gamma0_i))
        }

        return(- colSums(grad[, -1, drop = FALSE]))
    }
  }

  ## optimization
  opt <- optim(par = start, fn = cloglik, gr = agrad, method = "BFGS",
               hessian = hessian, control = list(reltol = reltol, maxit = maxit, ...))

  ## final estimates ...
  est <- opt$par
  names(est) <- if (is.null(colnames(y))) {
    paste("I", rep.int(which(ident_items), oj), "-C", unlist(lapply(oj, seq_len)), sep = "")[-1] 
  } else paste(rep(colnames(y)[which(ident_items)], oj), "-C", unlist(lapply(oj, seq_len)), sep = "")[-1]
  ## FIXME/Z: use parindex instead of unlist(lapply(oj, seq_len))
  ## at least if nullcats = "keep"?
  
  ## ... and (if requested) esf of these ...
  if (full) {
    if (any_nl_cat_item && nullcats == "keep") {
      all_par[est_par] <- est
      parx <- all_par
    } else {
      parx <- c(0, est)
    }
    esf <- if (any_na) {
      esf_par <- lapply(par_i, function (x) parx[x])
      esf_par <- mapply(split, esf_par, mapply(function (x, y) rep(1:x, y), m_i, oj_max_i))
      mapply(elementary_symmetric_functions, par = esf_par, MoreArgs = list(order = 1, diff = diff), SIMPLIFY = FALSE)
    } else {
      elementary_symmetric_functions(par = split(parx, rep.int(mv, oj_max)), order = 1, diff = diff)
    }

    ## ... as well as variance-covariance matrix
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
              categories = lapply(oj_vec, "[", -1),
              n = sum(weights_org > 0),
              n_org = n_org,
              weights = if (identical(as.vector(weights_org), rep.int(1L, n_org))) NULL else weights_org,
              na = any_na,
              nullcats = if (any_nl_cat_item && nullcats == "keep") lapply(nl_cats, "[", -1) else NULL,
              esf = esf,
              loglik = -opt$value,
              df = length(est),
              code = opt$convergence,
              iterations = tail(na.omit(opt$counts), 1L),
              reltol = reltol)
  class(res) <- "pcmodel"
  return(res)
}

print.pcmodel <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  cat("PCM item category parameters:\n")
  print(coef(x), digits = digits)
  invisible(x)
}

coef.pcmodel <- function (object, ...) object$coefficients

vcov.pcmodel <- function (object, ...) object$vcov

logLik.pcmodel <- function (object, ...) structure(object$loglik, df = object$df, class = "logLik")

weights.pcmodel <- function (object, ...) if (is.null(object$weights)) rep.int(1L, object$n_org) else object$weights

summary.pcmodel <- function (object, vcov. = NULL, ...) {
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
  class(object) <- "summary.pcmodel"
  return(object)
}

print.summary.pcmodel <- function (x, digits = max(3, getOption("digits") - 3),
                                   signif.stars = getOption("show.signif.stars"), ...)
{
  if (is.null(x$call)) {
    cat("\nPartial credit model\n\n")  
  } else {
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  }

  if (any(!x$items)) cat("Excluded items:",
                         paste(names(x$items)[!x$items], collapse = ", "), "\n\n")

  cat("Item category parameters:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)

  cat("\nLog-likelihood:", format(signif(x$loglik, digits)),
      "(df =", paste(x$df, ")", sep = ""), "\n")
  cat("Number of iterations in BFGS optimization:", x$iterations, "\n\n")
  invisible(x)
}

plot.pcmodel <-  function (x, type = c("regions", "profile", "curves", "information", "piplot"), ...)
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

predict.pcmodel <- function (object, newdata = NULL, type = c("probability",
  "cumprobability", "mode", "median", "mean", "category-information",
  "item-information", "test-information"), ref = NULL, ...)
{
  ## check type, process newdata, if NULL, use person parameters of given model object
  type <- match.arg(type)
  if (is.null(newdata)) {
    rs <- rowSums(object$data, na.rm = TRUE)
    rs <- rs[0 < rs & rs < sum(apply(object$data, 2, max))]
    newdata <- personpar(object, vcov = FALSE)[rs]
    names(newdata) <- NULL
  }
  nms <- names(newdata)

  ## threshold parameters, labels and raw probabilities
  tp <- threshpar(object, type = "mode", ref = ref, vcov = FALSE)
  lbs <- unlist(mapply(paste, names(tp), lapply(tp, function (j) c("C0", names(j))), sep = "-", SIMPLIFY = FALSE), use.names = FALSE)
  probs <- ppcm(theta = newdata, delta = tp)

  ## if requested: compute test/item/category information (see Muraki, 1993, for details and further references)
  if (grepl("information", type)) {
    m <- length(tp)
    clnms <- names(tp)
    if (type == "category-information") {
      info <- matrix(NA, nrow = length(newdata), ncol = length(lbs))
      colnames(info) <- lbs
    } else {
      info <- matrix(NA, nrow = length(newdata), ncol = m)
      colnames(info) <- clnms
    }
    for (j in 1:m) {
      ojm <- matrix(rep(1:ncol(probs[[j]]) - 1, each = length(newdata)), nrow = length(newdata), ncol = ncol(probs[[j]]))
      idx <- grepl(clnms[j], lbs)
      iteminfo <- rowSums((ojm - rowSums(ojm * probs[[j]]))^2 * probs[[j]])
      if (type == "category-information") {
        info[, idx] <- probs[[j]] * iteminfo
      } else {
        info[, j] <- iteminfo
      }
    }
  }

  ## return as requested in type, for RM mode, median, mean is the same
  switch(type,
         "probability" = {
           ## create (named) matrix with probabilities
           probs <- do.call("cbind", probs)
           dimnames(probs) <- list(nms, lbs)
           probs
           },
         "cumprobability" = {
           ## create (named) matrix with probabilities
           probs <- lapply(probs, function (probsj) t(apply(probsj, 1, function (x) rev(cumsum(rev(x))))))
           probs <- do.call("cbind", probs)
           dimnames(probs) <- list(nms, gsub("(.*)-(.*)", "\\1>=\\2", lbs))
           probs
           },
         "mode" = {
           ## create (named) matrix with probabilities
           probs <- sapply(probs, apply, 1, which.max) - 1
           dimnames(probs) <- list(nms, unique(gsub("(.*)-C[[:digit:]]+", "\\1", lbs)))
           probs
         },
         "median" = {
           ## create (named) matrix with probabilities
           probs <- sapply(probs, function (probsj) apply(probsj, 1, which.max) - 1)
           dimnames(probs) <- list(nms, unique(gsub("(.*)-C[[:digit:]]+", "\\1", lbs)))
           probs
         },
         "mean" = {
           ## create (named) matrix with probabilities
           probs <- round(sapply(probs, function (probsj) apply(probsj, 1, function (j) sum(j * 0:(length(j) - 1)))))
           dimnames(probs) <- list(nms, unique(gsub("(.*)-C[[:digit:]]+", "\\1", lbs)))
           probs
         },
         "category-information" = info,
         "item-information" = info,
         "test-information" = matrix(rowSums(info), ncol = 1))
}

itempar.pcmodel <- function (object, ref = NULL, alias = TRUE, vcov = TRUE, ...)
{
  ## extract estimated item category parameters and labels, include restricted parameter
  cf <- c(0.00, coef(object))
  b <- length(cf)
  m <- sum(object$items)
  lbs <- if (!is.null(colnames(object$data))) colnames(object$data) else paste0("I", 1:m)
  oj <- sapply(object$categories, length)
  ojc <- cumsum(oj)

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
    stopifnot(nrow(ref) == m & ncol(ref) == m)
  } else stop("Argument 'ref' is misspecified (see ?itempar for possible values).")

  ## if not given, specify contrast matrix
  if (is.matrix(ref)) {
    D <- ref
  } else {      
    D <- diag(m)
    D[, ref] <- D[, ref] - 1/length(ref)
  }

  ## transform estimated item category parameters (= sum of threshold parameters)
  ## to mean absolute threshold parameters ...
  C <- matrix(0, nrow = m, ncol = b)
  for (j in 1:m) C[j, sum(oj[1:j])] <- 1/oj[j]
  cf <- as.vector(C %*% cf)

  ## if requested: create adjusted vcov
  if (vcov) {
    vc <- C %*% rbind(0, cbind(0, vcov(object))) %*% t(C)
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
  rv <- structure(cf, class = "itempar", model = "PCM", ref = ref, alias = alias, vcov = vc)
  return(rv)
}

threshpar.pcmodel <- function (object, type = c("mode", "median", "mean"), ref = NULL,
                               alias = TRUE, relative = FALSE, cumulative = FALSE, vcov = TRUE, ...)
{
  ## check input
  type <- match.arg(type)

  ## extract relevant informations
  m <- sum(object$items)
  tp <- c(0.00, coef(object)) ## cumualtive absolute item threshold parameters
  b <- length(tp)
  lbs <- names(tp)
  oj <- sapply(object$categories, length)
  ojc <- cumsum(oj)
  lbs[1] <- if (lbs[2] == "I1-C2") "I1-C1" else if (lbs[2] == "I2-C1") "I1-C1" else paste0(colnames(object$data)[1], "-C1")
  ilbs <- unique(gsub("(.*)\\-(.*)", "\\1", lbs))
  clbs <- lapply(oj, function (oj) paste0("C", 1:oj))

  ## process argument relative
  if (relative) {
    
    if (type == "mode") {

      ## process ref
      if (!is.list(ref)) {
        if (is.null(ref)) {
          ref <- lapply(oj, seq)
        } else if (is.vector(ref) && is.character(ref)) {
          stopifnot(ref %in% unlist(clbs))
          ref <- lapply(clbs, function (j) which(j %in% ref))
        } else if (is.vector(ref) && is.numeric(ref)) {
          ref <- as.integer(ref)
          stopifnot(all(ref %in% 1:min(oj)))
          ref <- split(rep(ref, m), rep(1:m, length(ref)))
        } else if (is.matrix(ref) && is.numeric(ref)) {
          stopifnot(ncol(ref) == b && nrow(ref) == b)
          ref2 <- vector(mode = "list", length = m)
          for (i in 1:m) ref2[[i]] <- ref
          ref <- ref2
        } else stop("Argument 'ref' is misspecified (see ?threshpar for possible values).")
      } else {
        if (length(ref) < m) stop("Not enough restrictions provided in argument 'ref'.")
        else {
          for (j in 1:m) {
            if (is.null(ref[[j]])) {
              ref[[j]] <- 1:oj[j]
            } else if (is.vector(ref[[j]]) && is.character(ref[[j]])) {
              stopifnot(ref %in% clbs[[j]])
              ref[[j]] <- which(clbs[[j]] %in% ref[[j]])
            } else if (is.vector(ref[[j]]) && is.numeric(ref[[j]])) {
              ref[[j]] <- as.integer(ref[[j]])
              stopifnot(ref[[j]] %in% 1:oj[j])
            } else if (is.matrix(ref[[j]]) && is.numeric(ref[[j]])) {
              stopifnot(ncol(ref[[j]]) == oj[j] && nrow(ref[[j]]) == oj[j])
            } else stop("Argument 'ref' is misspecified (see ?threshpar for possible values).")
          }
        }
      }

      ## transform estimated cumulative absolute item threshold parameters
      ## into relative item threshold parameters
      C <- diag(b)
      for (j in 1:m) {
        for (k in (2:oj[j])) C[c(0, ojc)[j] + k, c(0, ojc)[j] + k - 1] <- -1
        C[c(0, ojc)[j] + 1:oj[j], ojc[j]] <- C[c(0, ojc)[j] + 1:oj[j], ojc[j]] - 1/oj[j]
      }
      tp <- as.vector(C %*% tp)
      ## if vcov requested: create adjusted vcov
      if (vcov) vc <- C %*% rbind(0, cbind(0, vcov(object))) %*% t(C)

      ## if not given, specify contrast matrix (block diagonal matrix with item specific contrasts)
      D <- diag(b)
      if (is.matrix(ref[[1]])) {
        for (j in 1:m) D[c(0, ojc)[j] + 1:oj[j], c(0, ojc)[j] + 1:oj[j]] <- ref[[j]]
      } else {
        for (j in 1:m) D[c(0, ojc)[j] + 1:oj[j], c(0, ojc)[j] + ref[[j]]] <- D[c(0, ojc)[j] + 1:oj[j], c(0, ojc)[j] + ref[[j]]] - 1/length(ref[[j]])
      }

      ## apply ref
      ## if vcov requested: adjust vcov too
      tp <- as.vector(D %*% tp)
      tp <- split(tp, rep.int(1:m, oj))
      if (vcov) {
        vc <- D %*% vc %*% t(D)
      } else { # else return NA matrix
        vc <- matrix(NA, nrow = b, ncol = b)
      }

      ## if cumulative relative item threshold parameters are requested: transform 
      ## relative item threshold parameters to cumulative relative item threshold parameters
      if (cumulative) {
        tp <- lapply(tp, cumsum)
        C <- matrix(0, nrow = b, ncol = b)
        for (j in 1:m) {
          for (k in 1:oj[j]) C[c(0, ojc)[j] + k, c(0, ojc)[j] + 1:k] <- 1
        }
        if (vcov) vc <- C %*% vc %*% t(C)
      }

    } else stop("Relative threshold parameters not implemented for types other than mode.")

  } else {

    if (type == "mode") {

      ## process ref
      if (is.null(ref)) {
        ref <- 1:b
      } else if (is.vector(ref) && is.character(ref)) {
        stopifnot(all(ref %in% lbs))
        ref <- which(lbs %in% ref)
      } else if (is.vector(ref) && is.numeric(ref)) {
        ref <- as.integer(ref)
        stopifnot(all(ref %in% 1:b))
      } else if (is.matrix(ref) && is.numeric(ref)) {
        stopifnot(nrow(ref) == b && ncol(ref) == b)
      } else stop("Argument 'ref' can only be a character vector with threshold parameter labels or a numeric vector with threshold parameter indices.")

      ## transform estimated cumulative absolute item threshold
      ## parameters to absolute item threshold parameters
      C <- diag(b)
      for (j in 1:m) {
        for (k in (2:oj[j])) C[c(0, ojc)[j] + k, c(0, ojc)[j] + k - 1] <- -1
      }
      tp <- as.vector(C %*% tp)
      if (vcov) vc <- C %*% rbind(0, cbind(0, vcov(object))) %*% t(C)

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
      tp <- split(tp, rep.int(1:m, oj))
      if (vcov) {
        vc <- D %*% vc %*% t(D)
      } else { # else return NA matrix
        vc <- matrix(NA, nrow = b, ncol = b)
      }

      ## if cumulative absolute item threshold paramters are requested
      ## cumulate calculated absolute item threshold parameters and adjust vcov if requested
      if (cumulative) {
        tp <- lapply(tp, cumsum)
        if (vcov) {
          C <- matrix(0, nrow = b, ncol = b)
          for (j in 1:m) {
            for (k in 1:oj[j]) C[c(0, ojc)[j] + k, c(0, ojc)[j] + 1:k] <- 1
          }
          vc <- C %*% vc %*% t(C)
        }
      }

    } else {

      ## process ref, create tp list, setup NA vcov
      if (!is.null(ref)) warning("Argument 'ref' is not processed for types other than mode.")
      tp <- split(tp, rep.int(1:m, oj))
      vc <- matrix(NA, nrow = b, ncol = b)

      if (type == "median") {

        ## function to find locations on theta axis
        zmedian <- function (theta = NULL, delta = NULL, geq = NULL, ncat = NULL) {
          rowSums(ppcm(theta = theta, delta = delta)[, (geq + 1):ncat, drop = FALSE]) - 0.5
        }
        
        ## loop though items and find locations by means of zmedian() and uniroot()
        for (j in seq_along(tp)) tp[[j]] <- sapply(1:oj[j], function (geq) uniroot(f = zmedian, interval = c(-10, 10), delta = tp[[j]], geq = geq, ncat = oj[j] + 1)$root)
      }

      if (type == "mean") {
        ## function to find locations on theta axis
        xpct <- lapply(oj, function (oj) 1:oj - 0.5)
        zexpct <- function (theta = NULL, delta = NULL, expct = NULL) ppcm(theta = theta, delta = delta) %*% 0:length(delta) - expct
        
        ## loop though items and find locations by means of zexpct() and uniroot()
        for (j in seq_along(tp)) tp[[j]] <- sapply(xpct[[j]], function (xp) uniroot(f = zexpct, interval = c(-10, 10), delta = tp[[j]], expct = xp)$root)
      }

      ## if cumulative parameter are requsted, just cumulate ...
      if (cumulative) tp <- lapply(tp, cumsum)
    }
  }

  ## set labels
  names(tp) <- ilbs
  rownames(vc) <- colnames(vc) <- lbs
  for (i in 1:m) names(tp[[i]]) <- paste0("C", 1:oj[i])

  ## process argument alias
  if (!alias && type == "mode") {
    if (is.matrix(ref)) {
      ## FIXME: Implement alias when ref is a specific constrast matrices -> detect linear dependent columns?
      stop("Processing of argument 'alias' not implemented with a contrast matrix given in argument 'ref'.")
    } else {
      if (relative) {
        i <- split(1:b, rep(1:m, oj))
        alias <- vector(mode = "list", length = m)
        for (j in 1:m) {
          tp[[j]] <- tp[[j]][-ref[[j]][1]]
          i[[j]] <- i[[j]][-ref[[j]][1]]
          alias[[j]] <- ref[[j]][1]
        }
        i <- unlist(i)
        vc <- vc[i, i]
        names(alias) <- ilbs
      } else {
        ref1 <- ref[1]
        i <- split(1:b, rep(1:m, oj))
        item <- which(sapply(i, function (j) ref1 %in% j))
        tp[[item]] <- tp[[item]][-which(ref1 == i[[item]])]
        vc <- vc[-ref1, -ref1]
        alias <- paste0("I", item, "-C", which(ref1 == i[[item]]))
        names(alias) <- ilbs[item]
      }
    }
  }

  ## setup and return result object
  rv <- structure(tp, class = "threshpar", model = "PCM", type = type, ref = ref, relative = relative, alias = alias, vcov = vc)
  return(rv)
}

discrpar.pcmodel <- function (object, ref = NULL, alias = TRUE, vcov = TRUE, ...)
{
  ## check input
  if (!is.null(ref)) warning("Argument 'ref' is currently not processed.")  ## FIXME: Implement argument ref for discrimination parameters

  ## extract labels and number of items
  lbs <- c("", names(coef(object)))
  lbs[1] <- if (lbs[2] == "I1-C2") "I1-C1" else if (lbs[2] == "I2-C1") "I1-C1" else paste0(colnames(object$data)[1], "-C1")
  lbs <- unique(sapply(strsplit(lbs, "-"), "[[", 1))
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
  rv <- structure(dp, .Names = if (is.logical(alias)) lbs, class = "discrpar", model = "PCM", ref = ref, alias = alias, vcov = vc)
  return(rv)
}

personpar.pcmodel <- function (object, ref = NULL, vcov = TRUE, interval = NULL, tol = 1e-8, ...)
{
  ## extract item parameters and relevant informations
  tp <- threshpar(object, type = "mode", ref = ref, vcov = FALSE)
  tp <- lapply(tp, cumsum)
  m <- length(tp)
  oj <- sapply(tp, length)
  ojvl <- lapply(oj, seq)
  ojv <- unlist(ojvl)
  rng <- 1:(sum(oj) - 1)

  ## iterate over raw scores (from 1 until rmax-1)
  if(is.null(interval)) interval <- c(-1, 1) * qlogis(1/m * 1e-3) #FIXME: 1e3 enough?
  pp <- sapply(rng, function(rawscore) {
    uniroot(function(pp) rawscore - sum(ojv * exp(ojv * pp - unlist(tp)) / unlist(sapply(1:m, function (j) rep.int(1 + sum(exp(ojvl[[j]] * pp - tp[[j]])), oj[j])))),
      interval = interval, tol = tol)$root
  })

  if (vcov) {
    ## if vcov is requested, fetch relevant data for loglik
    y <- object$data[weights(object) > 0, , drop = FALSE]
    rs <- rowSums(y)
    rf <- tabulate(rs, nbins = sum(oj) - 1)
    cs <- sapply(1:m, function (j) tabulate(y[, j], nbins = oj[j]))

    ## remove unidentified parameters
    rs <- rs[rf != 0]
    rng <- rng[rf != 0]
    pp <- pp[rf != 0]
    rf <- rf[rf != 0]

    ## obtain Hessian from objective function: joint log likelihood
    cloglik <- function (pp) {
      ppx <- lapply(ojvl, function (x) outer(x, pp)) # l * theta_i
      - sum(rf * rng * pp) + sum(unlist(mapply("*", cs, tp))) + sum(rf * rowSums(log(1 + mapply(function (x, y) colSums(exp(x - y)), ppx, tp))))
    }
    vc <- solve(optim(pp, fn = cloglik, hessian = TRUE, method = "BFGS", control = list(reltol = tol, maxit = 0, ...))$hessian)    
  } else {      
    vc <- matrix(NA_real_, nrow = length(rng), ncol = length(rng))
  }
  colnames(vc) <- rownames(vc) <- rng

  ## setup and return result object
  rv <- structure(pp, .Names = rng, class = "personpar", model = "PCM", vcov = vc)
  return(rv)
}

nobs.pcmodel <- function (object, ...)
{
  return(object$n)
}

bread.pcmodel <- function(x, ...) x$vcov * x$n

estfun.pcmodel <- function (x, ...) {
  ## get relevant informations
  dat <- x$data                    # completely cleaned (downcoded, null cats treatment, weights) data.
  weights_org <- weights(x)
  weights <- weights_org[weights_org > 0]
  n <- nrow(dat)
  m <- ncol(dat)
  oj_vec <- x$categories
  oj <- sapply(x$categories, length)
  npar_all <- sum(oj)
  npar_ident <- x$df
  ptot <- rowSums(dat, na.rm = TRUE) + 1 # +1 because gamma of score 0 is in row 1.

  ## helper variables
  parindex <- unlist(oj_vec)
  itemindex <- rep.int(1:m, oj)

  ## calculate gradient
  if (!x$na) {

    ## select gamma zero and first derivatives with ptot
    gamma0 <- x$esf[[1]][ptot]
    gamma1 <- apply(x$esf[[2]], 2, "[", ptot)

    ## construct data matrix ('selection' matrix, 0/1, cols = parameters)
    if (!is.null(x$nullcats)) { ## null cats & strategy == 'keep', remove column of unidentified par
      est_par <- !unlist(x$nullcats)
      gamma1 <- gamma1[, est_par, drop = FALSE]
    }
    xmat <- matrix(FALSE, nrow = n, ncol = npar_all)
    for (i in 1:n) xmat[i, ] <- dat[i, itemindex] == parindex 

    ## calculate gradient
    agrad <- weights * (- xmat + (gamma1 / gamma0))

  } else {

    ## return value & helper variables
    agrad <- matrix(0, nrow = n, ncol = npar_all)
    mv <- 1:m

    ## observed NA patterns 
    na_patterns <- factor(apply(is.na(dat), 1, function(z) paste(which(z), collapse = "\r")))

    ## loop through na patterns, select derivatives and calculate gradient
    for(i in seq_len(nlevels(na_patterns))) {

      ## parse NA patterns and setup necessary stuff for gradient calculation
      lev_i <- levels(na_patterns)[i]
      na_i <- which(na_patterns == lev_i)
      n_na_i <- length(na_i)
      mv_i <- as.integer(strsplit(lev_i, "\r")[[1]])
      mv_i <- if(length(mv_i) < 1) mv else mv[-mv_i]
      oj_i <- oj[mv_i]
      oj_vec_i <- oj_vec[mv_i]
      weights_i <- weights[na_i]
      ptot_i <- ptot[na_i]
      dat_i <- dat[na_i, , drop = FALSE]

      ## select gamma zero and first derivatives with ptot_i
      gamma0_i <- x$esf[[i]][[1]][ptot_i]
      gamma1_i <- apply(x$esf[[i]][[2]], 2, "[", ptot_i)
      if (!is.matrix(gamma1_i)) gamma1_i <- matrix(gamma1_i, nrow = 1)

      ## construct data matrix ('selection' matrix, 0/1, cols = parameters) for NA group i
      parindex_i <- unlist(oj_vec_i)
      itemindex_i <- rep.int(mv_i, oj_i)

      if (!is.null(x$nullcats)) { ## null cats & strategy == 'keep', remove column of unidentified par
        est_par_i<- !unlist(x$nullcat[mv_i])
        gamma1_i <- gamma1_i[, est_par_i, drop = FALSE]
      }

      xmat_i <- matrix(FALSE, nrow = n_na_i, ncol = sum(oj_i))
      for (i in 1:n_na_i) xmat_i[i, ] <- dat_i[i, itemindex_i] == parindex_i

      ## finally: the gradient for NA group i
      agrad[na_i, itemindex %in% mv_i] <- weights_i * (- xmat_i + (gamma1_i/ gamma0_i))
    }

  }  

  ## collect and return matrix of initial size with gradients plugged in.
  grad <- matrix(0, ncol = npar_ident, nrow = length(weights_org))
  grad[weights_org > 0, ] <- agrad[, -1, drop = FALSE]
  return(grad)
}


### misc. internal functions

## ppcm: calculate response probabilities for given thetas and deltas under the PCM.
ppcm <- function(theta = NULL, delta = NULL)
{
  ## check input
  stopifnot(!is.null(theta) && !is.null(delta))
  
  ## if list input, recurse...
  if (is.list(theta)) return(lapply(theta, ppcm, delta = delta))
  if (is.list(delta)) return(lapply(delta, ppcm, theta = theta))

  ## calculate probabilities
  num <- cbind(0, outer(theta, delta, "-")) # all possible differences, 0 for category zero (\sum_0^0 \def 0)
  num <- t(exp(apply(num, 1, cumsum)))      # numerator: all possible cumulative sums
  denom <- rowSums(num)                     # denominator: sum over cumulative sums
  return(num/denom)
}

## rpcm: calculate response matrices for given thetas and deltas under the PCM.
rpcm <- function(theta = NULL, delta = NULL, nullcats = FALSE, return_setting = TRUE)
{
  ## check input
  stopifnot(!is.null(theta) && !is.null(delta))
  
  ## if list input, recurse... (for deltas: one list means several items, two means several groups of items)
  if (is.list(theta)) return(lapply(theta, rpcm, delta = delta, nullcats = nullcats, return_setting = return_setting))
  if (is.list(delta) && is.list(delta[[1]])) return(lapply(delta, rpcm, theta = theta, nullcats = nullcats, return_setting = return_setting))

  ## calculate response probabilities
  probs <- ppcm(theta = theta, delta = delta)
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
    return(list(delta = delta, theta = theta, data = res))
  else
    return(res)
}
