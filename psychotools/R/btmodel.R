## workhorse function
btmodel <- function(y, weights = NULL, type = c("loglin", "logit"), ref = NULL,
  undecided = NULL, position = NULL, start = NULL, vcov = TRUE, estfun = FALSE, ...)
{
  ## main arguments
  if(missing(y)) stop("response missing")
  stopifnot(inherits(y, "paircomp"))  
  addargs <- list(...)
  if ("estfun" %in% names(addargs)) warning("The argument 'estfun' is deprecated and not longer used.
    The individual score contributions can be computed with the model-specific estfun() method.")

  ## basic paircomp properties
  lab <- labels(y)
  nsubj <- length(y)
  nobj <- length(lab)
  npc <- nobj * (nobj - 1)/2
  mscale <- mscale(y)
  if(max(abs(mscale)) > 1L) stop("comparisons on likert scales not yet implemented")
  has_ties <- mscale[2L] == 0L
  ix <- which(upper.tri(diag(nobj)), arr.ind = TRUE)
  
  ## weights
  if(is.null(weights)) weights <- 1L
  weights <- rep(weights, length.out = nsubj)
  nsubj <- sum(weights > 0L)

  ## further arguments
  type <- match.arg(type, c("loglin", "logit"))
  if(is.null(ref)) ref <- nobj
  if(is.character(ref)) ref <- match(ref, lab)
  if(is.null(undecided)) undecided <- has_ties
  if(undecided & type != "loglin") stop("only log-linear model can handle ties")  
  if(!has_ties & undecided) stop("data have no ties") ## FIXME: set coef to NA and logprob to -Inf
  npar <- nobj - !undecided
  if(is.null(position)) position <- attr(y, "ordered")
  if(position) stop("position effects not yet implemented")
  
  ## basic aggregation quantities
  ytab <- summary(y, weights = weights)
  if(has_ties) ytab <- ytab[, c(1L, 3L, 2L), drop = FALSE]
  if(!undecided) ytab <- ytab[, 1L:2L, drop = FALSE]
    
  ## set up auxiliary model
  if(type == "loglin") {
    famaux <- poisson()
    yaux <- as.vector(t(ytab))
    xaux <- matrix(0, nrow = 3L * npc, ncol = nobj + npc)
    for(i in 1L:nrow(ix)) {
      ## xaux[i*3L - (2:1), ix[i, 1L:2L]] <- c(1, -1, -1, 1) ## DHK parametrization
      xaux[i*3L - (2L:0L), ix[i,1L:2L]] <- c(1, 0, 0.5, 0, 1, 0.5)
      xaux[i*3L - (2L:0L), nobj + i] <- 1
    }
    xaux[,1L:nobj] <- xaux[,c((1L:nobj)[-ref], ref)]    
    xaux[,nobj] <- rep(c(0, 0, 1), npc)
    if(!undecided) {
      xaux <- xaux[,-nobj]
      xaux <- xaux[-((1L:npc) * 3L),]
    }
  } else {
    famaux <- binomial(link = type)
    yaux <- ytab
    xaux <- matrix(0, nrow = npc, ncol = nobj)
    for(i in 1L:nrow(ix)) xaux[i, ix[i, 1L:2L]] <- c(1, -1)
    xaux <- xaux[,-ref]    
  }

  ## fit auxiliary model and extract information
  fm <- suppressWarnings(glm.fit(xaux, yaux, family = famaux, control = glm.control(...), start = start))
  par <- fm$coefficients[1L:npar]
  nam <- c(lab[-ref], if(undecided) "(undecided)" else NULL)
  names(par) <- nam
  if(vcov) {
    vc <- summary.glm(fm, correlation = FALSE)$cov.unscaled[1L:npar, 1L:npar, drop = FALSE]
    rownames(vc) <- colnames(vc) <- nam
  } else {
    vc <- matrix(NA_real_, nrow = npar, ncol = npar, dimnames = list(nam, nam))
  }

  ## log-probabilities and log-likelihood
  par2logprob <- switch(type,
    "loglin" = function(i) {
      p <- rep(0, nobj)
      p[-ref] <- if(undecided) par[-npar] else par
      p <- p[ix[i,]]
      if(undecided) p <- c(p, par[npar] + mean(p))
      p - log(sum(exp(p)))
    },
    "logit" = function(i) {
      p <- rep(0, nobj)
      p[-ref] <- par
      plogis(c(-1, 1) * diff(p[ix[i,]]), log.p = TRUE)
    }
  )
  logp <- t(sapply(1L:npc, par2logprob))
  loglik <- sum(logp * ytab)

  ## collect, class, and return
  rval <- list(
    y = y,
    coefficients = par,
    vcov = vc,
    loglik = loglik,
    df = npar,
    weights = if(isTRUE(all.equal(weights, rep(1, nsubj)))) NULL else weights,
    n = nsubj,
    nobs = nsubj,
    type = type,
    ref = lab[ref],
    undecided = undecided,
    position = position,
    labels = lab
  )
  if(estfun) rval$estfun <- estfun.btmodel(rval)
  class(rval) <- "btmodel"
  return(rval)
}

## standard methods
vcov.btmodel <- function(object, ...) object$vcov

logLik.btmodel <- function(object, ...) structure(object$loglik, df = object$df, class = "logLik")

deviance.btmodel <- function(object, ...) -2 * object$loglik

## more elaborate methods
print.btmodel <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("BT regression coefficients:\n")
  print(coef(x), digits = digits)
  invisible(x)
}

summary.btmodel <- function(object, vcov. = NULL, ...)
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
  class(object) <- "summary.btmodel"
  return(object)
}

print.summary.btmodel <- function(x, digits = max(3, getOption("digits") - 3), 
    signif.stars = getOption("show.signif.stars"), ...)
{
  if(is.null(x$call)) {
    cat("\nBradley-Terry regression model\n\n")  
  } else {
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  }

  cat("Parameters:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)

  cat("\nLog-likelihood:", format(signif(x$loglik, digits)),
    "(df =", paste(x$df, ")", sep = ""), "\n\n")
  invisible(x)
}


coef.btmodel <- function(object, all = TRUE, ref = !all, ...) {
  lab <- object$labels
  nobj <- length(lab)
  acf <- object$coefficients
  cf <- structure(rep(0, nobj), .Names = lab)
  ocf <- acf[1:(nobj-1)]
  cf[names(ocf)] <- ocf
  cf <- c(cf, acf[-(1:(nobj-1))])
  if(!all) cf <- cf[1:nobj]
  if(!ref) cf <- cf[-match(object$ref, lab)]
  return(cf)
}

worth <- function(object, ...) UseMethod("itempar")

plot.btmodel <- function(x, 
  worth = TRUE, index = TRUE, names = TRUE, ref = TRUE, abbreviate = FALSE,
  type = NULL, lty = NULL, xlab = "Objects", ylab = NULL, ...)
{
  ## parameters to be plotted
  if(worth) {
    cf <- worth(x)
  } else {
    cf <- coef(x, all = FALSE, ref = TRUE)
    if(is.character(ref) | is.numeric(ref)) {
      reflab <- ref
      ref <- TRUE
    } else {
      reflab <- x$ref
    }
    if(is.character(reflab)) reflab <- match(reflab, x$labels)
    cf <- cf - cf[reflab]
  }

  ## labeling
  cf_ref <- if(!worth) 0 else 1/length(cf)
  if(is.character(names)) {
    names(cf) <- names
    names <- TRUE
  }
  if(is.null(ylab)) ylab <- if(worth) "Worth parameters" else "Parameters"
    
  ## abbreviation
  if(is.logical(abbreviate)) {
    nlab <- max(nchar(names(cf)))
    abbreviate <- if(abbreviate) as.numeric(cut(nlab, c(-Inf, 1.5, 4.5, 7.5, Inf))) else nlab
  }
  names(cf) <- abbreviate(names(cf), abbreviate)

  ## raw plot
  ix <- if(index) seq(along = cf) else rep(0, length(cf))
  plot(ix, cf, xlab = xlab, ylab = ylab, type = "n", axes = FALSE, ...)
  if(ref) abline(h = cf_ref, col = "lightgray")
  axis(2)
  box()  

  ## actual data
  if(index) {
    if(is.null(type)) type <- "b"
    if(is.null(lty)) lty <- 2  
    lines(ix, cf, type = type, lty = lty)
    axis(1, at = ix, labels = if(names) names(cf) else TRUE)
  } else {
    if(is.null(type)) type <- "p"
    if(names) text(names(cf), x = ix, y = cf, ...) else lines(ix, cf, type = type, lty = lty)
  }
}

itempar.btmodel <- function (object, ref = NULL, alias = TRUE, vcov = TRUE, log = FALSE, ...) {

  ## extract cf and labels, include restricted (ref = TRUE) but not (if existent) undecided parameter (all = FALSE)
  cf <- coef(object, ref = TRUE, all = FALSE)
  ncf <- length(cf)
  ncfv <- 1:ncf
  lbs <- names(cf)

  ## process ref
  if (is.null(ref)) {
    ref <- ncfv
  } else if (is.vector(ref) && is.character(ref)) {
    stopifnot(all(ref %in% lbs))
    ref <- which(lbs %in% ref)
  } else if (is.vector(ref) && is.numeric(ref)) {
    ref <- as.integer(ref)
    stopifnot(all(ref %in% ncfv))
  } else if (is.matrix(ref)) {
    stop("Handling of contrast matrices in argument 'ref' currently not implemented for itempar.btmodel().")
  } else stop("Argument 'ref' is misspecified (see ?itempar for possible values).")

  ## specify contrast matrix
  D <- diag(ncf)
  D[, ref] <- D[, ref] - 1/length(ref)

  ## impose contrast
  cf <- as.vector(D %*% cf)

  ## create adjusted vcov if requested
  if (vcov) {
    vc <- vcov(object)
    vc <- vc[setdiff(rownames(vc), "(undecided)"), setdiff(colnames(vc), "(undecided)")]
    vc <- D %*% rbind(0, cbind(0, vc)) %*% t(D)
  } else {
    vc <- matrix(NA, nrow = ncf, ncol = ncf)
  }

  ## if not on log scale ... 
  if (!log) {

    ## transform parameters as requested
    cf <- exp(cf)
    scf <- sum(cf[ref])
    cf <- cf/scf

    ## create adjusted vcov if requested
    if (vcov) {
      ## specify contrast matrix
      D <- matrix(0, nrow = ncf, ncol = ncf) # i != j and j not in ref -> 0
      for (i in ncfv) {
        if (i %in% ref) {
          D[i, i] <- (scf - cf[i])/(scf^2) # i = j and j in ref
          D[-i, i] <- -cf[-i]/(scf^2)      # i != j and j in ref
        } else {
          D[i, i] <- 1/scf                 # i = j and j not in ref
        }
      }

      ## adjust vcov
      vc <- diag(cf) %*% vc %*% t(diag(cf))
      vc <- D %*% vc %*% t(D)
    }

  }

  ## set labels
  names(cf) <- rownames(vc) <- colnames(vc) <- lbs

  ## process argument alias
  if (!alias) {
    cf <- cf[-ref[1]]
    vc <- vc[-ref[1], -ref[1]]
    alias <- ref[1]
    names(alias) <- lbs[ref[1]]
  }
  
  ## setup and return result object
  rv <- structure(cf, class = "itempar", model = "btmodel", ref = ref, alias = alias, vcov = vc)
  return(rv)
}

nobs.btmodel <- function (object, ...)
{
  return(object$n)
}

estfun.btmodel <- function(x, ...)
{
  ## setup relevant informations
  stopifnot(inherits(x$y, "paircomp"))
  y <- x$y
  ymat <- as.matrix(y)
  lab <- x$lab
  ref <- which(lab == x$ref)
  undecided <- x$undecided
  type <- x$type
  nobj <- length(lab)
  npc <- nobj * (nobj - 1)/2
  par <- x$coef
  npar <- nobj - !undecided
  ix <- which(upper.tri(diag(nobj)), arr.ind = TRUE)
  nam <- c(lab[-ref], if(undecided) "(undecided)" else NULL)
  weights <- if(is.null(x$weights)) 1L else x$weights
  ytab <- summary(y, weights = x$weights)
  mscale <- mscale(y)
  has_ties <- mscale[2L] == 0L
  if(has_ties) ytab <- ytab[, c(1L, 3L, 2L), drop = FALSE]
  if(!undecided) ytab <- ytab[, 1L:2L, drop = FALSE]

  ## log-probabilities and log-likelihood
  par2logprob <- switch(type,
    "loglin" = function(i) {
      p <- rep(0, nobj)
      p[-ref] <- if(undecided) par[-npar] else par
      p <- p[ix[i,]]
      if(undecided) p <- c(p, par[npar] + mean(p))
      p - log(sum(exp(p)))},
    "logit" = function(i) {
      p <- rep(0, nobj)
      p[-ref] <- par
      plogis(c(-1, 1) * diff(p[ix[i,]]), log.p = TRUE)
    })
  logp <- t(sapply(1L:npc, par2logprob))
  loglik <- sum(logp * ytab)

  ## estimating functions
  if(!undecided) logp <- cbind(logp, -Inf) ## ties impossible
  gradp <- matrix(0, nrow = npc * 3, ncol = nobj)
  cf <- -matrix(c(1, 0, 0, 0, 1, 0, 0.5, 0.5, 1), ncol = 3)
  ct <- matrix(c(1, 0, 0.5, 0, 1, 0.5, 0, 0, 1), ncol = 3)
  for(i in 1L:npc) gradp[i*3 - (2:0), c(ix[i,], nobj)] <- t(t(ct) + as.vector(cf %*% exp(logp[i,])))
  ef <- t(sapply(1L:length(y), function(i) {
    wi <- (0:(npc - 1)) * 3 + c(2, 3, 1)[ymat[i,] + 2]
    colSums(gradp[wi,, drop = FALSE], na.rm = TRUE)
  }))
  if(!undecided) ef <- ef[, -nobj, drop = FALSE]
  dimnames(ef) <- list(names(y), nam)
  if(!is.null(weights)) ef <- ef * weights  

  return(ef)
}
