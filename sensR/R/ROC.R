SDT <- function(tab, method = c("probit", "logit")) {
  ## Performs the "emperical probit or logit transform" on a 2 x J table of
  ## observations. Also called the Signal Detection Theory computation
  ## of d-prime
  tmp <- function(x) {
    cs <- cumsum(x)
    l <- length(cs)
    if(method == "probit") {
      cp <- cs/cs[l]
      qnorm(cp[-l])}
    else if(method == "logit")
      log( cs[-l] / (cs[l]-cs[-l]))
    else stop("method", method, "not recognized")
  }
  m <- match.call(expand.dots=FALSE)
  m$method <- NULL
  m[[1]] <- as.name("list")
  m <- eval.parent(m)
  mat <- m$tab
  stopifnot(is.matrix(mat))
  method <- match.arg(method)

  rn <- rownames(mat)
  sl <- apply(mat, 1, tmp)
  if(!is.matrix(sl)) sl <- matrix(sl, nrow=1)
  ## name <- paste(rn[1], "-", rn[2])
  d.sl <- cbind(sl, sl[,1]-sl[,2])
  colnames(d.sl) <- c("z(Hit rate)", "z(False alarm rate)", "d-prime")
  rownames(d.sl) <- 1:(ncol(mat)-1)
  d.sl
}


ROC.default <-
  function (object, se.d, scale = 1, length = 1000, fig = TRUE,
            se.type = c("CI", "SE"), CI.alpha = 0.05, ...)
{
    m <- match.call(expand.dots = FALSE)
    m$se.type <- NULL
    m[[1]] <- as.name("list")
    eval.parent(m)
    d <- object
    se.type <- match.arg(se.type)
    fpr <- seq(0, 1, length.out = length)
    q.fpr <- qnorm(fpr)
    tpr <- pnorm((q.fpr + d)/scale)
    res <- list(ROCx = fpr, ROCy = tpr)
    if(fig) {
      plot.default(fpr, tpr, type = "l", xlim = c(0, 1),
                   ylim = c(0, 1),
                   xlab = "True positive ratio",
                   ylab = "False positive ratio",
                   main = "", ...)
      lines(fpr, fpr, lty = 2)
      if (!missing(se.d)) {
        if (se.type == "CI")
          tol <- se.d * qnorm(1 - CI.alpha/2)
        else if (se.type == "SE")
          tol <- se.d
        lower <- pnorm((q.fpr + d - tol)/scale)
        upper <- pnorm((q.fpr + d + tol)/scale)
        lines(fpr, lower, lty = 3, ...)
        lines(fpr, upper, lty = 3, ...)
        res$lower <- lower
        res$upper <- upper
      }
    }
    invisible(res)
}

ROC.anota <-
  function(object, length = 1000, fig = TRUE,
           se.type = c("CI", "SE"), CI.alpha = .05, ...)
{
  stopifnot(object$test == "A-Not A")
  ROC.default(object$coef, object$se, , length, fig, se.type,
              CI.alpha, ...)
}

ROC <- function(object, ...) {
  UseMethod("ROC")
}

AUC <- function(d, ...) {
  UseMethod("AUC")
}

AUC.default <- function(d, se.d, scale = 1, CI.alpha = .05, ...) {
    stopifnot(is.numeric(d),
              length(d) == 1L)
### NOTE: We allow negative d-primes here. Confidence intervals for
### d-prime also work for negative d-primes. Note that AnotA can also
### report negative d-primes.
    stopifnot(is.numeric(scale),
              length(scale) == 1L,
              scale > 0)
    stopifnot(is.numeric(CI.alpha),
              length(CI.alpha) == 1L,
              CI.alpha > 0,
              CI.alpha < 1)
    if(!missing(se.d) && !is.null(se.d)) {
        stopifnot(is.numeric(se.d),
                  length(se.d) == 1L,
                  se.d >= 0)
    }
    ## Compute AUC:
    Scale <- 1 + scale^2
    res <- list(value = pnorm(d / sqrt(Scale)))
    if(!missing(se.d) && !is.null(se.d)) {
        tol <- se.d * qnorm(1 - CI.alpha/2)
        res$lower <- pnorm((d - tol)/sqrt(Scale))
        res$upper <- pnorm((d + tol)/sqrt(Scale))
        res$CI.alpha = CI.alpha
    }
    class(res) <- "AUC"
    res
}

print.AUC <- function(x, digits = getOption("digits"), ...){
  cat(paste("AUC:", signif(x$value, digits)), "\n")
  if(!is.null(x$lower)) {
    cat(paste(1-x$CI.alpha,"% CI: [", signif(x$lower, digits),
                ", ", signif(x$upper, digits), "]", sep=""), "\n")
  }
  invisible()
}

AUC.anota <- function(d, CI.alpha = .05, ...) {
  ## stopifnot(d$test == "A-Not A")
  AUC.default(d=d$coef, se.d=d$se, CI.alpha=CI.alpha, ...)
}
