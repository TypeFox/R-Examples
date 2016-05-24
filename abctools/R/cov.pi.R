getabcout <- function(x) {
  if (x$method == "rejection") {
    y <- x$unadj.values
  } else {
    y <- x$adj.values
  }
  if (!is.matrix(y)) {
    y <- matrix(y, nrow=1)
  }
  return(y)
}

covstats.pi <-
function(raw, diagnostics=c("KS", "CGR"), nacc.min=20) {
  usetol <- "tol" %in% colnames(raw)
  tosplit <- ifelse(usetol, "tol", "eps")
  diag <- NULL
  ##COMPUTE DIAGNOSTICS
  if ("KS" %in% diagnostics) {
    diag <- dlply(raw, tosplit, ##I should be able to use ddply but it produces errors...
                  function(x) {
                    xx <- NULL
                    if (any(x$nacc<nacc.min) || sum(complete.cases(x)) < nrow(x)) {
                      xx <- rep(NA, ncol(x)-3)
                    }
                    x <- x[,-3] ##Remove nacc entries
                    if (is.null(xx)) {
                      xx <- apply(x[,-(1:2), drop=FALSE], 2,
                                  function(y){
                                    suppressWarnings(ks.test(y, punif, exact=FALSE))$p.value ##Warnings of ties suppressed
                                  })
                    }
                    return(data.frame(tol=x[1,2], parameter=colnames(x)[-(1:2)], pvalue=xx))
                  })
    diag <- rbind.fill(diag)
    diag <- cbind(diag, test="KS")
  }
  if ("CGR" %in% diagnostics) {
    temp <- dlply(raw, tosplit, ##I should be able to use ddply but it produces errors...
                  function(x) {
                    xx <- NULL
                    if (any(x$nacc<nacc.min) || sum(complete.cases(x)) < nrow(x)) {
                      xx <- rep(NA, ncol(x)-3)
                    }
                    x <- x[,-3] ##Remove nacc entries
                    if (is.null(xx)) {
                      xx <- apply(x[,-(1:2), drop=FALSE], 2,
                                  function(y){
                                    z <- qnorm(y)
                                    z <- sum(z^2)
                                    p1 <- pchisq(z, df=length(y)) ##p-value for one-tailed test
                                    1-2*abs(p1-0.5) ##p-value for two-tailed test
                                  })
                    }
                    return(data.frame(tol=x[1,2], parameter=colnames(x)[-(1:2)], pvalue=xx))
                  })
    temp <- rbind.fill(temp)
    temp <- cbind(temp, test="CGR")
    if (is.null(diag)) {
      diag <- temp
    } else {
      diag <- rbind(diag, temp)
    }
  }
  if (!is.null(diag) && !usetol) {
    names(diag)[1] <- "eps"
  }
  return(diag)
}

cov.pi <-
function(param, sumstat, testsets, tol, eps, diagnostics=c(), multicore=FALSE, cores, method="rejection", nacc.min=20, ...) {
  ##CHECK INPUT
  if (!is.data.frame(param)) stop("param must be a data frame")
  if (!is.data.frame(sumstat)) stop("sumstat must be a data frame")
  if (nrow(param) != nrow(sumstat)) stop("Number of rows in sumstat and param must be equal")
  if (missing(tol) && missing(eps)) stop("Either tol or eps must be specified")
  ##INITIALISE SOME VALUES
  n <- nrow(param)
  usetol <- !missing(tol)
  ##PREPARE INPUT VALUES
  if (usetol) {
    xx <- expand.grid(testset=testsets, tol=tol, stringsAsFactors=FALSE)
  } else {
    xx <- expand.grid(testset=testsets, eps=eps, stringsAsFactors=FALSE)
  }
  ##FUNCTION TO DO ONE ITERATION IE ABC FOR ONE INPUT CHOICE
  doit <- function(i){
    set <- xx$testset[i]
    pseudo.theta0 <- param[set,,drop=TRUE]
    if (usetol) {
      mytol <- xx$tol[i]
    } else {
      temp <- abc(target=sumstat[set,,drop=FALSE], param=param[-set,,drop=FALSE], sumstat=sumstat[-set,,drop=FALSE], tol=1, method="rejection", ...)
      mytol <- mean(temp$dist <= xx$eps[i])
    }
    ##CHECK POST-PROCESSING POSSIBLE
    if (method != "rejection") {
      nacc <- floor(mytol * nrow(param))
      if (nacc < max(ncol(sumstat), ncol(param))) {
        out <- rep(NA, ncol(param))
        names(out) <- colnames(param)
        out <- c(nacc=nacc, out)
        return(out)
      }
    }
    temp <- abc(target=sumstat[set,,drop=FALSE], param=param[-set,,drop=FALSE], sumstat=sumstat[-set,,drop=FALSE], tol=mytol, method=method, ...)
    temp <- getabcout(temp)
    out <- sapply(1:length(pseudo.theta0),
                  function(j){
                    (1 + sum(temp[,j] < pseudo.theta0[j])) / (2 + nrow(temp)) ##A robust cdf estimate (equiv to mean of Bayesian binomial rate estimate with uniform prior)
                  })
    names(out) <- colnames(param)
    out <- c(nacc=nrow(temp), out)
    return(out)
  }
  if (multicore) {
    raw <- mclapply(1:nrow(xx), doit, mc.cores=cores)
    raw <- simplify2array(raw)
  } else {
    raw <- sapply(1:nrow(xx), doit)
  }
  if (!is.matrix(raw)) {
    temp <- matrix(raw, nrow=1)
    rownames(temp) <- names(raw)[1]
    raw <- temp
  }
  raw <- cbind(xx,t(raw))
  ##RETURN OUTPUT
  return(list(raw=raw, diag=covstats.pi(raw, diagnostics, nacc.min)))
}
