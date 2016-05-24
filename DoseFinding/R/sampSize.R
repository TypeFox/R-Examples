## function for sample size calculation and functions to evaluate
## performance metrics for different sample sizes

sampSize <- function (upperN, lowerN = floor(upperN/2),
                      targFunc, target, tol = 0.001, alRatio,
                      Ntype = c("arm", "total"), verbose = FALSE){
  ## target function to iterate
  func <- function(n){
    targFunc(n) - target
  }

  Ntype <- match.arg(Ntype)
  if (!missing(alRatio)) {
    if (any(alRatio <= 0)) {
      stop("all entries of alRatio need to be positive")
    } else {
      alRatio <- alRatio/sum(alRatio)
    }
    if(Ntype == "arm") {
      alRatio <- alRatio/min(alRatio)
    } 
  } else { ## by default assume
    stop("allocation ratios need to be specified")
  }
  
  ## first call
  upper <- func(round(upperN*alRatio))
  if(length(upper) > 1)
    stop("targFunc(n) to evaluate to a vector of length 1.")
  if(!is.numeric(upper))
    stop("targFunc(n) needs to evaluate to a numeric.")

  ## bracket solution
  if (upper < 0)
    message("upper limit for sample size is raised")

  while (upper < 0) {
    upperN <- 2 * upperN
    upper <- func(round(upperN*alRatio))
  }
  
  lower <- func(round(lowerN*alRatio))
  
  if (lower > 0) 
    message("lower limit for sample size is decreased")

  while (lower > 0) {
    lowerN <- round(lowerN/2)
    if (lowerN == 0) 
      stop("cannot find lower limit on n")
    lower <- func(round(lowerN*alRatio))
  }

  ## now start bisection
  if (verbose) {
    cat("Upper N:", upperN, "Upper value", round(upper+target, 4), "\n")
    cat("Lower N:", lowerN, "Lower value", round(lower+target, 4), "\n\n")
  }
  
  current <- tol+1
  niter <- 0
  ## bisect sample size until tolerance is achieved
  while (abs(current) > tol & (upperN > lowerN + 1)) {
    currN <- round((upperN + lowerN)/2)
    current <- func(round(currN * alRatio))
    if (current > 0) {
      upperN <- currN
    } else {
      lowerN <- currN
    }
    niter <- niter + 1
    if (verbose) {
      cat("Iter: ", niter, ", N = ", currN, ", current value = ",
          round(current+target, 4), "\n", sep = "")
    }
  }
  ## increase sample size so that the obtained value is larger than the target
  while (current < 0) {
    currN <- currN + 1
    current <- func(round(currN * alRatio))
  }

  res <- list(samp.size = round(currN * alRatio),
              target = round(current+target, 4))
  attr(res, "alRatio") <- round(alRatio/min(alRatio), 4)
  attr(res, "target") <- target
  attr(res, "Ntype") <- Ntype
  class(res) <- "sampSize"
  res
}

print.sampSize <- function(x, ...){
  cat("Sample size calculation\n\n")
  cat("alRatio:", attr(x, "alRatio"), "\n")
  cat("Total sample size:", sum(x$samp.size), "\n")
  cat("Sample size per arm:", x$samp.size, "\n")
  cat("targFunc:", x$target,"\n")
}

sampSizeMCT <- function(upperN, lowerN = floor(upperN/2),
                        ...,
                        power, sumFct = mean,
                        tol = 0.001, alRatio, Ntype = c("arm", "total"), verbose = FALSE){
  ## function to calculate sample size for multiple contrast test
  ## if S is specified this needs to be the (hypothetical) covariance matrix
  ## for a total sample size of 1 patient
  Ntype <- match.arg(Ntype)
  args <- list(...)
  namargs <- names(args)
  if(is.element("placAdj", namargs)){
    if(args$placAdj)
      stop("placAdj needs to be FALSE for sampSizeMCT.
  Use sampSize directly in placebo-adjusted case.")
  }
  if(is.element("S", namargs)){
    S <- args[["S"]]
    if(Ntype == "arm"){
      Ntype <- "total"
      message("Only Ntype == \"total\" possible if S is specified")
    }
    if(is.element("df", namargs)){
      if(is.finite(args$df))
        message("df argument set to Inf, if S is specified.
Use sampSize directly in case exact df are required.")
    }
    args$df <- Inf
    tFunc <- function(n){
      N <- sum(n)
      Sn <- 1/N*S
      args$S <- Sn
      powVals <- do.call("powMCT", args)
      sumFct(powVals)
    }
  } else {
    if(is.element("n", namargs))
      stop("n is not allowed to be specified for sample size calculation")
    if(!is.element("sigma", namargs))
      stop("need sigma if S is not specified")
    tFunc <- function(n){
      powVals <- powMCT(n=n, ...)
      sumFct(powVals)
    }
  }
  sampSize(upperN, lowerN, targFunc = tFunc, target = power,
           alRatio = alRatio, Ntype = Ntype, verbose = verbose)
}

targN <- function(upperN, lowerN, step, targFunc,
                  alRatio, Ntype = c("arm", "total"), sumFct = c("min", "mean", "max")){

  if(!is.character(sumFct))
    stop("sumFct needs to be a character vector")
  Ntype <- match.arg(Ntype)
  if (!missing(alRatio)) {
    if (any(alRatio <= 0)) {
      stop("all entries of alRatio need to be positive")
    } else {
      alRatio <- alRatio/sum(alRatio)
    }
    if(Ntype == "arm") {
      alRatio <- alRatio/min(alRatio)
    } 
  } else { ## by default assume 
    stop("allocation ratios need to be specified")
  }
  
  nseq <- seq(lowerN, upperN, by=step)
  out <-t(sapply(nseq, function(x){
    targFunc(round(x * alRatio))
  }))
  if(nrow(out) == 1 & length(nseq) > 1){
    out <- t(out)
    colnames(out) <- ""
  }
  out2 <- out
  for(i in 1:length(sumFct)){
    out2 <- cbind(out2, apply(out, 1, sumFct[i]))
  }
  dimnames(out2) <- list(nseq, c(colnames(out), sumFct))
  attr(out2, "alRatio") <- alRatio
  attr(out2, "sumFct") <- sumFct
  attr(out2, "Ntype") <- Ntype
  class(out2) <- "targN"
  out2
}

powN <- function(upperN, lowerN, step,
                 ...,
                 alRatio, Ntype = c("arm", "total"), sumFct = c("min", "mean", "max")){
  args <- list(...)
  namargs <- names(args)
  if(is.element("placAdj", namargs)){
    if(args$placAdj)
      stop("placAdj needs to be FALSE for powN.
  Use targN directly in placebo-adjusted case.")
  }
  Ntype <- match.arg(Ntype)
  if(is.element("S", namargs)){
    S <- args[["S"]]
    if(Ntype == "arm"){
      Ntype <- "total"
      message("Only Ntype == \"total\" possible if S is specified")
    }
    if(is.element("df", namargs)){
      if(is.finite(args$df))
        message("df argument set to Inf, if S is specified.
Use sampSize directly in case exact df are required.")
    }
    args$df <- Inf
    tFunc <- function(n){
      N <- sum(n)
      Sn <- 1/N*S
      args$S <- Sn
      do.call("powMCT", args)
    }
  } else {
    if(is.element("n", namargs))
      stop("n is not allowed to be specified for sample size calculation")
    if(!is.element("sigma", namargs))
      stop("need sigma if S is not specified")
    tFunc <- function(n)
      powMCT(n=n, ...)
  }
  targN(upperN=upperN, lowerN=lowerN, step=step, targFunc=tFunc,
        alRatio=alRatio, Ntype = Ntype, sumFct = sumFct)
}

## Produces Trellis plot of targN object
plot.targN <- function(x, superpose = TRUE, line.at = NULL, 
                       xlab = NULL, ylab = NULL, ...){
  nSeq <- as.integer(dimnames(x)[[1]])
  alRatio <- attr(x, "alRatio")
  unbN <- (length(unique(alRatio)) > 1)
  if (is.null(xlab)) {
    if(attr(x, "Ntype") == "total" | unbN){
      xlab <- "Overall sample size"
      nSeq <- sapply(nSeq, function(x){
        sum(round(x*alRatio))
      })
    } else {
      xlab <- "Sample size per dose (balanced)"
    }
  }
  nams <- dimnames(x)[[2]]
  ## separating model data from summary data
  x <- as.data.frame(unclass(x))
  nams <- names(x)
  nC <- ncol(x)
  pMatTr <- data.frame(targ = as.vector(unlist(x)), n = rep(nSeq, nC),
                       type = factor(rep(nams, each = length(nSeq)), levels = nams))
  if(superpose){
    panelFunc1 <- function(x, y, subscripts, groups, lineAt, ...) {
      panel.grid(h = -1, v = -1, col = "lightgrey", lty = 2)
      if(!is.null(line.at))
        panel.abline(h = lineAt, lty = 3, ..., col = "red")
      panel.superpose(x, y, subscripts, groups, ...)
    }
    trLn <- trellis.par.get("superpose.line")[c("col", "lwd", "lty")]
    for(i in seq(along = trLn)) {
      if(length(trLn[[i]]) > nC) trLn[[i]] <- trLn[[i]][1:nC]
    }
    ltplot <- xyplot(targ ~ n, pMatTr, groups = pMatTr$type, subscripts = TRUE,
                     panel = panelFunc1, type = "l", lineAt = line.at,
                     xlab = xlab, ylab = ylab,
                     key = list(lines = trLn, text = list(lab = nams), transparent = TRUE, 
                       columns = ifelse(nC < 5, nC, min(4,ceiling(nC/min(ceiling(nC/4),3))))), ...)
  } else {                              # models in different panels
    panelFunc2 <- function(x, y, lineAt, ...) {
      panel.grid(h = -1, v = -1, col = "lightgrey", lty = 2)
      if(!is.null(line.at))
        panel.abline(h = lineAt, lty = 3, ..., col = "red") ## used 2 for consistency with above
      panel.xyplot(x, y, ...)
    }
    ltplot <- xyplot(targ ~ n | type, pMatTr, panel = panelFunc2,
                     type = "l", lineAt = line.at,
                     xlab = xlab, ylab = ylab, 
                     strip = function(...) strip.default(..., style = 1), ...)
  }
  print(ltplot)
}
