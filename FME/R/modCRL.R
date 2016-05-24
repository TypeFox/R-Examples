## -----------------------------------------------------------------------------
## Monte Carlo runs
## -----------------------------------------------------------------------------

modCRL <- function(func, parms = NULL, sensvar = NULL, dist = "unif",
                   parInput = NULL, parRange = NULL, parMean = NULL,
                   parCovar = NULL, num = 100, ...) {

  if (is.null(parms) & ! is.null(parInput))
    parms <- parInput[1,]
  if (is.null(parms))
    parms <- parMean
  if (is.null(parms)) {
    if (is.vector(parRange))
      parms <- mean(parRange)
    else parms <- rowMeans(parRange)
  }
  if (is.null(parms))
    stop ("'parms' not known")
  if (is.matrix(parms) && nrow(parms)>1)
    stop ("'parms' should be a vector")

  Solve <- function(parms) func(parms, ...)

  if (! is.null(parInput)) {
    dist <- "input"
    nr   <- nrow(parInput)
    num  <- min(num,nr)
    if (num == nr)
      ii <- 1:nr
    else ii <- sample(1:nr, size = num, replace = FALSE)

    parset <- as.matrix(parInput[ii,])
    if(is.null(parms)) parms <- parInput[1,]
  }

  Parms <- parms

  ## reference run
  yRef <- func(parms, ...)

  if (is.matrix(yRef) | is.data.frame(yRef))
    if(nrow(yRef)>1)
      stop("func should return a vector or a matrix/data.frame with one row")

  sens <- sensRange(func, parms, sensvar, dist, parInput, parRange, parMean,
                  parCovar, map = NULL, num = num, ...)

  class(sens) <- c("modCRL", "data.frame")
  return (sens)
}

## -----------------------------------------------------------------------------
## S3 methods of modCRL
## -----------------------------------------------------------------------------

summary.modCRL <- function(object, ...) {
  SumSens <- summary.sensRange(object, ...)
  class(SumSens) <- c("summary.modCRL", "data.frame")
  return(SumSens)
}

## -----------------------------------------------------------------------------
hist.modCRL <- function(x, which = 1:ncol(x), ask = NULL, ...) {
  hh <- list()
  hh$pars <- x
  hist.modMCMC(hh, which = which, Full = FALSE, ask = ask, ...)
}

## -----------------------------------------------------------------------------
pairs.modCRL <- function(x, which = 1:ncol(x), nsample = NULL, ...) {

  panel.main <- function(x, y, ...)
    points(x[ii], y[ii], ...)

  var   <- colnames(x)
  which <- selectvar(which, var, "x", Nall = TRUE)

  X <- x[,which]
  X <- as.matrix(X)

  if (is.null(nsample))
    ii <- 1:nrow(X) else
    ii <- sample((1:nrow(X)), nsample)

  labels <- colnames(X)
  dots <- list(...)

  dots$diag.panel  <- if(is.null(dots$diag.panel))  panel.hist else dots$diag.panel
  dots$lower.panel <- if(is.null(dots$lower.panel)) panel.cor  else dots$lower.panel
  dots$upper.panel <- if(is.null(dots$upper.panel)) panel.main else dots$upper.panel
  dots$gap <- if(is.null(dots$gap)) 0 else dots$gap
  dots$labels <- if(is.null(dots$labels)) labels else dots$labels

  do.call("pairs", c(alist(X), dots))
}

## -----------------------------------------------------------------------------
plot.modCRL<-function(x, which = NULL, trace = FALSE, ask = NULL, ...) {

  vars <- attr(x,"var")

  var <- colnames(x)
  NP  <- attr(x,"npar")
  parnames <- var[1:NP]
  varnames <- var[-(1:NP)]

  ## which selection has to stay the way it is: selecting from TWO sets...
  if (!is.null(which)) {
    if (! is.numeric(which)) {
      ln <- length(which)
      Select <- NP + which(varnames %in% which)
      SelPar <- which(parnames %in% which)
      if(length(Select) + length(SelPar) != ln)
        stop("not all variables or parameters in 'which' are in 'x'")
    } else {
      Select <- which[which > NP]
      SelPar <- which[which <= NP]
      if (max(which) > length(var))
        stop("index in 'which' too large")
    }
  } else {
    Select <- NP+(1:length(varnames))
    SelPar <- 1:NP
  }
  nvar <- length(Select)
  np   <- length(SelPar)
  if (np == 0)
    stop("cannot plot Monte Carlo results: select at least one parameter")
  if (nvar == 0)
    stop("cannot plot Monte Carlo results: select at least one variable")

  if (np > 1 & nvar > 1)
    stop("cannot plot Monte Carlo results for more than 1 parameter and more than one variable")

  pars   <- attr(x, "names")
  dots   <- list(...)
  nmdots <- names(dots)

  nv <- max(nvar, np)
  ## Set par mfrow and ask.
  ask <- setplotpar (nmdots, dots, nv, ask)

  ## interactively wait if there are remaining figures
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  Main <- is.null(dots$main)
  Xlab <- is.null(dots$xlab)
  dots$ylab <- if(is.null(dots$ylab))  ""  else dots$ylab

  if (np == 1) {
    ip <- SelPar[1]
    dots$xlab <- if(is.null(dots$xlab)) pars[ip] else dots$xlab

    for(i in Select) {
      if (Main) dots$main <- var[i]
      do.call("plot", c(alist(x[,ip], x[,i]), dots))
      if (trace) lines(lowess(x[,ip], x[,i]), lwd = 2)
    }
  } else {  # nvar==1
    iv <- Select[1]
    for(i in SelPar) {
      if (Xlab) dots$xlab <- pars[i]
      if (Main) dots$main <- var[iv]
      do.call("plot", c(alist(x[,i], x[,iv]), dots))
      if (trace) lines(lowess(x[,i], x[,iv]), lwd = 2)
    }
  }
}
