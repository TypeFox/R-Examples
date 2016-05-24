
## -----------------------------------------------------------------------------
## Sensitivity functions
## -----------------------------------------------------------------------------

sensFun <- function(func, parms, sensvar = NULL, senspar = names(parms),
                    varscale = NULL, parscale = NULL,
                    tiny = 1e-8, map = 1, ...) {
  ## 1. The solver
  Solve <- function(parms) func(parms, ...)

  yRef  <- Solve(parms)
  Type <- 1
  if("modCost" %in% class(yRef))  {
    Res    <- yRef$residuals
    ynames <- Res$name
    YNAMES <- as.matrix(unique(ynames))                                                 # changed on 01-12-2011 (Tom): rearrangement of residual matrix to put costs of same var together
    iRes <- unlist(apply(YNAMES,MARGIN=1,function(x) which(Res[,"name"]==x)))              # (cf. modcost(... cost=cost) where residuals for different cost calculations are appended but not 
    yRef <- cbind(Res$x[iRes], Res$mod[iRes])                                           # sorted by variable)
    ynames <- ynames[iRes]
    names(yRef) <- ynames
    Type <- 2
    sensvar <- NULL   # input of sensvar not allowed for modCost type

    Solve <- function(parms) {
      Res<- func(parms, ...)$residuals
      YNAMES <- as.matrix(unique(Res[,"name"]))                                           # changed on 01-12-2011 (Tom): for same reason as above
      iRes <- unlist(apply(YNAMES,MARGIN=1,function(x) which(Res[,"name"]==x)))              # 
      cbind(Res$x[iRes], Res$mod[iRes])
    }

  }
  ## if a data.frame or a vector is returned, make it a matrix
  if (is.data.frame(yRef)) yRef <- as.matrix(yRef)
  if (is.vector(yRef)) {
    ynames <- names(yRef)
    yRef <- matrix(data=yRef, nrow = 1)
    colnames(yRef) <- ynames
  }

  ## 2. sensitivity variables
  if (is.null(sensvar)) {
    ivar <- 1:ncol(yRef)
    if (! is.null(map))
      ivar <- ivar[-map]
    sensvar <- colnames(yRef)[ivar]
    if(is.null(sensvar))
      sensvar <- ivar
  } else {
    ivar  <- findvar(yRef[1,], sensvar, "variables")
    if (! is.character(sensvar)) { # try to create names rather than nrs
      sv <- sensvar
      sensvar <- colnames(yRef)[ivar]
      if (is.null(sensvar))
        sensvar <- sv
    }
  }
  if (is.null(map)) {
    map   <- 1:nrow(yRef)
    mname <- "x"
    }
  else {
    mname <- colnames(yRef)[map]
    if (is.null(mname)) mname <- "x"
    map <- yRef[, map]
  }

  nout  <- length(ivar)
  ndim  <- nrow(yRef)
  if (Type == 1)
    grvar <- expand.grid(map, sensvar)
  else grvar<-data.frame(x = map, var = ynames)

  if (ndim ==1)
    svar <- sensvar
  else svar <- paste(grvar[, 2], grvar[, 1], sep = "")

  yRef <- as.vector(yRef[, ivar])

  if (is.null(senspar))
    senspar <- 1:length(parms)

  ## 3. sensitivity parameters/
  npar  <- length(senspar)
  if (npar == 0)
    stop ("cannot proceed: there are no sensitivity parameters")

  ipar <- findvar(parms, senspar, "parameters")
  pp    <- unlist(parms)[ipar]

  ## 4. perturbed parameters - perturbations always positive
  dp          <- abs(pp*tiny)
  dp[dp == 0] <- tiny
  ii          <- which (dp < tiny)
  dp[ii]      <- tiny

  if (is.null(parscale))
    parscale <- pp
  else parscale<-rep(parscale, length.out=npar)              #changed 21/10/2011 (Karline, Tom)

  if (is.null(varscale))
    varscale <- yRef
  else varscale <- rep (varscale, length.out=length(yRef))   #changed 21/10/2011 (Karline, Tom)

  ## 0 is set equal to a very small number
  varscale[varscale == 0] <- tiny*1e-12
  parscale[parscale == 0] <- tiny*1e-12

  Sens    <- matrix(data=NA, nrow = length(yRef), ncol = npar)

  ## 5. Loop over all parameters
  for (i in 1:length(ipar)) {
    dval    <- pp[i] + dp[i]
    parms[ipar[i]] <- dval
    Yres    <- Solve(parms)
    if (is.vector(Yres))
      yPert <- Yres[ivar]
    else yPert <- as.vector(unlist(Yres[, ivar]))
    Sens[, i] <- (yPert-yRef)/dp[i] * parscale[i] /varscale
    parms[ipar[i]] <- pp[i]
  }

  ## 6. Finally
  colnames(Sens) <- names(pp)
  Sens <- data.frame(x = grvar[, 1], var = as.character(grvar[, 2]), Sens)
  attr(Sens, "class") <- c("sensFun", "data.frame")
  attr(Sens, "pars") <- pp
  attr(Sens, "parscale") <- parscale
  attr(Sens, "varscale") <- varscale
  if (Type == 2) {
    attr(Sens, "var") <- as.vector(unique(ynames))
    attr(Sens, "nx")  <- c(0,cumsum(as.vector(table(ynames))))  #start of each var
  } else {
    attr(Sens, "var") <- sensvar
    attr(Sens, "nx")  <- length(map)
  }
  attr(Sens, "xname") <- mname
  attr(Sens, "x")     <- map
  attr(Sens, "Type" ) <- Type  # type 1: modCost
  return(Sens)
}

## -----------------------------------------------------------------------------
## S3 methods of sensFun
## -----------------------------------------------------------------------------

summary.sensFun <- function(object, vars=FALSE, ...) {

  pp       <- attributes(object)$pars
  parscale <-  attributes(object)$parscale
  Sens <- as.matrix(object[, -(1:2)])       # unwanted type conversion when only one parameter is analysed : matrix --> vector
  nout <- nrow(Sens)
  if (vars) { # summaries per variable
    Vars <- object[, 2]
    out <- data.frame(
      L1   = unlist(aggregate(abs(Sens), by = list(Vars), FUN = mean)[, -1]),
      L2   = unlist(aggregate(Sens*Sens, by = list(Vars), FUN = sum)[, -1]),
      Mean = unlist(aggregate(Sens, by = list(Vars), FUN = mean)[, -1]),
      Min  = unlist(aggregate(Sens, by = list(Vars), FUN = min)[, -1]),
      Max  = unlist(aggregate(Sens, by = list(Vars), FUN = max)[, -1]),
      N    = unlist(aggregate(Sens, by = list(Vars), FUN = length)[, -1])
    )
    out$L2 <- sqrt(out$L2/out$N)
    out$var <- levels(Vars)                                                     # changed 31-01-2012 (Tom): used to be unique(vars) => problem: aggregate orders according to levels, not first appearance (unique does)
    np <- length(pp)
    nv <- nlevels(Vars)                                                         # follows from levels(Vars)
    out <- data.frame(cbind(value = rep(pp, times=rep(nv,np)),
                      scale=rep(parscale, times=rep(nv,np)),out))
  } else {  # global summaries
    L1   <- colMeans(abs(Sens))
    L2   <- sqrt(colSums(Sens*Sens)) / nout
    Mean <- colMeans(Sens)
    Min  <- apply(Sens, 2, min)
    Max  <- apply(Sens, 2, max)
    N    <- apply(Sens, 2, length)
    out  <- data.frame(cbind(value = pp, scale = parscale, L1, L2, Mean, Min, Max, N))
    rownames(out) <- names(pp)
  }
  class(out) <- c("summary.sensFun", "data.frame")
  return(out)
}

## -----------------------------------------------------------------------------

print.summary.sensFun<-function(x, ...)
  print(format(x, digits = 2))

## -----------------------------------------------------------------------------

pairs.sensFun <- function (x, which = NULL, ...) {

  dots <- list(...)
  if(is.null(dots$pch)) dots$pch <- 16

  if (!is.null(which)) {
    nx  <-attr(x, "nx")
    var <-attr(x, "var")
    TYP <-attr(x, "Type")

    Select <- selectvar(which, var, Nall = TRUE)
    Nr <- length(var)
    ii <- NULL
    ij <- NULL
    if (TYP == 1)
      for (i in Select) {
        In <- (((i-1) * nx+1)):(i*nx)
        ij <- c(ij, rep(i, length(In)))
        ii <- c(ii, In)
      }
    else
      for (i in Select) {
        In <- (nx[i] + 1):nx[i + 1]
        ij <- c(ij,rep(i,length(In)))
        ii <- c(ii,In)
      }
    if(!is.null(dots$bg)) dots$bg <- dots$bg[ij]
    if(is.null(dots$col)) dots$col <- (1:Nr)[ij]
    else dots$col <- dots$col[ij]
  } else {
    ii <- 1:nrow(x)
    ij <- rep(1, nrow(x))
  }

  if (colnames(x)[1] == "x" && colnames(x)[2] == "var")
    X <- x[ii, -(1:2)]
  else
    X <- x[ii,]

  dots$diag.panel  <- if(is.null(dots$diag.panel)) NULL else dots$diag.panel
  dots$lower.panel <- if(is.null(dots$lower.panel)) panel.cor else dots$lower.panel
  dots$gap <- if(is.null(dots$gap)) 0 else dots$gap
  do.call("pairs", c(alist(as.matrix(X)), dots))

}

## -----------------------------------------------------------------------------

plot.sensFun<- function(x, which = NULL, legpos = "topleft", ask = NULL, ...) {
  nx    <- attr(x,"nx")
  xname <- attr(x,"xname")
  var   <- attr(x,"var")
  TYP   <- attr(x,"Type")

  dots   <- list(...)
  nmdots <- names(dots)

  nc <- ncol(x) - 2

  Main      <- is.null(dots$main)
  dots$ylab <- if(is.null(dots$ylab)) "sensitivity" else dots$ylab
  dots$type <- if(is.null(dots$type)) "l" else dots$type
  dots$col  <- if(is.null(dots$col)) 1:nc else dots$col
  if (dots$type %in% c("p","b")) 
    dots$pch  <- if(is.null(dots$pch)) 15 else dots$pch
  else
    dots$pch  <- NULL
    
  dots$lwd  <- if(is.null(dots$lwd)) 1 else dots$lwd
  dots$xlab <- if(is.null(dots$xlab)) xname else dots$xlab
  Allvars   <- FALSE
  Ylim      <- is.null(dots$ylim)

  ## Find selected variables
  Select <- selectvar(which, var, Nall = TRUE)

  ## Set par mfrow and ask.
  if (! is.null(which)) {
    ask <- setplotpar(nmdots, dots, length(Select), ask)

  } else {
    dots$ylim <- if(is.null(dots$ylim)) range(x[,-(1:2)],na.rm=T)
    Ylim      <- FALSE
    Allvars   <- TRUE
    if (is.null(ask))
      ask <- prod(par("mfrow")) < length(which) && dev.interactive()
  }

  ## interactively wait if there are remaining figures
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  Lty <- is.null(dots$lty)
  st <- 1
  if(Lty) dots$lty <- st

  ## xlim
  if (is.null(dots$xlim)) {
    is <- NULL
    for (i in Select){
      if (TYP == 1)
        ii <- ((i-1) * nx + 1):(i*nx)
      else
        ii <- (nx[i] + 1):nx[i + 1]       # sensFun calculates total amount of observations per variable in nx => plot now assumes they are grouped together,
                                          # which is not the case when using modcost with additional costs => sens object contains appended matrices for each
                                          # modcost calculation (residuals are also appended in one matrix in modcost function)
      is <- c(is, ii)
    }
    dots$xlim <- range(x[is, 1],na.rm=T)
  }
  for (i in Select){
    if (TYP == 1)
      ii <- ((i - 1)*nx + 1):(i*nx)       # changed (addition of +1; by Karline & Tom 21/10/2011)
    else
      ii <- (nx[i] + 1):nx[i + 1]
    if (Main)
      if (! Allvars) dots$main <- var[i] else dots$main <- "All variables"

    sens<- x[ii,]
    if (Ylim)  dots$ylim <- range(sens[-(1:2)],na.rm=T)

    if (st==1)
      do.call("matplot",c(alist(sens$x, as.matrix( sens[, -(1:2)])), dots))
    else
      do.call("matlines",c(alist(sens$x, as.matrix( sens[, -(1:2)])), dots))
    if (Allvars) st <- st + 1
  }
  if (! is.na(legpos))
    legend(legpos, names(x)[-(1:2)], col = dots$col, lty = dots$lty,   #  changed Tom 05/12/2011:  used to read   "legend(legpos, names(x[,-(1:2)]), col = dots$col, lty = dots$lty,"
           lwd = dots$lwd, pch = dots$pch)                             #                                                          _________________
                                                                       #  type conversion of data.frame to vector in case of only 1 analysed parameter invalidates this code
}

## -----------------------------------------------------------------------------

plot.summary.sensFun<- function(x, which = 1:nrow(x), ...) {
  dots   <- list(...)
  nmdots <- names(dots)

  mf <- par (mfrow = c(2, 3))
  Names <- names(x)
  X <- as.matrix(x[, 1:8])

  ii <- selectvar(which, rownames(x),Nall = TRUE)
  setnames <- is.null(dots$main)
  for (i in 3:7)  {
    dots$main <- if(setnames) Names[i] else dots$main
    dots$pch <- if(is.null(dots$pch)) 16 else dots$pch
    dots$col <- if(is.null(dots$col)) "black" else dots$col

    do.call("dotchart", c(alist(X[ii, i]), dots))
    abline(v = 0, lty = 3)
  }

  par (mfrow = mf)
}
