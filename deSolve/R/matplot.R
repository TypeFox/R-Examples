## =============================================================================
##  matplot methods  - it is not an S3 generic...
## =============================================================================


#matplot <- function (x, ...) UseMethod("matplot")
#matplot.default <- function (x, ...) {
#if ("deSolve" %in% class (x))
#  matplot.deSolve(x,...)
#else
#  graphics::matplot(x,...)
#  #NextMethod()
#}

matplot.deSolve <- function(x, ..., select = NULL, which = select,
    obs = NULL, obspar = list(), subset = NULL,
    legend = list(x = "topright")) {      # legend can be a list

    t       <- 1     # column with independent variable "times"

    # Set the observed data
    obs <- SetData(obs)

    # variables to be plotted and their position in "x"
    varnames <- colnames(x)
    xWhich <- NULL
    lW <- length(which)

    WhichVar <- function(xWhich, obs, varnames) {

     if (is.null(xWhich) & is.null(obs$dat))  # All variables plotted
      Which <- 2 : length(varnames)

     else if (is.null(xWhich)) {     # All common variables in x and obs plotted
       Which <- which(varnames %in% obs$name)
       Which <- Which [Which > 1]
     } else if (is.character(xWhich)) {
       Which <- which(varnames %in% xWhich)
       if (length(Which) != length(xWhich))
         stop ("unknown variable", paste(xWhich, collapse = ","))
     }
     else
       Which <- xWhich + 1
     return(Which)
    }

    if (lW & is.list(which))
      xWhich <- lapply(which, FUN = function (x) WhichVar(x, obs, varnames))
    else if (lW)
      xWhich <- list(WhichVar(which, obs, varnames))
    else
      xWhich <- list(2:length(varnames))

    vn  <- lapply(xWhich, FUN = function(x) paste(varnames[x], collapse = ","))
    vn2 <- unlist(lapply(xWhich, FUN = function(x) paste(varnames[x])))

    np <- length(xWhich)            # number of y-axes
    nx <- length(unlist(xWhich))    # number of y-variables

    # add Position of variables to be plotted in "obs"
    obs <- updateObs2 (obs, varnames, unlist(xWhich))

    # The ellipsis
    ldots  <- list(...)
    Dots   <- splitdots(ldots, varnames)

    if (Dots$nother > 1)                                             
      stop ("can plot only one deSolve output object at a time with matplot")

    Dotmain <- setdots(Dots$main, np)

    # these are different from the default
    Dotmain$xlab <- expanddots(ldots$xlab, varnames[t]                , np)
    Dotmain$ylab <- expanddots(ldots$ylab, vn                         , np)
    Dotmain$main <- expanddots(ldots$main, as.character(substitute(x)), np)

    # ylim and xlim can be lists and are at least two values
    yylim  <- expanddotslist(ldots$ylim, np)
    xxlim  <- expanddotslist(ldots$xlim, np)

    Dotpoints <- setdots(Dots$points, nx)   # expand all dots to nx values

    # these are different from default
    Dotpoints$type <- expanddots(ldots$type, "l", nx)
    Dotpoints$lty  <- expanddots(ldots$lty, 1:nx, nx)
    Dotpoints$pch  <- expanddots(ldots$pch, 1:nx, nx)
    Dotpoints$col  <- expanddots(ldots$col, 1:nx, nx)
    Dotpoints$bg   <- expanddots(ldots$bg,  1:nx, nx)

    if (! is.null(obs)) {

      ii <- which(unlist(xWhich) %in% unlist(obs$Which))
      ii <- ii[! is.na(ii)]
      if (is.null(obs$par))
       obs$par <- list()
      else
       obs$par <- lapply(obspar, repdots, obs$length)

      if (is.null(obs$par$pch))
        obs$par$pch <- Dotpoints$pch[ii]
      if (is.null(obs$par$cex))
        obs$par$cex <- Dotpoints$cex[ii]
      if (is.null(obs$par$col))
        obs$par$col <- Dotpoints$col[ii]
      if (is.null(obs$par$bg))
        obs$par$bg  <- Dotpoints$bg[ii]
    }
    if (!missing(subset)){
      e <- substitute(subset)
      r <- eval(e, as.data.frame(x), parent.frame())
      if (is.numeric(r)) {
        isub <- r
      } else {
        if (!is.logical(r))
          stop("'subset' must evaluate to logical or be a vector with integers")
        isub <- r & !is.na(r)
      }
    } else {
      isub <- TRUE
    }

    # LOOP for each (set of) output variables (and y-axes)
    if (np > 1)
      par(mar = c(5.1, 4.1, 4.1, 2.1) + c(0, (np-1)*4, 0, 0))

    ii <- 1
    for (ip in 1 : np) {

      ix  <- xWhich[[ip]]     # position of y-variables in 'x'
      iL  <- length(ix)
      iip <- ii:(ii+iL-1)     # for dotpoints
      ii  <- ii + iL
      io <- obs$Which[iip]

      # plotting parameters for matplot and axes
      dotmain   <- extractdots(Dotmain, ip)
      if (is.null(dotmain$axes))       dotmain$axes <- FALSE
      if (is.null(dotmain$frame.plot)) dotmain$frame.plot <- TRUE

      dotpoints <- extractdots(Dotpoints, iip)    # for all variables

      Xlog <- Ylog <- FALSE
      if (! is.null(dotmain$log)) {
        Ylog  <- length(grep("y",dotmain$log))
        Xlog  <- length(grep("x",dotmain$log))
      }

      SetRangeMat <- function(lim, x, isub, ix, obs, io, Log) {

       if ( is.null (lim)) {
        yrange <- Range(NULL, as.vector(x[isub, ix]), Log)
        if (! is.na(io[1])) yrange <- Range(yrange, as.vector(obs$dat[,io]), Log)
       } else
         yrange  <- lim

      return(yrange)
      }

      dotmain$ylim <- SetRangeMat(yylim[[ip]], x, isub, ix, obs, io, Ylog)
      dotmain$xlim <- SetRangeMat(xxlim[[ip]], x, isub,  t, obs, io, Xlog)

      Ylab <- dotmain$ylab
      dotmain$ylab <- ""
      if (ip > 1) {
        par(new = TRUE)
        dotmain$xlab <- dotmain$main <- ""
      }

      do.call("matplot", c(alist(x[isub, t], x[isub, ix]), dotmain, dotpoints))
      if (ip == 1)
        axis(1, cex = dotmain$cex.axis)

      cex <- ifelse (is.null(dotmain$cex.lab), 0.9, 0.9*dotmain$cex.lab)
      bL <- 4*(ip-1)
      axis(side = 2, line = bL, cex = dotmain$cex.axis)
      mtext(side = 2, line = bL+2, Ylab, cex = cex)
      if (! is.na(io[1]))
        for (j in 1: length(io)) {
          i <- which (obs$Which == io[j])
          if (length (i.obs <- obs$pos[i, 1]:obs$pos[i, 2]) > 0)
            do.call("points", c(alist(obs$dat[i.obs, 1], obs$dat[i.obs, io[j]]),
                  extractdots(obs$par, j) ))
        }
    }

    if (is.null(legend))
      legend <- list(x = "topright")

    if (is.list(legend)){  # can also be FALSE
      if (length(legend$legend))
        L <- legend$legend
      else
        L <- vn2
      legend$legend <- NULL
      if (is.null(legend$x))
        legend$x <- "topright"
      lty <- Dotpoints$lty
      pch <- Dotpoints$pch
      lty[Dotpoints$type == "p"] <- NA
      pch[Dotpoints$type == "l"] <- NA
      do.call ("legend", c(legend, alist(lty = lty, lwd = Dotpoints$lwd,
               pch =pch, col = Dotpoints$col, pt.bg =Dotpoints$bg,
               legend = L)))
    }
  }

### ============================================================================
### plotting 1-D variables as line plot, one for each time
### ============================================================================

matplot.1D <- function (x, select= NULL, which = select, ask = NULL,
                        obs = NULL, obspar = list(), grid = NULL,
                        xyswap = FALSE, vertical = FALSE, subset = NULL, ...) {

  ## Check settings of x
  att    <- attributes(x)
  nspec  <- att$nspec
  dimens <- att$dimens

  proddim <- prod(dimens)
  if (length(dimens) != 1)
    stop ("matplot.1D only works for models solved with 'ode.1D'")

  if ((ncol(x)- nspec*proddim) < 1)
    stop("ncol of 'x' should be > 'nspec' * dimens if x is a vector")

  # Set the observed data
  obs <- SetData(obs)

  # 1-D variable names
  varnames <-  if (! is.null(att$ynames))
    att$ynames else 1:nspec
  if (! is.null(att$lengthvar))
    varnames <- c(varnames, names(att$lengthvar)[-1])

  # variables to be plotted, common between obs and x
  Which <- WhichVarObs(which, obs, nspec, varnames, remove1st = FALSE)

  np <- length(Which)

  # Position of variables to be plotted in "x"
  Select <- select1dvar(Which, varnames, att)  # also start and end position
  xWhich  <- Select$Which

  # add Position of variables to be plotted in "obs"
  obs <- updateObs (obs, varnames, xWhich)
  obs$par <- lapply(obspar, repdots, obs$length)

  # the ellipsis
  ldots <- list(...)

  # number of figures in a row and interactively wait if remaining figures
  ask <- setplotpar(ldots, np, ask)
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  Dots <- splitdots(ldots, varnames)

  nother <- Dots$nother

  Dotpoints <- Dots$points
  Dotmain <- setdots(Dots$main, np) # expand all dots to np values (no defaults)

  # These are different from defaulst
  Dotmain$xlab <- expanddots(ldots$xlab,  "x", np)
  Dotmain$ylab <- expanddots(ldots$ylab,  "", np)
  Dotmain$main <- expanddots(ldots$main,  varnames[xWhich], np)

  # xlim and ylim are special:
  xxlim <- expanddotslist(ldots$xlim, np)
  yylim <- expanddotslist(ldots$ylim, np)

  xyswap <- rep(xyswap, length = np)
  vertical <- rep(vertical, length = np)

  if (!missing(subset)){

  e <- substitute(subset)
  r <- eval(e, as.data.frame(x), parent.frame())
  if (is.numeric(r)) {
     isub <- r
  } else {
    if (!is.logical(r))
      stop("'subset' must evaluate to logical or be a vector with integers")
    isub <- r & !is.na(r)
  }
  } else isub <- 1:nrow(x)

  grid <- expanddotslist(grid, np)

  for (ip in 1:np) {

    istart <- Select$istart[ip]
    istop  <- Select$istop[ip]
    io <- obs$Which[ip]

    out <- t(x[ isub, istart:istop])
    if (length (isub) > 1 & sum (isub) == 1)
      out <- matrix (out)

    Grid <- grid[[ip]]
    if (is.null(Grid))
      Grid <- 1:nrow(out)

    dotmain      <- extractdots(Dotmain, ip)

    Xlog <- Ylog <- FALSE
    if (! is.null(dotmain$log)) {
       Ylog  <- length(grep("y", dotmain$log))
       Xlog  <- length(grep("x", dotmain$log))
    }
    if (vertical[ip])  {  # overrules other settings; vertical profiles
      xyswap[ip] <- TRUE
      dotmain$axes <- FALSE
      dotmain$xlab <- ""
      dotmain$xaxs <- "i"
      dotmain$yaxs <- "i"
    }

    if (! xyswap[ip]) {
      if (! is.null(xxlim[[ip]]))
        dotmain$xlim <- xxlim[[ip]]
      dotmain$ylim <- SetRange(yylim[[ip]], x, NULL, isub, istart:istop, obs, io, Ylog)
    } else {
      if (! is.null(yylim[[ip]]))
        dotmain$ylim <- yylim[[ip]]
      dotmain$xlim <- SetRange(xxlim[[ip]], x, NULL, isub, istart:istop, obs, io, Xlog)
      if (is.null(yylim[[ip]]) & xyswap[ip])
        dotmain$ylim <- rev(range(Grid))    # y-axis
    }

    if (! xyswap[ip]) {
      do.call("matplot", c(alist(Grid, out), dotmain, Dotpoints))

      if (! is.na(io))
        plotObs(obs, io)
    } else {
      if (is.null(dotmain$xlab[ip]) | is.null(dotmain$ylab[ip])) {
         dotmain$ylab <- dotmain$xlab[ip]
         dotmain$xlab <- dotmain$ylab[ip]
      }

      do.call("matplot", c(alist(out, Grid), dotmain, Dotpoints))

      if (vertical[ip])
         DrawVerticalAxis(dotmain, min(out))
      if (! is.na(io))
         plotObs(obs, io, xyswap = TRUE)
    }
  }
}

## =============================================================================
## S3/S4 compatibility
## =============================================================================


##  make matplot an S4 method and then extend generic for class deSolve
##  but note that matplot.1D is not (yet) a generic, because .1D is just an
##  alternative way of plotting and not a well defined class
setGeneric("matplot", function(x, ...) graphics::matplot(x, ...))
setOldClass("deSolve")

setMethod("matplot", list(x = "deSolve"), matplot.deSolve)