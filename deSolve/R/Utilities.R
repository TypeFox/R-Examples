### ============================================================================
### ============================================================================
### S3 methods
### karline+Thomas: from version 1.9, also possible to plot multiple
###                 outputs and to add observations.
### ============================================================================
### ============================================================================

### ============================================================================
### first some common functions
### ============================================================================

## =============================================================================
## Update range, taking into account neg values for log transformed values
## =============================================================================

Range <- function(Range, x, log) {
   if ((log) & (!is.null(x)))
      x[x <= 0] <- min(x[x > 0])  # remove zeros
   return(range(Range, x, na.rm = TRUE) )
}

## =============================================================================
## Checking and expanding arguments in dots (...) with default
## =============================================================================

expanddots <- function (dots, default, n) {
  dots <- if (is.null(dots)) default else dots
  rep(dots, length.out = n)
}

# lists: e.g. xlim and ylim....
expanddotslist <- function (dots, n) {
  if (is.null(dots)) return(dots)
  dd <- if (!is.list(dots )) list(dots) else dots
  rep(dd, length.out = n)
}

## =============================================================================
## Expanding arguments in dots  (...)
## =============================================================================

repdots <- function(dots, n)
  if (is.function(dots)) dots else rep(dots, length.out = n)

setdots <- function(dots, n) lapply(dots, repdots, n)

## =============================================================================
## Extracting element 'index' from dots  (...)
## =============================================================================

extractdots <- function(dots, index) {
  ret <- lapply(dots, "[", index)
  ret <- lapply(ret, unlist) # flatten list
  return(ret)
}

## =============================================================================
## Merge two observed data files; assumed that first column = 'x' and ignored
## =============================================================================

# from 3-columned format (what, where, value) to wide format...
convert2wide <- function(Data) {
    cnames   <- as.character(unique(Data[,1]))

    MAT      <- Data[Data[,1] == cnames[1], 2:3]
    colnames.MAT <- c("x", cnames[1])

    for ( ivar in cnames[-1]) {
      sel <- Data[Data[,1] == ivar, 2:3]
      nt  <- cbind(sel[,1],
                   matrix(nrow = nrow(sel), ncol = ncol(MAT)-1, data = NA),
                   sel[,2])
      MAT <- cbind(MAT, NA)
      colnames(nt) <- colnames(MAT)
      MAT <- rbind(MAT, nt)
      colnames.MAT <- c(colnames.MAT, ivar)
    }
  colnames(MAT) <- colnames.MAT
  return(MAT)
}

# merge two observed data sets in one

mergeObs <- function(obs, Newobs) {

  if (! class(Newobs) %in% c("data.frame", "matrix"))
    stop ("the elements in 'obs' should be either a 'data.frame' or a 'matrix'")

  if (is.character(Newobs[, 1]) | is.factor(Newobs[, 1]))
    Newobs <- convert2wide(Newobs)

  obsname <- colnames(obs)

  ## check if some observed variables in NewObs are already in obs
  newname <- colnames(Newobs)[-1]    # 1st column = x-var and ignored
  ii <- which (newname %in% obsname)
  if (length(ii) > 0)
    obsname <- c(obsname, newname[-ii] )
  else
    obsname <- c(obsname, newname)

  ## padding with NA of the two datasets
  O1 <- matrix(nrow = nrow(Newobs), ncol = ncol(obs), data = NA)
  O1[ ,1] <- Newobs[, 1]

  for (j in ii) {   # observed data in common are put in correct position
    jj <- which (obsname == newname[j])
    O1[,jj] <- Newobs[, j+1]
  }
  O1 <- cbind(O1, Newobs[, -c(1, ii+1)] )
  colnames(O1) <- obsname

  nnewcol <- ncol(Newobs)-1 - length (ii)  # number of new columns
  if (nnewcol > 0) {
    O2 <- matrix(nrow = nrow(obs), ncol = nnewcol, data = NA)
    O2 <- cbind(obs, O2)
    colnames(O2) <- obsname
  } else O2 <- obs

  obs <- rbind(O2, O1)
  return(obs)
}

## =============================================================================
## Set the mfrow parameters and whether to "ask" for opening a new device
## =============================================================================

setplotpar <- function(ldots, nv, ask) {
  nmdots <- names(ldots)
  # nv = number of variables to plot
  if (!any(match(nmdots, c("mfrow", "mfcol"), nomatch = 0))) {
     nc <- min(ceiling(sqrt(nv)), 3)
     nr <- min(ceiling(nv/nc), 3)
     mfrow <- c(nr, nc)
  } else if ("mfcol" %in% nmdots)
     mfrow <- rev(ldots$mfcol)
  else mfrow <- ldots$mfrow

  if (! is.null(mfrow))  mf <- par(mfrow = mfrow)

  ## interactively wait if there are remaining figures
  if (is.null(ask))
    ask <- prod(par("mfrow")) < nv && dev.interactive()

  return(ask)
}

## =============================================================================
## find a variable
## =============================================================================

selectvar <- function (Which, var, NAallowed = FALSE) {
  if (!is.numeric(Which)) {
    ln <- length(Which)

    ## the loop is necessary so as to keep ordering...
    Select <- NULL
    for ( i in 1:ln) {
      ss <- which(Which[i] == var)
      if (length(ss) ==0 & ! NAallowed)
        stop("variable ", Which[i], " not in variable names")
      else if (length(ss) == 0)
        Select <- c(Select, NA)
      else
        Select <- c(Select, ss)
    }

  } else {
    Select <- Which + 1  # "Select" now refers to the column number
    if (max(Select) > length(var))
        stop("index in 'which' too large: ", max(Select)-1)
    if (min(Select) < 1)
        stop("index in 'which' should be > 0")
  }
  return(Select)
}

### ============================================================================
### print a deSolve object
### ============================================================================

print.deSolve <- function(x, ...)
  print(as.data.frame(x), ...)

### ============================================================================
### Create a histogram for a list of variables
### ============================================================================

hist.deSolve <- function (x, select = 1:(ncol(x)-1), which = select, ask = NULL,
                          subset = NULL, ...) {

  t        <- 1     # column with independent variable ("times")
  varnames <- colnames(x)
  Which    <- selectvar(which, varnames)

  np     <- length(Which)
  ldots  <- list(...)

  ## Set par mfrow and ask
  ask <- setplotpar(ldots, np, ask)
  if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
  }

  ## expand all dots to np values (no defaults)
  Dotmain  <- setdots(ldots, np)

  ## different from default settings
  Dotmain$main <- expanddots (ldots$main, varnames[Which], np)
  Dotmain$xlab <- expanddots (ldots$xlab, varnames[t],     np)
  # Dotmain$xlab <- expanddots (ldots$xlab, ""        ,     np)

  ## xlim and ylim are special: they are vectors or lists
  xxlim <- expanddotslist(ldots$xlim, np)
  yylim <- expanddotslist(ldots$ylim, np)

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
  } else isub <- TRUE

  ## plotting
  for (ip in 1:np) {
      ix <- Which[ip]
      dotmain <- extractdots(Dotmain, ip)
      if (! is.null(xxlim[[ip]])) dotmain$xlim <- xxlim[[ip]]
      if (! is.null(yylim[[ip]])) dotmain$ylim <- yylim[[ip]]
      do.call("hist", c(alist(x[isub, ix]), dotmain))
  }
}

### ============================================================================
### Image, filled.contour and persp plots
### ============================================================================

image.deSolve <- function (x, select = NULL, which = select, ask = NULL,
    add.contour = FALSE, grid = NULL, method = "image",
    legend = FALSE, subset = NULL, ...) {

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
  } else isub <- TRUE

  dimens <- attributes(x)$dimens
  if (is.null(dimens))
    stop("cannot make an image from deSolve output which is 0-dimensional")
  else if (length(dimens) ==1)  # 1-D
    plot.ode1D(x, which, ask, add.contour, grid, method=method,
        legend = legend, isub = isub, ...)
  else if (length(dimens) ==2)  # 2-D
    plot.ode2D(x, which, ask, add.contour, grid, method=method,
        legend = legend, isub = isub, ...)
  else
    stop("cannot make an image from deSolve output with more than 2 dimensions")
}

### ============================================================================
### Plot utilities for the S3 plot method, 0-D, 1-D, 2-D
### ============================================================================

## ============================================================================
## Observations cleanup
## ============================================================================

SetData <- function(obs) { ## check observed data
  nobs <- 0
  obs.pos <- NULL
  obsname <- NULL
  if (! is.null(obs)) {
    if (!is.data.frame(obs) & is.list(obs)) { # a list with different data sets
      Obs <- obs
      obs <- Obs[[1]]
      obs.pos <- matrix(nrow = 1, c(1, nrow(obs)))
      if (! class(obs) %in% c("data.frame", "matrix"))
        stop ("'obs' should be either a 'data.frame' or a 'matrix'")
      if (length(Obs) > 1)
        for ( i in 2 : length(Obs)) {
          obs <- mergeObs(obs, Obs[[i]])
          obs.pos <- rbind(obs.pos, c(obs.pos[nrow(obs.pos), 2] +1, nrow(obs)))
        }
      obsname <- colnames(obs)
    } else {                                 # a data.frame or matrix
      if (is.character(obs[, 1]) | is.factor(obs[, 1]))   # long format - convert
        obs <- convert2wide(obs)
      obsname <- colnames(obs)
      if (! class(obs) %in% c("data.frame", "matrix"))
        stop ("'obs' should be either a 'data.frame' or a 'matrix'")
      obs.pos <- matrix(nrow = 1, c(1, nrow(obs)))
    }
    DD <- duplicated(obsname)
    if (sum(DD) > 0)
      obs <- mergeObs(obs[,!DD], cbind(obs[, 1], obs[, DD]))
    nobs <- nrow(obs.pos)
  }
  return(list(dat = obs, pos = obs.pos, name = obsname, length = nobs))
}

## ============================================================================
## create several lists: x2:   other deSolve objects,
##                       dotmain, dotpoints: remaining (plotting) parameters
## ============================================================================

splitdots <- function(ldots, varnames){
  x2     <- list()
  dots   <- list()
  nd     <- 0
  nother <- 0
  ndots <- names(ldots)

  if (length(ldots) > 0)
    for ( i in 1:length(ldots))
      if ("deSolve" %in% class(ldots[[i]])) { # a deSolve object
        x2[[nother <- nother + 1]] <- ldots[[i]]
        names(x2)[nother] <- ndots[i]
        # a list of deSolve objects
      } else if (is.list(ldots[[i]]) & "deSolve" %in% class(ldots[[i]][[1]])) {
        for (j in 1:length(ldots[[i]])) {
          x2[[nother <- nother+1]] <- ldots[[i]][[j]]
          names(x2)[nother] <- names(ldots[[i]])[[j]]
        }
      } else if (! is.null(ldots[[i]])) {  # a graphical parameter
        dots[[nd <- nd+1]] <- ldots[[i]]
        names(dots)[nd] <- ndots[i]
      }

  nmdots <- names(dots)

  # check compatibility of all deSolve objects
  if (nother > 0) {
    for ( i in 1:nother) {
      if (min(colnames(x2[[i]]) == varnames) == 0)
        stop("'x' is not compatible with other deSolve objects - colnames not the same")
    }
  }

  # plotting parameters : split in plot parameters and point parameters
  plotnames <- c("xlab", "ylab", "xlim", "ylim", "main", "sub", "log", "asp",
                 "ann", "axes", "frame.plot", "panel.first", "panel.last",
                 "cex.lab", "cex.axis", "cex.main")

  # plot.default parameters
  ii <- names(dots) %in% plotnames
  dotmain <- dots[ii]

  # point parameters
  ip <- !names(dots) %in% plotnames
  dotpoints <- dots[ip]
  list(points = dotpoints, main = dotmain, nother = nother, x2 = x2)
}

## =============================================================================
## Which variable in common between observed and selected variables
## =============================================================================

WhichVarObs <- function(Which, obs, nvar, varnames, remove1st = TRUE) {

  if (is.null(Which) & is.null(obs$dat))  # All variables plotted
    Which <- 1 : nvar

  else if (is.null(Which)) {     # All common variables in x and obs plotted
    Which <- which(varnames %in% obs$name)
    if (remove1st) Which <- Which[Which != 1]  # remove first element (x-value)
      Which <- varnames[Which]                 # names rather than numbers
  }
  return(Which)
}

## =============================================================================
## Update Obs with position of observed variable in x
## =============================================================================

updateObs <- function (obs, varnames, xWhich) {
  if (obs$length > 0 ) {
    obs$Which <- selectvar(varnames[xWhich], obs$name, NAallowed = TRUE)
    obs$Which [ obs$Which > ncol(obs$dat)] <- NA
#    if (nrow(obs$pos) != length(obs$Which))
#      obs$pos <- matrix(nrow = length(obs$Which), ncol = ncol(obs$pos),
#        byrow = TRUE, data =obs$pos[1,])
  } else
    obs$Which <- rep(NA, length(xWhich))
  return(obs)
}
updateObs2 <- function (obs, varnames, xWhich) {
  if (obs$length > 0 ) {
    obs$Which <- selectvar(varnames[xWhich], obs$name, NAallowed = TRUE)
    obs$Which [ obs$Which > ncol(obs$dat)] <- NA
    if (nrow(obs$pos) != length(obs$Which))
      obs$pos <- matrix(nrow = length(obs$Which), ncol = ncol(obs$pos),
        byrow = TRUE, data =obs$pos[1,])
  } else
    obs$Which <- rep(NA, length(xWhich))
  return(obs)
}

## =============================================================================
## Set range of a plot, depending on deSolve object and data...
## =============================================================================

SetRange <- function(lim, x, x2, isub, ix, obs, io, Log) {

  nother <- length (x2)
  if ( is.null (lim)) {
    yrange <- Range(NULL, x[isub, ix], Log)
    if (nother>0)
      for (j in 1:nother)
        yrange <- Range(yrange, x2[[j]][isub,ix], Log)
      if (! is.na(io)) yrange <- Range(yrange, obs$dat[,io], Log)
  } else
     yrange  <- lim

  return(yrange)
}

## =============================================================================
## Add observed data to a plot
## =============================================================================

plotObs <- function (obs, io, xyswap = FALSE) {
  oLength <- min(nrow(obs$pos), obs$length)
  if (! xyswap) {
    for (j in 1: oLength) {
      i.obs <- obs$pos[j, 1] : obs$pos[j, 2]
      if (length (i.obs) > 0)
        do.call("points", c(alist(obs$dat[i.obs, 1], obs$dat[i.obs, io]),
                extractdots(obs$par, j) ))
       }
  } else {
    for (j in 1: oLength)
      if (length (i.obs <- obs$pos[j, 1]:obs$pos[j, 2]) > 0)
        do.call("points", c(alist(obs$dat[i.obs, io], obs$dat[i.obs, 1]),
                extractdots(obs$par, j) ))
  }
}

### ============================================================================
### Plotting 0-D variables
### ============================================================================


plot.deSolve <- function (x, ..., select = NULL, which = select, ask = NULL,
                          obs = NULL, obspar = list(), subset = NULL) {

  t       <- 1     # column with independent variable "times"

  # Set the observed data
  obs <- SetData(obs)

  # variables to be plotted
  varnames <- colnames(x)
  Which    <- WhichVarObs(which, obs, ncol(x) - 1, varnames)

  # Position of variables to be plotted in "x"
  xWhich  <- selectvar(Which, varnames)
  np      <- length(xWhich)

  # Position of variables in "obs" (NA = not observed)
  obs <- updateObs(obs, varnames, xWhich)
  obs$par <- lapply(obspar, repdots, obs$length)

  # The ellipsis
  ldots   <- list(...)

  # number of figures in a row and interactively wait if remaining figures
  ask <- setplotpar(ldots, np, ask)
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  Dots <- splitdots(ldots, varnames)

  nother <- Dots$nother
  x2     <- Dots$x2
  nx     <- nother + 1 # total number of deSolve objects to be plotted

  Dotmain <- setdots(Dots$main, np)  # expand to np for each plot

  # these are different from the default
  Dotmain$xlab <- expanddots(ldots$xlab, varnames[t]     , np)
  Dotmain$ylab <- expanddots(ldots$ylab, ""              , np)
  Dotmain$main <- expanddots(ldots$main, varnames[xWhich], np)

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

  # LOOP for each output variable (plot)

  for (ip in 1 : np) {

    ix <- xWhich[ip]      # position of variable in 'x'
    io <- obs$Which[ip]   # position of variable in 'obs'

    # plotting parameters for deSolve output 1 (opens a plot)
    dotmain   <- extractdots(Dotmain, ip)
    dotpoints <- extractdots(Dotpoints, 1)  # 1st dotpoints

    Xlog <- Ylog <- FALSE
    if (! is.null(dotmain$log)) {
      Ylog  <- length(grep("y",dotmain$log))
      Xlog  <- length(grep("x",dotmain$log))
    }

    dotmain$ylim <- SetRange(yylim[[ip]], x, x2, isub, ix, obs, io, Ylog)
    dotmain$xlim <- SetRange(xxlim[[ip]], x, x2, isub,  t, obs,  1, Xlog)

    # first deSolve object plotted (new plot created)
    do.call("plot", c(alist(x[isub, t], x[isub, ix]), dotmain, dotpoints))

    if (nother > 0)        # if other deSolve outputs
      for (j in 2:nx)
        do.call("lines", c(alist(x2[[j-1]][isub, t], x2[[j-1]][isub, ix]),
                extractdots(Dotpoints, j)) )

    if (! is.na(io)) plotObs(obs, io)   # add observed variables
  }
}


## =============================================================================
##  to draw a legend
## =============================================================================

drawlegend <- function (parleg, dots) {
  Plt <- par(plt = parleg)
  par(new = TRUE)
  usr <- par("usr")
  ix <- 1
  minz <- dots$zlim[1]
  maxz <- dots$zlim[2]
  binwidth <- (maxz - minz)/64
  iy <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
  iz <- matrix(iy, nrow = 1, ncol = length(iy))

  image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "",
        ylab = "", col = dots$col)

  do.call("axis", list(side = 4, mgp = c(3, 1, 0), las = 2))

  par(plt = Plt)
  par(usr = usr)
  par(new = FALSE)
}

## =============================================================================
## to drape a color over a persp plot.
## =============================================================================

drapecol <- function (A,
          col = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
              "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100),
              NAcol = "white", Range = NULL)
{
  nr <- nrow(A)
  nc <- ncol(A)
  ncol <- length(col)

  AA <- 0.25 * (A[1:(nr - 1), 1:(nc - 1)] + A[1:(nr - 1), 2:nc] +
        A[2:nr, 1:(nc - 1)] + A[2:nr, 2:nc])
  if (is.null(Range))
    Range <- range(A, na.rm = TRUE)
  else {
    AA[AA > Range[2]] <- Range[2]
    AA[AA < Range[1]] <- Range[1]
  }
  Ar <- Range
  rn <- Ar[2] - Ar[1]
  ifelse(rn != 0, drape <- col[1 + trunc((AA - Ar[1])/rn *
        (ncol - 1))], drape <- rep(col[1], ncol))
  drape[is.na(drape)] <- NAcol
  return(drape)
}

## =============================================================================
## Finding 1-D variables
## =============================================================================

select1dvar <- function (Which, var, att) {

  if (is.null(att$map))
    proddim <- prod(att$dimens)
  else
    proddim <- sum(!is.na(att$map))

  ln   <- length(Which)
  csum <- cumsum(att$lengthvar) + 2

  if (!is.numeric(Which)) {
    # loop used to keep ordering...
    Select <- NULL
    for ( i in 1 : ln) {
      ss <- which(Which[i] == var)
      if (length(ss) == 0)
        stop("variable ", Which[i], " not in variable names")
      Select <- c(Select, ss)
    }
  } else {
    Select <- Which  # "Select now refers to the column number
    if (max(Select) > length(var))
      stop("index in 'which' too large")
    if (min(Select) < 1)
      stop("index in 'which' should be > 0")
  }

  istart <- numeric(ln)
  istop  <- numeric(ln)
  for ( i in 1 : ln) {
    if (Select[i] <= att$nspec) {
      ii <- Select[i]
      istart[i] <- (ii-1)*proddim + 2
      istop[i]  <- istart[i] + proddim - 1
    } else {
      ii <- Select[i] - att$nspec
      istart[i] <- csum[ii]
      istop[i]  <- csum[ii+1]-1
    }
    if (istart[i] == istop[i])
      stop ("variable ",Which[i], " is not a 1-D variable")

  }
  return(list(Which = Select, istart = istart, istop = istop))
}

## =============================================================================
## Finding 2-D variables
## =============================================================================

select2dvar <- function (Which, var, att) {

  if (is.null(att$map))
    proddim <- prod(att$dimens)
  else
    proddim <- sum(!is.na(att$map))
  ln   <- length(Which)
  csum <- cumsum(att$lengthvar) + 2

  if (!is.numeric(Which)) {
     # loop to keep ordering...
     Select <- NULL
     for ( i in 1 : ln) {
       ss <- which(Which[i] == var)
       if (length(ss) == 0)
         stop("variable ", Which[i], " not in variable names")
       Select <- c(Select, ss)
     }
  } else {
    Select <- Which  # "Select now refers to the column number
    if (max(Select) > length(var))
        stop("index in 'which' too large")
    if (min(Select) < 1)
        stop("index in 'which' should be > 0")
  }

  istart <- numeric(ln)
  istop  <- numeric(ln)
  dimens <- list()
  for ( i in 1 : ln) {
    if (Select[i] <= att$nspec) {  # a state variable
      ii <- Select[i]
      istart[i] <- (ii-1)*proddim + 2
      istop[i]  <- istart[i] + proddim-1
      dimens[[i]] <- att$dimens
    } else {
      ii <- Select[i] - att$nspec
      istart[i] <- csum[ii]
      istop[i]  <- csum[ii+1]-1
      ij <- which(names(att$dimvar) == var[Select[i]])
      if (length(ij) == 0)
        stop("variable ",var[Select]," is not two-dimensional")
      dimens[[i]] <- att$dimvar[[ij]]

    }
  }
  return(list(Which = Select, istart = istart, istop = istop, dim = dimens))
}

## =============================================================================
## Adding a vertical axis to a plot
## =============================================================================

DrawVerticalAxis <- function (dot, xmin) {

  if (is.null(dot$xlim))
    v <- xmin
  else
    v <- dot$xlim[1]

  abline(h = dot$ylim[2])
  abline(v = v)
  axis(side = 2)
  axis(side = 3, mgp = c(3,0.5,0))
}


### ============================================================================
### plotting 1-D variables as line plot, one for each time
### ============================================================================

plot.1D <- function (x, ... , select= NULL, which = select, ask = NULL,
                     obs = NULL, obspar = list(), grid = NULL,
                     xyswap = FALSE, delay = 0, vertical = FALSE,
                     subset = NULL) {

## Check settings of x
  att     <- attributes(x)
  nspec   <- att$nspec
  dimens  <- att$dimens
  proddim <- prod(dimens)
  if (length(dimens) != 1)
    stop ("plot.1D only works for models solved with 'ode.1D'")

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

  Select <- select1dvar(Which, varnames, att)
  xWhich <- Select$Which

  # add Position of variables to be plotted in "obs"
  obs <- updateObs (obs, varnames, xWhich)
  obs$par <- lapply(obspar, repdots, obs$length)         # karline: small bug fixed here

  # the ellipsis
  ldots <- list(...)

  ## number of figures in a row and interactively wait if remaining figures
  ask <- setplotpar(ldots, np, ask)
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  Dots <- splitdots(ldots, colnames(x))
  # for time-moving figures; number of plots should = mfrow settings
  prodx <- prod(par("mfrow"))
  if (np < prodx) eplot <- prodx - np else eplot <- 0

  nother <- Dots$nother
  x2     <- Dots$x2
  nx     <- nother + 1 # total number of deSolve objects to be plotted

  Dotmain <- setdots(Dots$main, np)  # expand to np for each plot
  Dotpoints <- setdots(Dots$points, nx)

  # These are different from defaulst
  Dotmain$xlab <- expanddots(ldots$xlab,  "x", np)
  Dotmain$ylab <- expanddots(ldots$ylab,  varnames[xWhich], np)


  # xlim and ylim are special:
  xxlim <- expanddotslist(ldots$xlim, np)
  yylim <- expanddotslist(ldots$ylim, np)

  xyswap <- rep(xyswap, length = np)
  vertical <- rep(vertical, length = np)

  grid <- expanddotslist(grid, np)

  if (!missing(subset)){
    e <- substitute(subset)
    r <- eval(e, as.data.frame(x), parent.frame())
    if (is.numeric(r)) {
      isub <- r
    } else {
      if (!is.logical(r))
        stop("'subset' must evaluate to logical or be a vector with integers")
      isub <- which(r & !is.na(r))
    }
  } else {
    isub <- 1:nrow(x)
  }
  # allow individual xlab and ylab (vectorized)
  times <- x[isub,1]
  Dotsmain <- expanddots(Dotmain$main, paste("time", times), length(times))

  for (j in isub) {
    for (ip in 1:np) {
      istart <- Select$istart[ip]
      istop  <- Select$istop[ip]
      io <- obs$Which[ip]

      out <- x[j,istart:istop]
      Grid <- grid[[ip]]
      if (is.null(Grid))
        Grid <- 1:length(out)

      dotmain      <- extractdots(Dotmain, ip)
      dotpoints    <- extractdots(Dotpoints, 1)  # 1st one

      dotmain$main <- Dotsmain[j]
      if (vertical[ip])  {  # overrules other settings; vertical profiles
        xyswap[ip] <- TRUE
        dotmain$axes <- FALSE
        dotmain$xlab <- ""
        dotmain$xaxs <- "i"
        dotmain$yaxs <- "i"
      }
    Xlog <- Ylog <- FALSE
    if (! is.null(dotmain$log)) {
      Ylog  <- length(grep("y",dotmain$log))
      Xlog  <- length(grep("x",dotmain$log))
    }

    if (! xyswap[ip]) {
      if (! is.null(xxlim[[ip]]))
        dotmain$xlim <- xxlim[[ip]]
      dotmain$ylim <- SetRange(yylim[[ip]], x, x2, isub, istart:istop, obs, io, Ylog)
    } else {
      if (! is.null(yylim[[ip]]))
        dotmain$ylim <- yylim[[ip]]
      dotmain$xlim <- SetRange(xxlim[[ip]], x, x2, isub, istart:istop, obs, io, Xlog)
      if (is.null(yylim[[ip]]) & xyswap[ip])
        dotmain$ylim <- rev(range(Grid))    # y-axis
    }

      if (! xyswap[ip]) {
        do.call("plot", c(alist(Grid, out), dotmain, dotpoints))
        if (nother > 0)        # if other deSolve outputs
          for (jj in 2:nx)
            do.call("lines", c(alist(Grid, x2[[jj-1]][j,istart:istop]),
                    extractdots(Dotpoints, jj)) )
        if (! is.na(io))
          plotObs(obs, io)

      } else {
        if (is.null(Dotmain$xlab[ip]) | is.null(Dotmain$ylab[ip])) {
          dotmain$ylab <- Dotmain$xlab[ip]
          dotmain$xlab <- Dotmain$ylab[ip]
        }

        do.call("plot", c(alist(out, Grid), dotmain, dotpoints))
        if (nother > 0)        # if other deSolve outputs
          for (jj in 2:nx)
            do.call("lines", c(alist(x2[[jj-1]][j,istart:istop], Grid),
                    extractdots(Dotpoints, jj)) )

        if (vertical[ip]) DrawVerticalAxis(dotmain,min(out))
        if (! is.na(io))
         plotObs(obs, io, xyswap = TRUE)

      }
    } # end loop ip
      if (eplot > 0)
        for (i in 1:eplot) plot(0, type ="n", axes = FALSE, xlab="", ylab="")
    if (delay > 0) Sys.sleep(0.001 * delay)
  }
}

### ============================================================================

plot.ode1D <- function (x, which, ask, add.contour, grid,
   method = "image", legend, isub = 1:nrow(x), ...) {

  # Default color scheme
  BlueRed <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

  # if x is vector, check if there are enough columns ...
  att <- attributes(x)
  nspec <- att$nspec
  dimens <- att$dimens
  proddim <- prod(dimens)

  if ((ncol(x)- nspec * proddim) < 1)
    stop("ncol of 'x' should be > 'nspec' * dimens if x is a vector")

  # variables to be plotted
  if (is.null(which))
    Which <- 1 : nspec
  else
    Which <- which
  np <- length(Which)

  varnames <-  if (! is.null(att$ynames)) att$ynames else 1:nspec

  if (! is.null(att$lengthvar))
    varnames <- c(varnames, names(att$lengthvar)[-1])

  Select <- select1dvar(Which, varnames, att)
  Which <- Select$Which

  ldots <- list(...)

  # number of figures in a row and interactively wait if remaining figures
  ask <- setplotpar(ldots, np, ask)
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  Dotmain  <- setdots(ldots, np)   # expand dots to np values (no defaults)

  # different from the default
  Dotmain$main  <- expanddots(ldots$main, varnames[Which], np)
  Dotmain$xlab  <- expanddots(ldots$xlab, "times",  np)
  Dotmain$ylab  <- expanddots(ldots$ylab, "",       np)

  # colors - different if persp, image or filled.contour

  if (method == "persp")
    dotscol <- ldots$col

  else if (method == "filled.contour")  {
    dotscolorpalette <- if (is.null(ldots$color.palette))
      BlueRed else ldots$color.palette
    dotscol <- dotscolorpalette(100)
    add.contour <- FALSE
    legend <- FALSE
  } else
    if (is.null(ldots$col))
      dotscol <- BlueRed(100) else dotscol <- ldots$col

  Addcontour <- rep(add.contour, length = np)

  # xlim, ylim and zlim are special:
  xxlim <- expanddotslist(ldots$xlim, np)
  yylim <- expanddotslist(ldots$ylim, np)
  zzlim <- expanddotslist(ldots$zlim, np)

  times <- x[isub,1]

  if (legend) {
    parplt <- par("plt") - c(0,0.07,0,0)
    parleg <- c(parplt[2]+0.02, parplt[2]+0.05, parplt[3], parplt[4])
    plt.or <- par(plt = parplt)
#    on.exit(par(plt = plt.or))
  }
  # Check if grid is increasing...
  if (! is.null(grid))
    gridOK <- min(diff (grid)) >0
  else
    gridOK <- TRUE

  if (! gridOK) grid <- rev(grid)

  # for each output variable (plot)
  for (ip in 1:np) {
    # ix     <- Which[ip]
    istart <- Select$istart[ip]
    istop  <- Select$istop[ip]
    if (gridOK)
      out    <- x[isub ,istart:istop]
    else
      out    <- x[isub ,istop:istart]

    dotmain      <- extractdots(Dotmain, ip)
    if (! is.null(xxlim)) dotmain$xlim <- xxlim[[ip]]
    if (! is.null(yylim)) dotmain$ylim <- yylim[[ip]]
    if (! is.null(zzlim))
      dotmain$zlim <- zzlim[[ip]]
    else
      dotmain$zlim <- range(out, na.rm=TRUE)
    List <- alist(z = out, x = times)
    if (! is.null(grid)) List$y = grid

    if (method == "persp") {
      if (is.null(dotmain$zlim))  # this to prevent error when range = 0
        if (diff(range(out, na.rm=TRUE)) == 0)
          dotmain$zlim <- c(0, 1)
      if (is.null(dotscol))
        dotmain$col <- drapecol(out, col = BlueRed (100), Range = dotmain$zlim)
      else
        dotmain$col <- drapecol(out, col = dotscol, Range = dotmain$zlim)

    } else if (method == "filled.contour")
      dotmain$color.palette <- dotscolorpalette
    else
      dotmain$col <- dotscol

    do.call(method, c(List, dotmain))
    if (Addcontour[ip]) do.call("contour", c(List, add = TRUE))

    if (legend) {
      if (method == "persp")
         if (is.null(dotscol))
           dotmain$col <- BlueRed(100)
         else
           dotmain$col <- dotscol
      if (is.null(dotmain$zlim)) dotmain$zlim <- range(out, na.rm=TRUE)
      drawlegend(parleg, dotmain)
    }
  }
  if (legend) {
      par(plt = plt.or)
      par(mar = par("mar")) # TRICK TO PREVENT R FROM SETTING DEFAULTPLOT = FALSE
  }
}

### ============================================================================
### plotting 2-D variables
### ============================================================================

plot.ode2D <- function (x, which, ask, add.contour, grid, method = "image",
   legend = TRUE, isub = 1:nrow(x), ...) {

  # Default color scheme
  BlueRed <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

  # if x is vector, check if there are enough columns ...
  att <- attributes(x)
  nspec <- att$nspec
  dimens <- att$dimens
  proddim <- prod(dimens)

  Mask <- att$map
  map  <- (! is.null(Mask))

  if (!map & (ncol(x) - nspec*proddim) < 1)
    stop("ncol of 'x' should be > 'nspec' * dimens if x is a vector")

  # variables to be plotted
  if (is.null(which))
    Which <- 1:nspec
  else
    Which <- which
  np <- length(Which)

  varnames <-  if (! is.null(att$ynames)) att$ynames else 1:nspec
  if (! is.null(att$lengthvar))
     varnames <- c(varnames, names(att$lengthvar)[-1])

  Select <- select2dvar(Which,varnames,att)
  Which <- Select$Which

  ldots <- list(...)
  Mtext <- ldots$mtext
  ldots$mtext <- NULL


  # number of figures in a row and interactively wait if remaining figures
  Ask <- setplotpar(ldots, np, ask)

  # here ask is always true by default...
  if (is.null(ask)) ask <- TRUE
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  N <- np * nrow(x)

  if (method == "filled.contour") {
    add.contour <- FALSE
    legend <- FALSE
  }
  Dotmain  <- setdots(ldots, N)  # expand dots to np values (no defaults)

  # different from the default
  Dotmain$main <- expanddots(ldots$main, varnames[Which], N)
  Dotmain$xlab <- expanddots(ldots$xlab,  "x"  , N)
  Dotmain$ylab <- expanddots(ldots$ylab,  "y"  , N)

  if (method == "persp")
    dotscol <- ldots$col
  else if (method == "filled.contour") {
    dotscolorpalette <- if (is.null(ldots$color.palette))
      BlueRed else ldots$color.palette
    dotscol <- dotscolorpalette(100)
    add.contour <- FALSE
    legend <- FALSE
  }  else if (is.null(ldots$col))
    dotscol <- BlueRed(100) else  dotscol <-  ldots$col

  dotslim <- ldots$zlim

  xxlim <- expanddotslist(ldots$xlim, np)
  yylim <- expanddotslist(ldots$ylim, np)
  zzlim <- expanddotslist(ldots$zlim, np)

  Addcontour <- rep(add.contour, length = np)

  i <- 0
  if (legend) {
    parplt <- par("plt") - c(0, 0.05, 0, 0)
    parleg <- c(parplt[2] + 0.02, parplt[2] + 0.05, parplt[3], parplt[4])
    plt.or <- par(plt = parplt)
#    on.exit(par(plt = plt.or))
  }
  x <- x[isub,]
  if (length(isub) > 1 & sum (isub) == 1)
      x <- matrix (nrow = 1, data =x)

  if (! is.null(Mtext))
    Mtext <- rep(Mtext, length.out = nrow(x))

  for (nt in 1:nrow(x)) {
    for (ip in 1:np) {
      i       <- i+1
      istart <- Select$istart[ip]
      istop  <- Select$istop[ip]
      if (map) {
        out <- rep (NA, length = prod(Select$dim[[ip]]))
        ii <- which (! is.na(Mask))
        out[ii] <- x[nt, istart:istop]
      } else
        out <- x[nt, istart:istop]

      dim(out) <- Select$dim[[ip]]

      dotmain    <- extractdots(Dotmain, i)
      if (! is.null(xxlim)) dotmain$xlim <- xxlim[[ip]]
      if (! is.null(yylim)) dotmain$ylim <- yylim[[ip]]
      if (! is.null(zzlim))
        dotmain$zlim <- zzlim[[ip]]
      else  {
        dotmain$zlim <- range(out, na.rm=TRUE)
        if (diff(dotmain$zlim ) == 0 )
          dotmain$zlim[2] <- dotmain$zlim[2] +1
      }
       if (map) {
          if (is.null(dotmain$zlim))
            dotmain$zlim <- range(out, na.rm=TRUE)
          out[is.na(out)] <- dotmain$zlim[1] - 0.01*max(1e-18,diff(dotmain$zlim))
          dotmain$zlim [1] <- dotmain$zlim[1] - 0.01*max(1e-18,diff(dotmain$zlim))
        }

      List <- alist(z = out)
      if (! is.null(grid)) {
        List$x <- grid$x
        List$y <- grid$y
      }

      if (method == "persp") {
        if (is.null(dotmain$zlim))
          if (diff(range(out, na.rm = TRUE)) == 0)
            dotmain$zlim <- c(0, 1)

        if (is.null(dotscol))
          dotmain$col <- drapecol(out, col = BlueRed(100), Range = dotmain$zlim)
        else
          dotmain$col <- drapecol(out, col = dotscol, Range = dotmain$zlim)

      } else if (method == "image") {
        dotmain$col <- dotscol
        if (map)           dotmain$col <- c("black", dotmain$col)
      } else if (method == "filled.contour")
        dotmain$color.palette <- dotscolorpalette

      do.call(method, c(List, dotmain))
      if (! method %in% c("persp", "filled.contour")) box()
      if (add.contour) do.call("contour", c(List, add = TRUE))

      if (legend) {
        if (method == "persp")
          if (is.null(dotscol))
            dotmain$col <- BlueRed(100)
          else
            dotmain$col <- dotscol
        if (is.null(dotmain$zlim))
          dotmain$zlim <- range(out, na.rm=TRUE)

        drawlegend(parleg, dotmain)
      }
    }
    if (! is.null(Mtext))
      mtext(outer = TRUE, side = 3, Mtext[nt],
             cex = 1.5, line = par("oma")[3]-1.5)

  }
  if (legend) {
        par(plt = plt.or)
        par(mar = par("mar")) # TRICK TO PREVENT R FROM SETTING DEFAULTPLOT = FALSE
  }
  # karline: ???   removed that... make it an argument?

  #  if (sum(par("mfrow") - c(1, 1)) == 0 )
  #   mtext(outer = TRUE, side = 3, paste("time ", x[nt, 1]),
  #         cex = 1.5, line = -1.5)
}

### ============================================================================
### Summaries of ode variables
### ============================================================================

summary.deSolve <- function(object, select = NULL, which = select,
   subset = NULL, ...){

  att  <- attributes(object)
  svar <- att$lengthvar[1]   # number of state variables
  lvar <- att$lengthvar[-1]  # length of other variables
  nspec <- att$nspec          # for models solved with ode.1D, ode.2D
  dimens <- att$dimens

  if (is.null(svar)) svar <- att$dim[2]-1  # models solved as DLL

  # variable names: information for state and ordinary variables is different
  if (is.null(att$ynames))
    if (is.null(dimens))
      varnames <- colnames(object)[2:(svar+1)]
    else
      varnames <- 1:nspec
  else
    varnames <- att$ynames   # this gives one name for multi-dimensional var.

  if (length(lvar) > 0) {
    lvarnames <- names(lvar)
    if (is.null(lvarnames))
      lvarnames <- (length(varnames)+1):(length(varnames)+length(lvar))
    varnames <- c(varnames, lvarnames)
  }

  # length of state AND other variables
  if (is.null(dimens))                 # all 0-D state variables
    lvar <- c(rep(1, len = svar), lvar)
  else
    lvar <- c(rep(prod(dimens), nspec), lvar) # multi-D state variables

  if (!missing(subset)){

   e <- substitute(subset)
   r <- eval(e, as.data.frame(object), parent.frame())
   if (is.numeric(r)) {
       isub <- r
   } else {
      if (!is.logical(r))
          stop("'subset' must evaluate to logical or be a vector with integers")
      isub <- r & !is.na(r)
      object <- object[isub,]
    }
  }

  # summaries for all variables
  Summ <- NULL
  for (i in 1:length(lvar)) {
    if (lvar[i] > 1) {
      Select <- select1dvar(i, varnames, att)
      out <- as.vector(object[, Select$istart:Select$istop])
    } else {
      Select <- selectvar(varnames[i], colnames(object), NAallowed = TRUE)
      if (is.na(Select))   # trick for composite names, e.g. "A.x" rather than "A"
         Select <- cumsum(lvar)[i]
      out <- object[ ,Select]
    }
  Summ <- rbind(Summ, c(summary(out, ...), N = length(out), sd = sd(out)))
  }
  rownames(Summ) <- varnames  # rownames or an extra column?
  if (! is.null(which))
    Summ <- Summ[which,]
  data.frame(t(Summ))         # like this or not transposed?
}

### ============================================================================
### Subsets of ode variables
### ============================================================================

subset.deSolve  <- function(x, subset = NULL, select = NULL,
  which = select, arr = FALSE, ...) {

  Which <- which # for compatibility between plot.deSolve and subset

  if (arr & length(Which) > 1)
    stop("cannot combine 'arr = TRUE' when more than one variable is selected")

  if (missing(subset))
    r <- TRUE
  else {
    e <- substitute(subset)
    r <- eval(e, as.data.frame(x), parent.frame())
    if (is.numeric(r)) {
      isub <- r
    } else {
      if (!is.logical(r))
        stop("'subset' must evaluate to logical or be a vector with integers")
      r <- r & !is.na(r)
    }
  }

  if (is.numeric(Which))
    return(x[r ,Which+1])

  if (is.null(Which))
    return(x[r , -1])         # Default: all variables, except time

  att   <- attributes(x)
  svar  <- att$lengthvar[1]   # number of state variables
  lvar  <- att$lengthvar[-1]  # length of other variables
  nspec <- att$nspec          # for models solved with ode.1D, ode.2D

  dimens <- att$dimens
  if (arr & length(dimens) <= 1    )
    warning("does not make sense to have 'arr = TRUE' when output is not 2D or 3D")

  if (is.null(svar)) svar <- att$dim[2]-1  # models solved as DLL
  if(is.null(nspec)) nspec <- svar

  # variable names: information for state and ordinary variables is different
  if (is.null(att$ynames))
    if (is.null(dimens))
      varnames <- colnames(x)[2:(svar+1)]
    else
      varnames <- 1:nspec
  else
    varnames <- att$ynames   # this gives one name for multi-dimensional var.
  varnames <- c("time",varnames)
  if (length(lvar) > 0) {
    lvarnames <- names(lvar)
    if (is.null(lvarnames))
      lvarnames <- (length(varnames)+1):(length(varnames)+length(lvar))
    varnames <- c(varnames, lvarnames)
  }

  # length of state AND other variables
  if (is.null(dimens))                 # all 0-D state variables
    lvar <- c(rep(1, len = svar), lvar)
  else
    lvar <- c(rep(prod(dimens), nspec), lvar) # multi-D state variables

  cvar <- cumsum(c(1,lvar))

  # Add selected variables to Out
  Out <- NULL
  for (iw in 1:length(Which)) {
    i <- which (varnames == Which[iw])
    if (length(i) == 0) {
      i <- which (colnames(x) == Which[iw])
      if (length(i) == 0)
         stop ("cannot find variable ", Which[iw], " in output")
      Out <- cbind(Out, x[,i])
    } else {
      if (is.null(i))
         stop ("cannot find variable ", Which[iw], " in output")
      istart <- 1
      if (i > 1) istart <- cvar[i-1]+1
      istop <- cvar[i]
       Out <- cbind(Out, x[ ,istart:istop])
    }
  }
  if (length(Which) == ncol(Out)) colnames(Out) <- Which
  OO    <- Out[r, ]
  if(is.vector(OO)) OO <- matrix(ncol = ncol(Out), data = OO)
  times <- x[r,1]

  if (arr & length(dimens) > 1 & ncol(OO) == prod(dimens)) {
     Nr <- nrow(OO)
     OO <- array(dim = c(dimens, Nr) , data = t(OO))
  }
   attr(OO, "times") <- times
   return(OO)
}
