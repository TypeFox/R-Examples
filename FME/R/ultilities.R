## =============================================================================
## Set the mfrow parameters and whether to "ask" for opening a new device
## =============================================================================

setplotpar <- function(nmdots, dots, nv, ask) {
  if (!any(match(nmdots, c("mfrow", "mfcol"), nomatch = 0))) {
    nc <- min(ceiling(sqrt(nv)), 3)
    nr <- min(ceiling(nv/nc), 3)
    mfrow <- c(nr, nc)
  } else if ("mfcol" %in% nmdots)
    mfrow <- rev(dots$mfcol)
  else
    mfrow <- dots$mfrow

  if (! is.null(mfrow)) mf <- par(mfrow = mfrow)

  ## interactively wait if there are remaining figures
  if (is.null(ask))
    ask <- prod(par("mfrow")) < nv && dev.interactive()

  return(ask)
}

## =============================================================================
## Panels to be used in pairs plots...
## =============================================================================

panel.cor <- function(x, y,...)
  text(x = mean(range(x)), y = mean(range(y)),
       labels = format(cor(x, y), digits = 2))

panel.hist <- function(x,...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 2))
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks
  nB <- length(breaks)
  y <- h$counts
  y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey")
}

## =============================================================================
## Find a certain variable
## =============================================================================

findvar <- function(var1, var2, str = "var") {
  if (is.character(var2[[1]])){
    ivar <- sapply(var2, FUN= function(x) which (names(var1) %in% x))
#    ivar  <- which (names(var1) %in% var2)  #returns sorted list
    if (length(ivar)!= length(var2))
      stop(paste("cannot proceed: not all sensitivity", str, "are known"))
    return(ivar)
  } else {
  if (max(var2) > length(var1))
    stop (paste("cannot proceed: index to sensitivity ", str, "too large"))
  return(var2)
  }
}

## =============================================================================
## Selecting numbers from "which"...
## =============================================================================

selectvar <- function (which, var, nm = "x", Nall = FALSE, NAallowed = FALSE) { # var = list from which to select...
  if (!is.null(which)) {
    if (! is.numeric(which)) {
      ln <- length(which)
      Select <- NULL
      for (i in which) {  # use loop rather than which(...%in%) to keep ordering of "which"
        ii <- which (var == i)
        if (length(ii) ==0 & ! NAallowed) 
          stop("variable ", i, " not in variable names")
        else if (length(ii) == 0)
          Select <- c(Select, NA)
        Select <- c(Select, ii)
      }
    } else {     # index
      Select <- which
      if (max(Select) > length(var))
        stop("index in 'which' too large")
      if (min(Select) < 1)
        stop("index in 'which' should be > 0")
    }
  } else if (is.null(which))
    if (Nall)
      Select <- 1:length(var)
    else
      Select <- NULL

  return(Select)
}
### ============================================================================
### first some common functions
### ============================================================================
# Update range, taking into account neg values for log transformed values
Range <- function(Range, x, log) {
   if (log) 
      x[x <= 0] <- min(x[x>0])  # remove zeros
   return( range(Range, x, na.rm = TRUE) )
}

## =============================================================================
## function for checking and expanding arguments in dots (...) with default
## =============================================================================

expanddots <- function (dots, default, n) {
  dots <- if (is.null(dots)) default else dots
  rep(dots, length.out = n)
}

# ks->Th: for xlim and ylim....
expanddotslist <- function (dots, n) {
  if (is.null(dots)) return(dots)
  dd <- if (!is.list(dots )) list(dots) else dots
  rep(dd, length.out = n)
}

## =============================================================================
## functions for expanding arguments in dots  (...)
## =============================================================================

repdots <- function(dots, n) 
  if (is.function(dots)) dots else rep(dots, length.out = n)

setdots <- function(dots, n) lapply(dots, repdots, n)

## =============================================================================
## function for extracting element 'index' from dots  (...)
## =============================================================================

extractdots <- function(dots, index) {
  ret <- lapply(dots, "[", index)
  ret <- lapply(ret, unlist) ## thpe: flatten list (experimental)
  return(ret)
}

