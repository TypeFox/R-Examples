## =============================================================================
## =============================================================================
## Matrix plotting - most of these function were slightly modified from 
## similar functions in the deSolve and rootSolve packages.
## =============================================================================
## =============================================================================

## ============================================================================
## Finds common variables in a set of matrices
## ============================================================================

Mcommon <- function (M, ..., verbose = FALSE) {     

# creates a list with subsets of matrices that have only common variables

  Colnames <- function(M) 
    if (is.null (cn <- colnames(M))) return (1:ncol(M)) else return(cn)

  LL <- list(...)
  if (length(LL) == 0)
    return (list(M))

  if (is.list(M)) {   
    if (! class(M[[1]]) %in% c("matrix", "data.frame"))
      stop ("elements in list 'M' should be either a 'matrix' or 'data.frame' ")
    LL <- c(M, LL)
  } else {
    if (! class(M) %in% c("matrix", "data.frame"))
      stop ("'M' should be either a 'matrix' or 'data.frame' or a 'list'")
    LL <- c(list(M), LL)
   #dirty trick to get ALL names of ellipsis 
    NN <- deparse(substitute(x(...)))
    NN <- gsub("x(","",NN,fixed=TRUE)
    NN <- gsub(")","",NN)
    NN <- gsub(" ","",NN)    
    dotnames <- unlist(strsplit(NN, ","))
    names(LL) <- c(deparse(substitute(M)), dotnames)                
    
  }
  cn <- Colnames(LL[[1]])
  for (i in 2:length(LL))  {
    if (! class(LL[[i]]) %in% c("matrix", "data.frame"))
      stop ("elements in '...' should be either a 'matrix' or 'data.frame'")
    cn <- cn[cn %in% Colnames(LL[[i]])]
  }
  if (verbose) 
    print(paste("common variable names: ", paste(cn, collapse = ", ")))
  if (length(cn) ==  0) {
    if (verbose) 
      warning("No variable names in common - returning 'NULL'")
    LL <- NULL
  } else {
    if (length(cn) ==  1 & verbose) 
      warning("Only one variable name in common - returning list of vectors")
    for (i in 1:length(LL))  
      LL[[i]] <- LL[[i]][, cn]
  }
  return(LL)
}

## ============================================================================
## Splits a matrix according to values in 'split'; the result is a list
## ============================================================================

Msplit <- function(M, 
                   split = 1, 
                   subset = NULL) {
  LL <- list()
  
  # quick and dirty
  if (!missing(subset)){
    e <- substitute(subset)
    r <- eval(e, as.data.frame(M), parent.frame())
    if (!is.logical(r))
      stop("'subset' must evaluate to logical")
    isub <- r & !is.na(r)
    M <- M[isub, ]
  }  

  ux <- unique(M[ , split])
  if (length(split) == 1) 
    ux <- matrix(ncol = 1, data = ux)
  lux <- nrow(ux)

  isel <- 1 : ncol(M)
  if (is.numeric(split)) 
    isel <- isel[-split]
  else
    isel <- isel[ -which(colnames(M) %in% split)]

  for (i in 1:lux) {
    Sel <- M 
    for (j in 1:length(split))
      Sel <- Sel[Sel[ ,split[j]] == ux[i,j], ]
 
    LL[[i]] <- as.data.frame(Sel[,isel] )
  }  

  lnames <- ux[,1]
  if (length(split) > 1) 
    for (i in 2:length(split))  
      lnames <- paste(lnames, ux[,i])

  names(LL) <- lnames

  LL
}

## =============================================================================
## Plot a (list of) matrices
## =============================================================================

Mplot <- function (M, ..., 
                   x = 1, 
                   select = NULL, which = select, 
                   subset = NULL, ask = NULL, 
                   legend = list(x = "center"),
                   pos.legend = NULL,
                   xyswap = FALSE, rev = "") {

  getnames <- function(x) {
    if (is.null (cn <- colnames(x))) return (1:ncol(x)) else return(cn)
  }

  plotlegend <- function () {
    if (nolegend) return()
       
    # Add legend, if legend not equal to NULL or not equal to FALSE
    if (! is.list(legend)) {
      if (legend[1] == FALSE) 
        legend <- list()
        else if (legend[1] == TRUE)  
          legend <- list(x = "top")
    }
  
    if (length(legend) > 0) {
      if (!is.list(legend)) 
        stop ("'legend' should be a list or NULL")

      if (is.null(legend$col))
        legend$col <- Dotpoints$col

      if (is.null(legend$pt.bg))
        legend$pt.bg <- Dotpoints$pt.bg

      if (is.null(legend$lwd))
        legend$lwd <- Dotpoints$lwd

      if (is.null(legend$pch)) {
        legend$pch <- Dotpoints$pch
        legend$pch[Dotpoints$type == "l"] <- NA
      }
    
      if (is.null(legend$pt.cex))
        legend$pt.cex <- Dotpoints$cex

      if (is.null(legend$lty)) {
        legend$lty <- Dotpoints$lty
        legend$lty[Dotpoints$type == "p"] <- NA
        if (all(is.na(legend$lty)))
          legend$lty <- NULL
      }

      if (is.null(legend$legend))
        legend$legend <- names(x2)
      if (is.null(legend$x))
        legend$x <- "center"
    
      do.call("legend", legend)
    }
  }

                  
  # The ellipsis
  ldots   <- list(...)
  
  mtext <- ldots$mtext
  ldots$mtext <- NULL
  
  Dots    <- splitdots(ldots, M)
  x2      <- Dots$x2
  nother  <- Dots$nother
  nx      <- nother + 1 # total number of objects to be plotted
  varnames <- getnames(x2[[1]])
  
  nolegend <- (nother == 0 & is.null(pos.legend))

 # x-variable  
  xPos <- vector()
  for (i in 1: length(x2))
    xPos[i] <- selectvar(x, getnames(x2[[i]]))
  xisfactor <- is.factor(x2[[1]][,xPos[1]])
  xname <- varnames[xPos[1]]

 # variables to be plotted
  Which <- which
  if (is.null(Which)) {
    for (i in 1: length(x2))
      Which <- c(Which,getnames(x2[[i]])[- xPos[i]])
    Which <- unique(Which)
  }

  np      <- length(Which)
  if (np == 0)  
    stop ("M cannot be a (list of) vector(s)")
  if (!is.null(pos.legend)) {
    if (is.character(pos.legend)) 
      pos.legend <- which(pos.legend == Which)
    if (pos.legend > np)
      stop("'pos.legend' should be referring to a variable name or number to plot")
  } else
    pos.legend <- np
  
 # Position of variables to be plotted in "M" and other matrices
  xWhich <- list()

  for (i in 1: length(x2))
    xWhich[[i]] <- selectvar(Which, getnames(x2[[i]]))

  if (! is.character(Which)) 
    Which <- varnames[xWhich[[1]]]

  # number of figures in a row and interactively wait if remaining figures
  ask <- setplotpar(ldots, np + (pos.legend == 0), ask)                   
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  Dotmain <- setdots(Dots$main, np)  # expand to np for each plot

 # these are different from the default
  Dotmain$xlab <- expanddots(ldots$xlab, xname      , np)
  Dotmain$ylab <- expanddots(ldots$ylab, ""         , np)
  Dotmain$main <- expanddots(ldots$main, Which      , np)

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
    isub <- list()
    for (i in 1:nx) {
      e <- substitute(subset)
      r <- eval(e, as.data.frame(x2[[i]]), parent.frame())
      if (!is.logical(r))
        stop("'subset' must evaluate to logical")
      isub[[i]] <- r & !is.na(r)
    }  
  } else isub <- rep(TRUE, nx)

  xyswap <- rep(xyswap, length = np)
  rev <- rep(rev, length = np)

 # LOOP for each output variable (plot)
  for (ip in 1 : np) {

   # plotting parameters for matrix output 1 (opens a plot)
    dotmain   <- extractdots(Dotmain, ip)
    dotpoints <- extractdots(Dotpoints, 1)  # 1st dotpoints

    Xlog <- Ylog <- FALSE
    if (! is.null(dotmain$log)) {
      Ylog  <- length(grep("y",dotmain$log))
      Xlog  <- length(grep("x",dotmain$log))
    }

   # first object plotted (new plot created)
    ix <- xWhich[[1]][[ip]]      # position of variable in 'x'

    if (! xyswap[ip]) {
      if (is.null(yylim[[ip]]))
        dotmain$ylim <- SetRange(yylim[[ip]], x2, isub, xWhich, ip, Ylog)
      else
        dotmain$ylim <- yylim[[ip]]

      if (is.null(xxlim[[ip]])) {
        dotmain$xlim <- SetRange(xxlim[[ip]], x2, isub, xPos, 1, Xlog)
        if (xisfactor) 
          dotmain$xlim <- dotmain$xlim + c(-0.5, 0.5)
      } else
        dotmain$xlim <- xxlim[[ip]]
    } else {
      if (is.null(xxlim[[ip]]))
        dotmain$xlim <- SetRange(NULL, x2, isub, xWhich, ip, Ylog)
      else
        dotmain$xlim <- xxlim[[ip]]

      if (is.null(yylim[[ip]])) {
        dotmain$ylim <- SetRange(NULL, x2, isub, xPos, 1, Xlog)
        if (xisfactor) 
          dotmain$ylim <- dotmain$ylim + c(-0.5, 0.5)
      } else
        dotmain$ylim <- yylim[[ip]]
    }
    if (length(grep("x",rev[ip])))
      dotmain$xlim <- rev(dotmain$xlim)
    if (length(grep("y",rev[ip])))
      dotmain$ylim <- rev(dotmain$ylim)

    if (all(is.na(x2[[1]][isub[[1]], ix]))) {
      dotpoints$type <- dotmain$axes <- dotmain$xlab <- dotmain$ylab <-  NULL
      plot(x = 0.5, y = 0.5, type = "n", main = dotmain$main)
      text (0.5, 0.5, labels = "No data")  
    } else if (xyswap[ip]) {
      do.call("plot", c(alist(y = x2[[1]][isub[[1]], xPos[1]], 
        x = x2[[1]][isub[[1]], ix]), dotmain, dotpoints))
    } else
      do.call("plot", c(alist(x = x2[[1]][isub[[1]], xPos[1]], 
        y = x2[[1]][isub[[1]], ix]), dotmain, dotpoints))

    if (nother > 0)        # if other outputs
      for (j in 2:nx) {
        ix <- xWhich[[j]][[ip]]      # position of variable in 'x2'
        if (!is.na(ix)) {
          xx <- x2[[j]][isub[[j]], xPos[j]]
          yy <- x2[[j]][isub[[j]], ix]
          ii <- c(which(is.na(yy)), which(is.na(xx)))
          if (length(ii) > 0) {
            xx <- xx[-ii]
            yy <- yy[-ii]
          }
            
         if (xyswap[ip]) 
          do.call("lines", c(alist(y = xx, 
                x = yy), extractdots(Dotpoints, j)) )

         else
          do.call("lines", c(alist(x = xx, 
                y = yy), extractdots(Dotpoints, j)) )
        }
      }
     if (pos.legend == ip)
       plotlegend()
     
  }
  
  if (! is.null(mtext))
    mtext(outer = TRUE, side = 3, mtext, line = par()$oma[3]-1, 
          cex = par()$cex*1.5)
  if (pos.legend == 0) {
    plot.new()
    plotlegend()
  }   

}

## =============================================================================
## Update range, taking into account neg values for log transformed values
## =============================================================================

Range <- function(Range, x, log) {
  if (is.null(x)) 
    return(Range)
  if (log)
    x[x <= 0] <- min(x[x>0])  # remove zeros
   
   RR <- range(Range, as.double(x), na.rm = TRUE)
   RR[is.infinite(RR)]<- NA

   return( RR )
}


SetRange <- function(lim, x2, isub, xWhich, ip, Log) {

  nx <- length (x2)
  if ( is.null (lim)) {
    yrange <- NULL
      for (j in 1:nx){
        ix <- xWhich[[j]][ip]
        if (! all(is.na(x2[[j]][isub[[j]],ix])))
          yrange <- Range(yrange, x2[[j]][isub[[j]],ix], Log)
      }  
  } else
    yrange  <- lim

  return(yrange)
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
  else 
    mfrow <- ldots$mfrow

  if (! is.null(mfrow))  
    mf <- par(mfrow = mfrow)

  ## interactively wait if there are remaining figures
  if (is.null(ask))
    ask <- prod(par("mfrow")) < nv && dev.interactive()

  return(ask)
}

## =============================================================================
## find a variable  - and keep the ordering
## =============================================================================

selectvar <- function (Which, var, NAallowed = TRUE) {
  if (!is.numeric(Which)) {
    ln <- length(Which)

   # the loop is necessary so as to keep ordering...
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
    Select <- Which  # "Select" now refers to the column number
    if (max(Select) > length(var))
      stop("index in 'which' too large: ", max(Select)-1)
    if (min(Select) < 1)
      stop("index in 'which' should be > 0")
  }
  return(Select)
}

## ============================================================================
## create several lists: x2:   other matrix objects,
##                       dotmain, dotpoints: remaining (plotting) parameters
## ============================================================================

splitdots <- function(ldots, x){
  x2      <- list()
  nother <- 0
  islist <- (! is.data.frame(x) & is.list(x))
  
  if (! islist) {
    x2[[1]] <- x
    names(x2)[1] <-"M"
  } else {
    for(i in 1:length(x))
      x2[[i]] <- x[[i]]
    names(x2) <- names(x)
    nother <- length(x) - 1
  }

  dots   <- list()
  nd     <- 0
  ndots <- names(ldots)
    
  if (length(ldots) > 0)
    for ( i in 1:length(ldots))
      if ("matrix" %in% class(ldots[[i]]) | "data.frame" %in% class(ldots[[i]])) { 
        nother <- nother + 1        
        x2[[nother + 1]] <- ldots[[i]]
        if (is.null(ndots[i]))
          names(x2)[nother+1] <- nother 
        else 
          names(x2)[nother+1] <- ndots[i]
        # a list of matrix objects
      } else if (is.list(ldots[[i]]) & 
        ("matrix" %in% class(ldots[[i]][[1]]) | 
         "data.frame" %in% class(ldots[[i]][[1]]))) {
        for (j in 1:length(ldots[[i]])) {
          nother <- nother + 1        
          x2[[nother+1]] <- ldots[[i]][[j]]
          nn <- names(ldots[[i]])[[j]]
          if (is.null(nn)) 
            nn <- nother
          names(x2)[nother+1] <- nn
        }
      } else if (! is.null(ldots[[i]])) {  # a graphical parameter
        dots[[nd <- nd+1]] <- ldots[[i]]
        names(dots)[nd] <- ndots[i]
      }

  nmdots <- names(dots)

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
  list (points = dotpoints, main = dotmain, nother = nother, x2 = x2)
}
