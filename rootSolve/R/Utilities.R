### ============================================================================
### ============================================================================
### S3 methods
### also possible to plot multiple outputs and to add observations.
### ============================================================================
### ============================================================================

  is.compiled <- function (FF)
   (is.character(FF) || class(FF) == "CFunc")

  CheckFunc <- function (FF)
   (is.function(FF) || is.character(FF) || class(FF) == "CFunc")


### ============================================================================
### first some common functions
### ============================================================================
# Update range, taking into account neg values for log transformed values
Range <- function(Range, x, log) {
   if (log) 
      x[x<=0] <- min(x[x>0])  # remove zeros
   return( range(Range, x, na.rm = TRUE) )
}

## =============================================================================
## function for checking and expanding arguments in dots (...) with default
## =============================================================================

expanddots <- function (dots, default, n) {
  dots <- if (is.null(dots)) default else dots
  rep(dots, length.out = n)
}

# for xlim and ylim....
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

## =============================================================================
## Set the mfrow parameters and whether to "ask" for opening a new device
## =============================================================================

setplotpar <- function(nmdots, dots, nv, ask) {
    if (!any(match(nmdots, c("mfrow", "mfcol"), nomatch = 0))) {
      nc <- min(ceiling(sqrt(nv)),3)
      nr <- min(ceiling(nv/nc),3)
      mfrow <- c(nr, nc)
    }
    else if ("mfcol" %in% nmdots)
        mfrow <- rev(dots$mfcol)
    else mfrow <- dots$mfrow

    if (! is.null(mfrow)) {
      mf <- par(mfrow=mfrow)
    }

   ## interactively wait if there are remaining figures
    if (is.null(ask))
      ask <- prod(par("mfrow")) < nv && dev.interactive()

    return(ask)
}

## =============================================================================
## find a variable
## =============================================================================

selectstvar <- function (which, var, NAallowed = FALSE) {

    if (!is.numeric(which)) {
        ln <- length(which)
        # keep ordering...
        Select <- NULL
        for ( i in 1:ln) {
          ss <- which(which[i] == var)
          if (length(ss) == 0 & ! NAallowed)
            stop("variable ", which[i], " not in variable names")
          else if (length(ss) == 0)
            Select <- c(Select,NA)
          else
            Select <- c(Select,ss)
        }        
    }
    else {
        Select <- which   # "Select now refers to the column number
        if (max(Select) > length(var))
            stop("index in 'which' too large")
        if (min(Select) < 1)
            stop("index in 'which' should be > 0")
    }
  return(Select)
}
### ============================================================================
### Merge two observed data files; assumed that first column = 'x' and ignored
### ============================================================================

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


mergeObs <- function(obs, Newobs) {
      
  if (! class(Newobs) %in% c("data.frame","matrix"))
    stop ("the elements in 'obs' should be either a 'data.frame' or a 'matrix'")
      
  if (is.character(Newobs[,1]) | is.factor(Newobs[,1])) 
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
  O1[ ,1] <- Newobs[,1]
  for (j in ii) {   # obseerved data in common are put in correct position
    jj <- which (obsname == newname[j])
    O1[,jj] <- Newobs[,j+1]
  }
  O1 <- cbind(O1, Newobs[,-c(1,ii+1)] )
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

### ============================================================================
### S3 methods
### ============================================================================

plot.steady1D <- function (x, ..., which = NULL, grid = NULL, 
  xyswap = FALSE, ask = NULL, obs = NULL, obspar = list(),
  vertical = FALSE ) {

## if x is vector, check it; put in correct format  
    checkX <- function (x) {
      X <- x$y
      if (is.vector(X)) {
        nspec <- attributes(x)$nspec
        if (length(X)%%nspec != 0) 
          stop("length of 'x' should be a multiple of 'nspec' if x is a vector")
        x <- matrix(ncol = nspec, data = X)
      } else x <- X   # only state variables
      if (is.null(colnames(x))) colnames(x) <- 1:ncol(x)
      x
    }
    if (! is.null(grid))
      if (! is.vector(grid))
        stop("'grid' should be a vector")
    
    preparex <- function(Which, xother, x, xx) {
      ii <- which (xother %in% Which) 
      if (length (ii) == length(Which))
        xx <- NULL
      if (length(ii) > 0)  {
        xnew <-  matrix(ncol = length(ii), data = unlist(x[ii+ 1]))
        colnames(xnew) <- xother[ii]
        xx <- cbind (xx, xnew)
      }  
     return(xx)
    }

    xx <- checkX (x)

## check observed data
    nobs <- 0

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
           obs.pos <- rbind(obs.pos, c(obs.pos[nrow(obs.pos),2] +1, nrow(obs)))
         }
       obsname <- colnames(obs) 
      } else {
       if (is.character(obs[,1]) | is.factor(obs[,1])) 
          obs <- convert2wide(obs)
       obsname <- colnames(obs) 
       if (! class(obs) %in% c("data.frame", "matrix"))
         stop ("'obs' should be either a 'data.frame' or a 'matrix'")
       obs.pos <- matrix(nrow = 1, c(1, nrow(obs)))
      }                       
    DD <- duplicated(obsname)
    if (sum(DD) > 0)  
      obs <- mergeObs(obs[,!DD], cbind(obs[,1],obs[,DD]))

    nobs <- nrow(obs.pos)   
    }

## variables to be plotted
    varnames <- c(colnames(xx),names(x)[-1])
    if(is.null(varnames)) varnames <- 1:ncol(xx)

    xother <- names(x)[-1]              # other names

    Which <- which 
    
    if (is.null(Which) & is.null(obs))  # All variables plotted
      Which <- 1:ncol(xx)
    else if (is.null(Which)) {          # All common variables in xx and obs plotted
     Which <- which(varnames %in% obsname)
     Which <- varnames[Which]           # names rather than numbers
     if (length (Which) == 0)
       stop ("observed data and model output have no variables in common")
    } 

## Some variables may not be state variables (not in x$y, but in other list values)
    xx <- preparex (Which, xother, x, xx)
    varnames <- colnames(xx)

## Position of variables to be plotted in "x" 
     if (length (Which) == 0)
       stop ("nothing to plot")

    xWhich <- selectstvar(Which,varnames)
    np <- length(xWhich)

    ldots <- list(...)
    ndots <- names(ldots)

## number of figures in a row and interactively wait if remaining figures
    ask <- setplotpar(ndots, ldots, np, ask)
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }

## Position of variables to be plotted in observed data
    if (! is.null(obs)) {
      ObsWhich <- selectstvar(varnames[xWhich], obsname, NAallowed = TRUE)
      ObsWhich [ ObsWhich > ncol(obs)] <- NA
    } else 
      ObsWhich <- rep(NA, length(xWhich))

## create two lists: x2:   other deSolve objects, 
##                   dots: remaining (plotting) parameters
    x2 <- list()
    dots <- list()
    nd <- 0
    nother <- 0                          # number of other steady1D instances
    
    if (length(ldots) > 0) 
     for ( i in 1:length(ldots))    
      # an element of type steady1D
      if ("steady1D" %in% class(ldots[[i]]))  {    
        x2[[nother <- nother+1]] <- ldots[[i]]  
        names(x2)[nother] <- ndots[i]       
      # a list of types steady1D
      } else if (is.list(ldots[[i]]) & "steady1D" %in% class(ldots[[i]][[1]])) {   
        for (j in 1:length(ldots[[i]])) {
          x2[[nother <- nother+1]] <- ldots[[i]][[j]]  
          names(x2)[nother] <- names(ldots[[i]])[[j]]
        }
      # a graphical parameter  
      } else if (! is.null(ldots[[i]])) {
        dots[[nd <- nd+1]] <- ldots[[i]]
        names(dots)[nd] <- ndots[i]
      }
      
## check compatibility of all steady1D objects    
    if (nother > 0) {
      for ( i in 1:nother) {             
        X <- checkX(x2[[i]])
        x2[[i]] <- preparex (Which, xother, x2[[i]], X)
        if (min(dim(x2[[i]]) - dim(xx) == c(0, 0)) == 0)   
          stop(" 'x2' and 'x' are not compatible - dimensions not the same")
        if (min(colnames(x2[[i]]) == varnames) == 0)
          stop(" 'x2' and 'x' are not compatible - colnames not the same")
      }
    } 
    
    nx <- nother + 1
    
    if (is.null(grid)) 
       grid <- 1:nrow(xx)
    if (length(grid) != nrow(xx)) 
      stop("length of grid (x-axis) should be = number of rows in 'x$y'")

## plotting parameters : split in plot parameters and point parameters
    plotnames <- c("xlab","ylab","xlim","ylim","main","sub","log","asp",
                   "ann","axes","frame.plot","panel.first","panel.last",
                   "cex.lab","cex.axis","cex.main")    
    
    # plot.default parameters
    ii <- names(dots) %in% plotnames
    dotmain <- dots[ii]
    dotmain <- setdots(dotmain, np)  # expand to np for each plot
    
    # these are different from the default
    dotmain$xlab <- expanddots(dots$xlab,  "x"      , np)
    dotmain$ylab <- expanddots(dots$ylab,  ""               , np)
    dotmain$main <- expanddots(dots$main,  varnames[xWhich] , np)

    yylim   <- expanddotslist(dots$ylim, np)
    xxlim   <- expanddotslist(dots$xlim, np)
    
    # point parameters
    ip <- !names(dots) %in% plotnames
    dotpoints <- dots[ip]
    dotpoints <- setdots(dotpoints, nx)   # expand all dots to nx values

    # these are different from default
    dotpoints$type <- expanddots(dots$type, "l",      nx)
    dotpoints$lty  <- expanddots(dots$lty, 1:nx,      nx)
    dotpoints$pch  <- expanddots(dots$pch, 1:nx,      nx)
    dotpoints$col  <- expanddots(dots$col, 1:nx,      nx)
    dotpoints$bg   <- expanddots(dots$bg,  1:nx,      nx)

    xyswap   <- rep(xyswap, length  = np)
    vertical <- rep(vertical, length = np)

    if (nobs > 0) 
      Obspar <- setdots(obspar, nobs)
    
## for each variable
    for (ip in 1:np) {
      i  <- xWhich[ip]
      io <- ObsWhich[ip] 

      # plotting parameters for deSolve output 1 (opens a plot)
      Dotmain   <- extractdots(dotmain, ip)
      Dotpoints <- extractdots(dotpoints, 1)

      Xlog <- Ylog <- FALSE
      if (! is.null(Dotmain$log)) { 
        Ylog  <- length(grep("y", Dotmain$log))
        Xlog  <- length(grep("x", Dotmain$log))
      }       
      if (vertical[ip])  {  # overrules other settings; vertical profiles
        xyswap[ip] = TRUE
        Dotmain$axes = FALSE
#        Dotmain$frame.plot = TRUE
        Dotmain$xlab = ""
        Dotmain$xaxs = "i"
        Dotmain$yaxs = "i"
      }
      if (! xyswap[ip]) {            

        if ( is.null (yylim[[ip]])){
          yrange <- Range(NULL, xx[, i], Ylog)
          if (nother>0) 
           for (j in 1:nother) 
             yrange <- Range(yrange, x2[[j]][,i], Ylog)
           if (! is.na(io)) 
             yrange <- Range(yrange, obs[,io], Ylog)
           Dotmain$ylim <- yrange
        } else 
          Dotmain$ylim <- yylim[[ip]]           
          if (Dotmain$xlab =="x" & Dotmain$ylab == "") {
              xl <- Dotmain$ylab
              Dotmain$ylab <- Dotmain$xlab
              Dotmain$xlab <- xl
          }
  
        if (is.null(xxlim[[ip]])) {
          xrange <- Range(NULL, grid, Xlog)
          if (! is.na(io)) 
            xrange <- Range(xrange, obs[,1], Xlog)
          Dotmain$xlim <- xrange

        } else 
          Dotmain$xlim <- xxlim[[ip]]  

         do.call("plot", c(alist(grid, xx[, i]), Dotmain, Dotpoints))

        if (nother>0)        # if other rootSolve outputs
          for (j in 2:nx) 
            do.call("lines", c(alist(grid, x2[[j-1]][, i]), 
                extractdots(dotpoints, j)) )

        if (! is.na(io))     # one or more observed data inputs
           for (j in 1: nobs) 
              if (length (i.obs <- obs.pos[j, 1]:obs.pos[j, 2]) > 0) 
                do.call("points", c(alist(obs[i.obs, 1], obs[i.obs, io]), 
                         extractdots(Obspar, j) ))        

      
      } else {  # xyswap
      
        if ( is.null (xxlim[[ip]])) {
          xrange <- Range(NULL, xx[, i], Xlog)
          if (nother>0) 
            for (j in 1:nother) 
              xrange <- Range(xrange, x2[[j]][,i], Xlog)

          if (! is.na(io)) 
            xrange <- Range(xrange, obs[,io], Xlog)
          
          Dotmain$xlim <- xrange
        
        } else 
          Dotmain$xlim <- xxlim[[ip]]
          
        if ( is.null (yylim[[ip]])){
          yrange <- Range(NULL, range(grid), Ylog)
          if (! is.na(io)) 
             yrange <- Range(yrange, obs[,1], Ylog)
          Dotmain$ylim <- rev(yrange)
        } else {
          Dotmain$ylim <- yylim[[ip]]  
          if (vertical[ip])
            Dotmain$ylim <- Dotmain$ylim[c(2,1)]
        }
        do.call("plot", c(alist(xx[, i], grid), Dotmain, Dotpoints))
        if (vertical[ip]) {
          abline(h=Dotmain$ylim[2])
          abline(v=Dotmain$xlim[1])
          axis(side = 2)
          axis(side = 3, mgp = c(3,0.5,0))
        }
        if (nother>0)       # if other rootSolve outputs
          for (j in 2:nx) 
            do.call("lines", c(alist(x2[[j-1]][, i],grid), 
                extractdots(dotpoints, j)))
        
        if (! is.na(io))   # one or more observed data inputs
           for (j in 1: nobs) 
              if (length (i.obs <- obs.pos[j, 1]:obs.pos[j, 2]) > 0) 
                 do.call("points", c(alist(obs[i.obs, io], obs[i.obs, 1]), 
                          extractdots(Obspar, j)))        

     } # end xyswap
   }
}


### ============================================================================
## to draw a legend
### ============================================================================

drawlegend <- function (parleg, dots) {
        Plt <- par(plt = parleg)
        usr <- par("usr")
        par(new = TRUE)
        ix <- 1
        minz <- dots$zlim[1]
        maxz <- dots$zlim[2]
        binwidth <- (maxz - minz)/64
        iy <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
        iz <- matrix(iy, nrow = 1, ncol = length(iy))

        image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
              ylab = "", col = dots$col)
      
        do.call("axis", list(side = 4, mgp = c(3,1,0), las=2))
      
        par(plt = Plt)
        par(usr = usr)
        par(new = FALSE)
}

### ============================================================================
## to drape a color over a persp plot.
### ============================================================================

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

### ============================================================================

image.steady2D <- function (x, which = NULL, 
    add.contour = FALSE, grid = NULL, ask = NULL, method = "image", 
    legend = FALSE, ...) {

## Default color scheme
   BlueRed <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

# if x is vector, check if there is more than one species...  
    X <- x$y
    out <- list()
    nspec <- attributes(x)$nspec
    dimens <- attributes(x)$dimens  
    if (is.vector(X)) {
      if (length(X) - nspec*prod(dimens) != 0) 
        stop("length of 'x' should be = 'nspec' * prod(dimens) if x is a vector")
      # x <- matrix(ncol = nspec, data = X)
      
      for ( i in 1:nspec){
        istart <- (i-1)*prod(dimens) 
        out[[i]] <- matrix(nrow=dimens[1], ncol=dimens[2], data =
          X[(istart+1):(istart+prod(dimens))])
      }
    } else 
        out <- X   # only one state variable

    map <- any(is.na(X))    # if mapping applied: some elements will be NA.
     
    varnames <- attributes(x)$ynames
    if (is.null(varnames)) varnames <- 1:nspec

    if (length(x) > 1) {
       for ( ii in 2:length(x)) {
        out[[i+ii-1]] <- x[[ii]]
        varnames <- c(varnames,names(x)[ii])
       }
      }          

# ADD NON-STATE VARIABLES...      
    if (is.null(which)) which <- 1:nspec
    which <- selectstvar(which,varnames)
    
    np <- length(which)

    dots <- list(...)
    nmdots <- names(dots)

    # number of figures in a row and 
    # interactively wait if there are remaining figures
   
    ask <- setplotpar(nmdots, dots, np, ask)
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    
   Dots  <- setdots(dots, np)   # expand dots to np values (no defaults)

    # different from the default
   Dots$main  <- expanddots(dots$main, varnames[which], np)
   Dots$xlab  <- expanddots(dots$xlab, "x",  np)
   Dots$ylab  <- expanddots(dots$ylab, "y",  np)

   # colors - different if persp, image or filled.contour
    
   if (method == "persp")
      dotscol <- dots$col

   else if (method == "filled.contour")  {
      dotscolorpalette <- if (is.null(dots$color.palette))
        BlueRed else dots$color.palette
      dotscol <- dotscolorpalette(100)
      add.contour <- FALSE
      legend <- FALSE
   } else
     if (is.null(dots$col))
       dotscol <- BlueRed(100) else dotscol <- dots$col

   Addcontour <- rep(add.contour, length = np)

## xlim, ylim and zlim are special:  
   xxlim <- expanddotslist(dots$xlim, np)
   yylim <- expanddotslist(dots$ylim, np)
   zzlim <- expanddotslist(dots$zlim, np)

   if (legend) {
      parplt <- par("plt") - c(0,0.07,0,0) 
      parleg <- c(parplt[2]+0.02, parplt[2]+0.05, parplt[3], parplt[4])
      plt.or <- par(plt = parplt)
#      on.exit(par(plt = plt.or))
   }  

    for (i in 1:np) {
        ii <- which[i]

        dots      <- extractdots(Dots, i)
        if (! is.null(xxlim)) dots$xlim <- xxlim[[i]]
        if (! is.null(yylim)) dots$ylim <- yylim[[i]]
        if (! is.null(zzlim)) 
          dots$zlim <- zzlim[[i]]
        else
          dots$zlim <- range(out[[ii]], na.rm=TRUE)
          
        List <- alist(z=out[[ii]])
        if (! is.null(grid)) {
          List$x <- grid[[1]]
          List$y <- grid[[2]]
        }

        
        if (method=="persp") {
          if (is.null(dots$zlim))  # this to prevent error when range = 0
            if (diff(range(out[[i]], na.rm=TRUE)) == 0) 
              dots$zlim <- c(0, 1)

          if(is.null(dotscol))
             dots$col <- drapecol(out[[i]], col = BlueRed (100), Range = dots$zlim)
          else
             dots$col <- drapecol(out[[i]], col = dotscol, Range = dots$zlim)

        } else if (method == "filled.contour")
          dots$color.palette <- dotscolorpalette
        else 
          dots$col <- dotscol 

        if (! is.null(map)) {
          if (is.null(dots$zlim))
            dots$zlim <- range(out[[i]], na.rm=TRUE)
          dots$col <- c("black", dots$col)
          out[[i]][is.na(out[[i]])] <- dots$zlim[1] - 0.01*diff(dots$zlim)
          dots$zlim [1] <- dots$zlim[1] - 0.01*diff(dots$zlim)
        }
        
        do.call(method, c(List, dots)) 
        if (Addcontour[i]) do.call("contour", c(List, add=TRUE))
        box(lwd = 2)  # Karline: added that          
      if (legend) {
        if (method == "persp") 
           if (is.null(dotscol))
             dots$col <- BlueRed(100)
           else
              dots$col <- dotscol
        if (is.null(dots$zlim)) dots$zlim <- range(out, na.rm=TRUE)

        drawlegend(parleg, dots)      
        
      }    
   }
  if (legend)  {
        par(plt = plt.or)  
      par(mar = par("mar")) # TRICK TO PREVENT R FROM SETTING DEFAULTPLOT = FALSE
  }
}

### ============================================================================

image.steady3D <- function (x, which = NULL, dimselect = NULL,
    add.contour = FALSE, grid = NULL, ask = NULL, 
    method="image", legend = FALSE, ...) {

## Default color scheme
   BlueRed <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
# if x is vector, check if there is more than one species...  
    X <- x$y
    out    <- list()
    nspec  <- attributes(x)$nspec
    dimens <- attributes(x)$dimens
    Nx <- dimens[1]
    Ny <- dimens[2]
    Nz <- dimens[3]
    if (is.vector(X)) {
      if (length(X) - nspec*prod(dimens) != 0) 
        stop("length of 'x' should be = 'nspec' * prod(dimens) if x is a vector")
      x <- matrix(ncol = nspec, data = X)
      
      for ( i in 1:nspec){
        istart <- (i-1)*prod(dimens) 
        out[[i]] <- array(dim = dimens, data =
          X[(istart+1):(istart+prod(dimens))])
      }
    } else 
        out <- X   # only state variables
      
    if (is.null(which)) which <- 1:nspec
    varnames <- 1:nspec
    which <- selectstvar(which,varnames)
    
    np <- length(which)

    dots <- list(...)
    nmdots <- names(dots)

    # number of figures in a row and 
    # interactively wait if there are remaining figures
   
    ask <- setplotpar(nmdots, dots, np, ask)
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    
   Dots  <- setdots(dots, np)   # expand dots to np values (no defaults)
   ismain <- !is.null(dots$main)
    # different from the default
   Dotsmain  <- expanddots(dots$main, varnames[which], np)
   Dots$xlab  <- expanddots(dots$xlab, "x",  np)
   Dots$ylab  <- expanddots(dots$ylab, "y",       np)

   # colors - different if persp, image or filled.contour
    
   if (method == "persp")
      dotscol <- dots$col

   else if (method == "filled.contour")  {
      dotscolorpalette <- if (is.null(dots$color.palette))
        BlueRed else dots$color.palette
      dotscol <- dotscolorpalette(100)
      add.contour <- FALSE
      legend <- FALSE
   } else
     if (is.null(dots$col))
       dotscol <- BlueRed(100) else dotscol <- dots$col

   Addcontour <- rep(add.contour, length = np)

## xlim, ylim and zlim are special:  
   xxlim <- expanddotslist(dots$xlim, np)
   yylim <- expanddotslist(dots$ylim, np)
   zzlim <- expanddotslist(dots$zlim, np)

   if (legend) {
      parplt <- par("plt") - c(0,0.07,0,0) 
      parleg <- c(parplt[2]+0.02, parplt[2]+0.05, parplt[3], parplt[4])
      plt.or <- par(plt = parplt)
#      on.exit(par(plt = plt.or))
   }  
  
    dselect <- 1:Nz
    sel <- 3
    Nn <- Nz
    if(!is.null(dimselect$z)) {
      dselect <- dimselect$z
    } else if (! is.null(dimselect$y)) {
      dselect <- dimselect$y
      sel <- 2
      Nn <- Ny
    } else if (! is.null(dimselect$x)) {
      dselect <- dimselect$x
      sel <- 1
      Nn <- Nx
    }
    
    if (max(dselect) > Nn) 
        stop("Numbers in 'dimselect' can not be larger than corresponding dimension")
    if (min(dselect) < 1) 
        stop("Numbers in 'zselect' can not be < 1")
    
    for (d in dselect) {
      for (i in 1:np) {
        ii <- which[i]
        dots      <- extractdots(Dots, i)
        if (ismain)
           dots$main <-  Dotsmain[i]
        else   
           dots$main <- paste("var", Dotsmain[i], "dim ",sel," = ", d)
        if (sel == 1)
          zdat <- out[[i]][d, , ]
        else if (sel == 2)
          zdat <- out[[i]][, d, ]
        else if (sel == 3)
          zdat <- out[[i]][, , d]

        if (! is.null(xxlim)) dots$xlim <- xxlim[[i]]
        if (! is.null(yylim)) dots$ylim <- yylim[[i]]

        if (! is.null(zzlim)) 
          dots$zlim <- zzlim[[i]]
        else
          dots$zlim <- range(out[[ii]], na.rm = TRUE)
        
        if (method=="persp") {
          if (is.null(dots$zlim))  # this to prevent error when range = 0
            if (diff(range(out[[i]], na.rm = TRUE)) == 0) 
              dots$zlim <- c(0, 1)

          if(is.null(dotscol))
             dots$col <- drapecol(zdat, col = BlueRed (100), Range = dots$zlim)
          else
             dots$col <- drapecol(zdat, col = dotscol, Range = dots$zlim)

        } else if (method == "filled.contour")
          dots$color.palette <- dotscolorpalette
        else 
          dots$col <- dotscol 
        
        List <- list()
        if (! is.null(grid)) {
          List$x <- grid[[1]]
          List$y <- grid[[2]]
        }
          List$z <- zdat
        
          do.call(method, c(List, dots)) 
          if (Addcontour[i]) do.call("contour", c(List, add = TRUE))

          if (legend) {
            if (method == "persp") 
            if (is.null(dotscol))
              dots$col <- BlueRed(100)
            else
               dots$col <- dotscol
            if (is.null(dots$zlim)) dots$zlim <- range(out, na.rm=TRUE)

            drawlegend(parleg, dots)      
            
          }   
        }
     }   
   if (legend) {
      par(plt = plt.or)  
      par(mar = par("mar")) # TRICK TO PREVENT R FROM SETTING DEFAULTPLOT = FALSE
   }
}







### ============================================================================

subset.steady2D <- function (x, which = NULL, ... ) {

# if x is vector, check if there is more than one species...  
    X <- x$y
    out <- list()
    nspec <- attributes(x)$nspec
    dimens <- attributes(x)$dimens  
    if (is.vector(X)) {
      if (length(X) - nspec*prod(dimens) != 0) 
        stop("length of 'x' should be = 'nspec' * prod(dimens) if x is a vector")
      # x <- matrix(ncol = nspec, data = X)
      
      for ( i in 1:nspec){
        istart <- (i-1)*prod(dimens) 
        out[[i]] <- matrix(nrow=dimens[1], ncol=dimens[2], data =
          X[(istart+1):(istart+prod(dimens))])
      }
    } else 
        out <- X   # only one state variable

    map <- any(is.na(X))    # if mapping applied: some elements will be NA.
     
    varnames <- attributes(x)$ynames
    if (is.null(varnames)) varnames <- 1:nspec

    if (length(x) > 1) {
       for ( ii in 2:length(x)) {
        out[[i+ii-1]] <- x[[ii]]
        varnames <- c(varnames,names(x)[ii])
       }
      }          

# ADD NON-STATE VARIABLES...      
    if (is.null(which)) which <- 1:nspec
    which <- selectstvar(which,varnames)

    if (length(which) > 1)
      stop("Can only select one variable")
    out[[which]] 
}
