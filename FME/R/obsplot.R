### ============================================================================
### ============================================================================
### Plotting observed data
### also possible to plot multiple outputs and to add observations.
### ============================================================================
### ============================================================================


### ============================================================================
### first some common functions
### ============================================================================
# Update range, taking into account neg values for log transformed values
Range <- function(Range, x, log) {
   if (log) 
      x[x<0] <- min(x[x>0])  # remove zeros
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
      nt  <- cbind(sel[,1],matrix(nrow = nrow(sel), ncol = ncol(MAT)-1, data = NA),sel[,2])
      MAT <- cbind(MAT, NA)
      colnames(nt) <- colnames(MAT)
      MAT <- rbind(MAT, nt)
      colnames.MAT <- c(colnames.MAT, ivar)
    }  
  colnames(MAT) <- colnames.MAT
  return(MAT)
}

### ============================================================================

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

obsplot <- function (x, ..., which = NULL, xyswap = FALSE, ask = NULL) {

## check observed data
    checkobs <- function (obs) {
      Obs <- obs
      obsname <- colnames(obs) 
        if (! class(obs) %in% c("data.frame", "matrix"))
          stop ("'obs' should be either a 'data.frame' or a 'matrix'")
      DD <- duplicated(obsname)
      if (sum(DD) > 0) {   # Karline: changed this to account for more columns 
        wD <- (1:ncol(obs))[DD]
        for (id in wD) {
          Add <- cbind(Obs[,1],Obs[,id]); colnames(Add) <- obsname[c(1,id)]
          obs <- mergeObs(obs[,!DD], Add)
        }  
      }        
      else  if (is.character(obs[,1]) | is.factor(obs[,1])) 
        obs <- convert2wide(obs)
      return(obs)
    } 

    obs <- checkobs(x)
       
    ldots <- list(...)
    ndots <- names(ldots)

## create two lists: x2:   other matrix.data.frame objects, 
##                   dots: remaining (plotting) parameters
    dots <- list()
    nd <- 0
    
    obs.pos <- matrix(nrow = 1, c(1, nrow(obs)))

    if (length(ldots) > 0) 
     for ( i in 1:length(ldots))
      if (class(ldots[[i]]) %in% c("matrix", "data.frame") ) {
        NewObs <- ldots[[i]]  
        obs    <- mergeObs(obs, NewObs)
        obs.pos   <- rbind(obs.pos, c(obs.pos[nrow(obs.pos),2] +1, nrow(obs)))
      } else if (! is.null(ldots[[i]])) {
        dots[[nd <- nd+1]] <- ldots[[i]]
        names(dots)[nd] <- ndots[i]
      }
    nx      <- nrow(obs.pos)         # number of observed data...
    obsname <- colnames(obs)
    
## variables to be plotted
    Which <- which 
    
    if (is.null(Which))        # All variables plotted
      Which <- 2:(length(obsname))
    else if (is.numeric(Which))
      Which <- Which + 1
## Position of variables to be plotted in observed data
    ObsWhich <- selectvar(Which, obsname)  
    np <- length(ObsWhich)

## number of figures in a row and interactively wait if remaining figures
    nmdots <- names(dots)    
    ask <- setplotpar(ndots, ldots, np, ask)
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }

    labs <- (is.null(dots$xlab) && is.null(dots$ylab))

## plotting parameters : split in plot parameters and point parameters
    plotnames <- c("xlab","ylab","xlim","ylim","main","sub","log","asp",
                   "ann","axes","frame.plot","panel.first","panel.last",
                   "cex.lab","cex.axis","cex.main")    
    
    # plot.default parameters
    ii <- names(dots) %in% plotnames
    dotmain <- dots[ii]
    dotmain <- setdots(dotmain, np)  # expand to np for each plot

    # these are different from the default
    dotmain$xlab <- expanddots(dots$xlab,  "x"                , np)
    dotmain$ylab <- expanddots(dots$ylab,  ""                 , np)
    dotmain$main <- expanddots(dots$main,  obsname[ObsWhich] , np)

    isylim <- !is.null(dots$ylim)
    yylim   <- expanddotslist(dots$ylim, np)
    
    isxlim <- !is.null(dots$xlim)
    xxlim   <- expanddotslist(dots$xlim, np)
    
    # point parameters
    ip <- !names(dots) %in% plotnames
    dotpoints <- dots[ip]
    dotpoints <- setdots(dotpoints, nx)   # expand all dots to nx values

    # these are different from default
    dotpoints$pch  <- expanddots(dots$pch, 1:nx,      nx)
    dotpoints$col  <- expanddots(dots$col, 1:nx,      nx)
    dotpoints$bg   <- expanddots(dots$bg,  1:nx,      nx)

    xyswap <- rep(xyswap, length  = np)
## for each variable
    iobs <- 0
    for (ip in 1:np) {
      io <- ObsWhich[ip] 

      # plotting parameters for deSolve output 1 (opens a plot)
      Dotmain   <- extractdots(dotmain, ip)
      Dotpoints <- extractdots(dotpoints, 1)

      Xlog <- Ylog <- FALSE
      if (! is.null(Dotmain$log)) { 
        Ylog  <- length(grep("y", Dotmain$log))
        Xlog  <- length(grep("x", Dotmain$log))
      }       

      if (! xyswap[ip]) {            
         if (!isylim)
            Dotmain$ylim <- Range(NULL, obs[,io], Ylog) 
         else 
            Dotmain$ylim <- yylim[[ip]]           
          
         Dotmain$xlim <- xxlim[[ip]]  
         if (is.null(Dotmain$xlim)) 
           Dotmain$xlim <- Range(NULL, obs[,1], Xlog)
         
         i.obs <- obs.pos[1,1]:obs.pos[1,2]
         if (length (i.obs) > 0)  
           do.call("plot", c(alist(obs[i.obs,1], obs[i.obs, io]), Dotmain, Dotpoints))
         else
           do.call("plot", c(NULL, Dotmain, Dotpoints))
      
         if (nx > 1) 
           for (j in 2:nx)  
             if (length (i.obs <- obs.pos[j,1]:obs.pos[j,2]) > 0)  
               do.call("points", c(alist(obs[i.obs,1], obs[i.obs,io]),  
                   extractdots(dotpoints, j)))

      } else {  # xyswap
          if (Dotmain$xlab =="x" & Dotmain$ylab == "") {
            xl <- Dotmain$ylab
            Dotmain$ylab <- Dotmain$xlab
            Dotmain$xlab <- xl
          }

         if (!isxlim) {
           xrange <- Range(NULL, obs[,io], Xlog)
           Dotmain$xlim <- xrange
         } else 
           Dotmain$xlim <- xxlim[[ip]]
          
         if (is.null(yylim[[ip]])) {
            Dotmain$ylim <- rev(Range(NULL, range(obs[,1]), Ylog))    # y-axis reversed
         } else
            Dotmain$ylim <- yylim[[ip]]  
            
            i.obs <- obs.pos[1,1]:obs.pos[1,2]
            if (length (i.obs) > 0) 
              do.call("plot", c(alist(obs[i.obs , io], obs[i.obs ,1]), Dotmain, Dotpoints))
             else
              do.call("plot", c(NULL, Dotmain, Dotpoints))
             
            if (nx > 1) 
              for (j in 2:nx) {
                i.obs <- obs.pos[j,1]:obs.pos[j,2]
                if (length (i.obs) > 0)  {
                  Dotpoints <- extractdots(dotpoints, j)
                  do.call("points", c(alist(obs[i.obs , io], obs[i.obs ,1]), Dotpoints))
                }  
              }

      }
   }
}
