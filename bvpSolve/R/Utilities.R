### ============================================================================
### S3 methods
### ============================================================================

## An approximation function for bvpSolve objects that were solved with \
## bvpcol
approx <- function (x, ...) UseMethod("approx")

approx.default <- function (x, ...) {
if ("bvpSolve" %in% class (x))
  approx.bvpSolve(x,...)
else  
  stats::approx(x,...)
#nextmethod()  
}

### ============================================================================

approx.bvpSolve <- function(x, xout=NULL, ...){
 
  Attr <- attributes(x)
  if (Attr$name != "bvpcol")
    stop("can only use 'approx.bvpSolve' if problem was solved with 'bvpcol'")
  istate <- Attr$istate[-(1:6)]   # first 6 elements have nothing to do with continuation
  rstate <- Attr$rstate
  il <- length(xout)
  if (il <= 0)
    stop ("'approx' requires at least one value to approximate")   

  z <- rep(1, istate[4])    #### WAS 4
  if (Attr$colmod)
    appone <- function(x)
      .Fortran("mappsln", as.double(x),
            result = as.double(z), as.double(rstate), as.integer(istate))$result
  else if (! Attr$bspline)
    appone <- function(x)
      .Fortran("appsln", as.double(x),
            result = as.double(z), as.double(rstate), as.integer(istate))$result
  else 
    appone <- function(x)
      .Fortran("sysappsln", as.double(x), 
            result = as.double(z), as.double(rstate), as.integer(istate))$result
  
  Out <- NULL
  for (i in 1:il)
    Out <- rbind(Out,c(xout[i],appone(xout[i])))
  colnames(Out) <- colnames(x)
  class (Out) <- c( "bvpSolve", "matrix" ) 
  attr(Out, "name") <-  "approx"
  dimnames(Out) <- dimnames(x)
  Out
}

### ============================================================================

print.bvpSolve <- function(x, ...)
   print(as.data.frame(x), ... )

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

selectvar <- function (which, var, NAallowed = FALSE) {

    if (!is.numeric(which)) {
        ln <- length(which)
        # keep ordering...
        Select <- NULL
        for ( i in 1:ln) {
          ss <- which(which[i]==var)
          if (length(ss) ==0 & ! NAallowed) 
            stop("variable ", which[i], " not in var")
          else if (length(ss) == 0)
            Select <- c(Select,NA)
          else
            Select <- c(Select,ss)
        }
    }
    else {
        Select <- which + 1  # "Select now refers to the column number
        if (max(Select) > length(var))
            stop("index in 'which' too large")
        if (min(Select) < 1)
            stop("index in 'which' should be > 0")
    }
  return(Select)
}

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

### ============================================================================
### Plotting bvpSolve objects
### ============================================================================

plot.bvpSolve <- function (x, ..., which = NULL, ask = NULL, obs = NULL, 
    obspar = list()) {

## check observed data - can be a list
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
       if (is.character(obs[,1]) | is.factor(obs[,1]))   # long format - convert
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
    varnames <- colnames(x)
    Which   <- which 
    
    if (is.null(Which) & is.null(obs))  # All variables plotted
      Which <- 1 : (ncol(x)-1)
      
    else if (is.null(Which)) {          # All common variables in x and obs plotted
     Which <- which(varnames %in% obsname)
     Which <- Which[Which != 1]         # remove first element (x-value)
     Which <- varnames[Which]           # names rather than numbers
    } 

## Position of variables to be plotted in "x" 
    t       <- 1     # column with "times" 
    xWhich  <- selectvar(Which, varnames)
    np      <- length(xWhich)

    ldots   <- list(...)
    ndots   <- names(ldots)

## number of figures in a row and interactively wait if remaining figures
    ask <- setplotpar(ndots, ldots, np, ask)
    if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }

## create two lists: x2:   other bvpSolve objects, 
##                   dots: remaining (plotting) parameters
    x2     <- list()
    dots   <- list()
    nd     <- 0
    nother <- 0            

    if (length(ldots) > 0) 
     for ( i in 1:length(ldots))
      if ("bvpSolve" %in% class(ldots[[i]])) {
        x2[[nother <- nother + 1]] <- ldots[[i]]  
        names(x2)[nother] <- ndots[i]
      } else if (! is.null(ldots[[i]])) {
        dots[[nd <- nd+1]] <- ldots[[i]]
        names(dots)[nd] <- ndots[i]
      }
    nmdots <- names(dots)

## check compatibility of all bvpSolve objects    
    if (nother > 0) {
      for ( i in 1:nother) {            
        if (min(colnames(x2[[i]]) == varnames) == 0)
          stop("'x' is not compatible with other bvpSolve objects - colnames not the same")
      }
    } 

    nx <- nother + 1 # total number of bvpSolve objects to be plotted

## Position of variables in "obs" (NA = not observed)
    if (nobs > 0) {
      ObsWhich <- selectvar(varnames[xWhich], obsname, NAallowed = TRUE)
      ObsWhich [ ObsWhich > ncol(obs)] <- NA  
      Obspar <- setdots(obspar, nobs)
    } else 
      ObsWhich <- rep(NA, np)

## plotting parameters : split in plot parameters and point parameters
    plotnames <- c("xlab","ylab","xlim","ylim","main","sub","log","asp",
                   "ann","axes","frame.plot","panel.first","panel.last",
                   "cex.lab","cex.axis","cex.main")    
    
    # plot.default parameters
    ii <- names(dots) %in% plotnames
    dotmain <- dots[ii]
    dotmain <- setdots(dotmain, np)  # expand to np for each plot

    # these are different from the default
    dotmain$xlab <- expanddots(dots$xlab, varnames[t]     , np)
    dotmain$ylab <- expanddots(dots$ylab, ""              , np)
    dotmain$main <- expanddots(dots$main, varnames[xWhich], np)

    # ylim and xlim can be lists and are at least two values
    yylim  <- expanddotslist(dots$ylim, np)
    xxlim  <- expanddotslist(dots$xlim, np)

    # point parameters
    ip <- !names(dots) %in% plotnames
    dotpoints <- dots[ip]
    dotpoints <- setdots(dotpoints, nx)   # expand all dots to nx values

    # these are different from default
    dotpoints$type <- expanddots(dots$type, "l", nx)
    dotpoints$lty  <- expanddots(dots$lty, 1:nx, nx)
    dotpoints$pch  <- expanddots(dots$pch, 1:nx, nx)
    dotpoints$col  <- expanddots(dots$col, 1:nx, nx)
    dotpoints$bg   <- expanddots(dots$bg,  1:nx, nx)
 
## for each output variable (plot)
    iobs <- 0
    for (i in 1 : np) {
      ii <- xWhich[i]     # position of variable in 'x'
      io <- ObsWhich[i]   # position of variable in 'obs'
      
      # plotting parameters for bvpSolve output 1 (opens a plot)
      Dotmain   <- extractdots(dotmain, i)
      Dotpoints <- extractdots(dotpoints, 1)
      
      Xlog <- Ylog <- FALSE
      if (! is.null(Dotmain$log)) { 
        Ylog  <- length(grep("y",Dotmain$log))
        Xlog  <- length(grep("x",Dotmain$log))
      }       
      
      # ranges
      if ( is.null (yylim[[i]])) {
        yrange <- Range(NULL, x[, ii], Ylog)
        if (nother>0) 
         for (j in 1:nother) 
           yrange <- Range(yrange, x2[[j]][,ii], Ylog)
        if (! is.na(io)) yrange <- Range(yrange, obs[,io], Ylog)
          Dotmain$ylim <- yrange
      } else  
        Dotmain$ylim  <- yylim[[i]]
       

      if ( is.null (xxlim[[i]])) {
        xrange <- Range(NULL, x[, t], Xlog)
        if (nother>0) 
         for (j in 1:nother) 
           xrange <- Range(xrange, x2[[j]][,t], Xlog)
        if (! is.na(io)) xrange <- Range(xrange, obs[,1], Xlog)
          Dotmain$xlim <- xrange
      } else  
        Dotmain$xlim  <- xxlim[[i]]
      
      # first bvpSolve object plotted (new plot created)
      do.call("plot", c(alist(x[, t], x[, ii]), Dotmain, Dotpoints))
      
      # if other bvpSolve outputs
      if (nother > 0) 
        for (j in 2:nx)   
          do.call("lines", c(alist(x2[[j-1]][, t], x2[[j-1]][, ii]), 
                  extractdots(dotpoints, j)) )
         
      # if observed variables: select correct pars
      if (! is.na(io))   
           for (j in 1: nobs) 
              if (length (i.obs <- obs.pos[j, 1]:obs.pos[j, 2]) > 0) 
                do.call("points", c(alist(obs[i.obs, 1], obs[i.obs, io]), 
                         extractdots(Obspar, j) ))        
            
    }
}


### ============================================================================

diagnostics.bvpSolve<- function(obj, ...) {
  if (!"bvpSolve" %in% class(obj)) return(NULL)
  Attr <- attributes(obj)
  istate <- Attr$istate
  rstate <- Attr$rstate
  idid <- istate[1]

  if (is.null(istate) || is.null (rstate)) return(NULL)
  cat("\n--------------------\n")
  cat(paste( "solved with ",Attr$name))
  cat("\n--------------------\n")

  if (Attr$name == "bvpshoot") {
    cat("\n---------------------------------------------\n")
    cat("diagnostics of BVP solver ")
    cat("\n---------------------------------------------\n")
    df <- c( "The number of function evaluations              :", 
             "The number of jacobian evaluations +LU decomp   :",	
             "The number of steps                             :",
             "The number of calls to the ivp solver           :")
    printmessage(df, Attr$istate2)

    cat("\n---------------------------------------------\n")
    cat("diagnostics of the last run of the IVP solver ")
    cat("\n---------------------------------------------\n")
      diagnostics.deSolve(obj)
    
  } else if (Attr$name == "bvptwp") {
  
    if (idid ==0)  cat("  Integration was successful.\n") else
			{if (idid < 0 &  (Attr$acdc == FALSE))	cat("  Integration was successful but conditioning parameters NOT stabilized.\n")  else
       cat("  Integration was NOT successful\n")}
    df <- c( "The return code                               :",   #1
             "The number of function evaluations            :",
             "The number of jacobian evaluations            :",
             "The number of boundary evaluations            :",
             "The number of boundary jacobian evaluations   :",
             "The number of steps                           :",
             "The number of mesh resets                     :",
             "The maximal number of mesh points             :",   #2
             "The actual number of mesh points              :",
             "The size of the real work array               :",
             "The size of the integer work array            :"
             )

    printmessage(df, istate[c(1:7,10:13)])

    cat("\n--------------------\n")
    cat(paste( "conditioning pars"))
    cat("\n--------------------\n")

    df <- c( "kappa1  :",
             "gamma1  :",
             "sigma   :",
             "kappa   :",
             "kappa2  :")
    printmessage(df, rstate)
    
    if (Attr$acdc == TRUE)
      cat(paste("    The problem was solved for final eps equal to :",
               Attr$eps[1]," \n"))
  } else if (Attr$name == "bvpcol") {
    if (idid ==1)  cat("  Integration was successful.\n") else
       cat("  Integration was NOT successful\n")
    if (! Attr$colmod)  {
     df <- c( "The return code                                   :",   #1
              "The number of function evaluations                :",
              "The number of jacobian evaluations                :",
              "The number of boundary evaluations                :",
              "The number of boundary jacobian evaluations       :",
              "The number of steps                               :",
              "The actual number of mesh points                  :",
              "The number of collocation points per subinterval  :",
              "The number of equations                           :",
              "The number of components (variables)              :")
       printmessage(df, istate[1:10])
    } else {
     df <- c( "The return code                                   :",   #1
              "The number of function evaluations                :",
              "The number of jacobian evaluations                :",
              "The number of boundary evaluations                :",
              "The number of boundary jacobian evaluations       :",
              "The number of continuation steps                  :",
              "The number of succesfull continuation steps       :",
              "The actual number of mesh points                  :",
              "The number of collocation points per subinterval  :",
              "The number of equations                           :",
              "The number of components (variables)              :")
       printmessage(df, c(idid,Attr$icount[1:6],istate[7:10]))
    cat(paste("    The problem was solved for final eps equal to     :",Attr$eps[1]," \n"))
    }
  }

}

### ============================================================================
## internal helper functions for printing solver return code messages
## these functions are not exported

printmessage <-function(message1, state, message2 = NULL, Nr = 1:length(message1)) {
  if (is.null(message2)) {
    cat("\n", paste(formatC(Nr, "##", width = 2), message1,
              signif(state, digits = getOption("digits")), "\n"), "\n")
  } else {
    cat("\n", paste(formatC(Nr, "##", width = 2), message1,
              signif(state, digits = getOption("digits")), message2, "\n"), "\n")

  }
}
