###########################################
#### AUTHOR:     Arnost Komarek        ####
####             (2005)                ####
####                                   ####
#### FILE:       bayesGspline.R        ####
####                                   ####
#### FUNCTIONS:  bayesGspline          ####
####             plot.bayesGspline     ####
###########################################

### ======================================
### bayesGspline
### ======================================
##
## 02/11/2004
## 02/02/2005
bayesGspline <- function(dir = getwd(),
                         extens = "",
                         extens.adjust= "_b",
                         grid1,
                         grid2,
                         skip = 0,
                         by = 1,
                         last.iter,                         
                         nwrite,
                         only.aver = TRUE,
                         standard = FALSE,
                         version = 0)
{
  thispackage = "bayesSurv"
  #thispackage = NULL

  if (version < 30) standard <- FALSE
  
  if (missing(grid1)){
    stop("grid1 must be given")
  }    
  else{
    nx1 <- length(grid1)
    if (nx1 <= 0) stop("grid1 must have a positive length")
  }    
  if (missing(grid2)){
    dim <- 1
    grid2 <- 0
    nx2 <- 0
  }    
  else{
    dim <- 2
    nx2 <- length(grid2)
    if (nx2 <= 0) stop("grid2 must have a positive length")    
  }
  ngrid <- ifelse(dim == 1, nx1, nx1*nx2)
  
  ## Check whether needed files are available
  ## * further,  whether at least first row has correct number of elements
  ## * and determine the MC sample size
  ## =======================================================================
  filesindir <- dir(dir)                                                      ## character vector with available files
  if (!length(filesindir)) stop("Empty directory with simulated values?")  
  if (sum(!is.na(match(filesindir, paste("mweight", extens, ".sim", sep=""))))){
    mix <- read.table(paste(dir, "/mweight", extens, ".sim", sep = ""), nrows = 1)
    k.max <- length(mix)
  }
  else{
    stop("File with simulated values of mixture weights not found.")
  }  

  if (sum(!is.na(match(filesindir, paste("mmean", extens, ".sim", sep=""))))){
    mix <- read.table(paste(dir, "/mmean", extens, ".sim", sep = ""), nrows = 1)
    kmax2 <- length(mix)/dim
    if (k.max != kmax2) stop("Different total_length of the G-spline indicated by files mweight.sim and mmean.sim.")
  }
  else{
    stop("File with simulated values of mixture means indeces not found.")
  }    
   
  if (sum(!is.na(match(filesindir, paste("gspline", extens, ".sim", sep=""))))){
    mix <- read.table(paste(dir, "/gspline", extens, ".sim", sep = ""), nrows = 1)     
    lmix <- length(mix)
    if (lmix != 5*dim) stop(paste("You indicate that dimension is ", dim, " so that file gspline.sim must have ", 5*dim, " columns", sep=""))
  }
  else{
    stop("File with simulated values of gamma/sigma/delta/intercept/scale not found.")
  }  

  if (sum(!is.na(match(filesindir, paste("mixmoment", extens, ".sim", sep=""))))){
    mix <- read.table(paste(dir, "/mixmoment", extens, ".sim", sep = ""), header = TRUE)
    if (missing(last.iter)) M <- dim(mix)[1]
    else{
      M <- last.iter
      if (last.iter > dim(mix)[1]) M <- dim(mix)[1]
      if (last.iter <= 0)          M <- dim(mix)[1]
    }      
  }
  else{
    stop("File mixmoment.sim not found.")
  }

  if (version >= 30){
    if (version != 30 & version != 31 & version != 32){
      stop("version argument must be either 30 or 31 or 32")
    }
    if (version == 30 | version == 31){
      if (!sum(!is.na(match(filesindir, paste("mixmoment", extens.adjust, ".sim", sep=""))))){
        stop(paste("File mixmoment", extens.adjust, ".sim not found.", sep=""))      
      }
    }
  }
    
  if (missing(skip)) skip <- 0
  else{
    if (skip > M) stop("You ask to skip more rows from the file than available.")
    if (skip < 0) skip <- 0
  }
  if (missing(by)) by <- 1
  else{
    if (by <= 0) by <- 1
  }    

  lvalue <- ifelse(only.aver, ngrid, ngrid* (1 + (M-skip-1) %/% by))
  if (missing(nwrite)) nwrite <- M
  
  mcmc <- .C("bayesGspline", average = double(ngrid), value = double(lvalue),           M.now = integer(1),     as.integer(only.aver),
                             as.character(dir),       as.character(extens),             as.character(extens.adjust),
                             x1 = as.double(grid1),   x2 = as.double(grid2),
                             as.integer(k.max),       as.integer(M),                    as.integer(skip),       as.integer(by),
                             as.integer(nwrite),      nx1 = as.integer(nx1),            nx2 = as.integer(nx2),
                             as.integer(version),     standard = as.integer(standard),  err = integer(1),
               PACKAGE = thispackage)

  if (mcmc$err) stop("No results produced, something is wrong.")
  
  if (dim == 1){
    if (!only.aver){
      val <- matrix(mcmc$value, nrow = mcmc$nx1, ncol = mcmc$M.now)
      rownames(val) <- paste(grid1)
      colnames(val) <- 1:mcmc$M.now
    }
    aver <- mcmc$average
    to.return <- data.frame(grid = grid1, average = aver)
    rownames(to.return) <- 1:mcmc$nx1
  }
  else{     ## dim == 2
    if (!only.aver){
      val <- matrix(mcmc$value, nrow = ngrid, ncol = mcmc$M.now)
      rnames1 <- rep(grid1, mcmc$nx2)
      rnames2 <- rep(grid2, rep(mcmc$nx1, mcmc$nx2))
      rnames <- paste(rnames1, ", ", rnames2, sep = "")      
      rownames(val) <- rnames
      colnames(val) <- 1:mcmc$M.now
    }      
    aver <- matrix(mcmc$average, nrow = mcmc$nx1, ncol = mcmc$nx2)
    to.return <-  list(grid1 = grid1, grid2 = grid2, average = aver)
  }

  attr(to.return, "sample.size") <- mcmc$M.now
  if (!only.aver) attr(to.return, "sample") <- val
  
  class(to.return) <- "bayesGspline"
  return(to.return)
}


### ======================================
### plot.bayesGspline
### ======================================
plot.bayesGspline <- function(x,
                              add = FALSE,
                              type = "l",
                              lty=1,                         
                              bty = "n",
                              xlab,
                              ylab,
                              main,
                              sub,
                              ...)
{

  if (is.null(x$grid2)) dim <- 1
  else                  dim <- 2

  if (missing(main)) main <- "McMC average of the density"
  if (missing(sub)) sub <- paste("M = ", attr(x, "sample.size"), sep="")
  
  if (dim == 1){
    if (missing(xlab)) xlab <- "x"
    if (missing(ylab)) ylab <- "g(x)"
    if (add){ lines(x$grid, x$average, lty=lty, ...) }
    else{     plot(x$grid, x$average, type=type, bty=bty, lty=lty, xlab=xlab, ylab=ylab, ...)
              title(main=main, sub=sub)
        }              
  }
  else{
    if (missing(xlab)) xlab <- "x1"
    if (missing(ylab)) ylab <- "x2"
    if (add){ contour(x$grid1, x$grid2, x$average, lty=lty, add=TRUE, ...) }
    else{     contour(x$grid1, x$grid2, x$average, type=type, bty=bty, lty=lty, xlab=xlab, ylab=ylab, ...)
              title(main=main, sub=sub)
        }
  }    

  return(invisible(x))
}  
                              


