##################################################
#### AUTHOR:     Arnost Komarek               ####
####             (2005)                       ####
####                                          ####
#### FILE:       marginal.bayesGspline.R      ####
####                                          ####
#### FUNCTIONS:  marginal.bayesGspline        ####
####             plot.marginal.bayesGspline   ####
####                                          ####
##################################################

### ======================================
### marginal.bayes.Gspline
### ======================================
##
## 27/11/2005
marginal.bayesGspline <- function(dir = getwd(), extens = "", K, grid1, grid2, skip = 0, by = 1, last.iter, nwrite, only.aver = TRUE)
{
  thispackage = "bayesSurv"
  #thispackage = NULL
  standard <- FALSE
  dim <- 2
  n.in.mixmoment <- 1 + 2 + 3   ## number of columns in mixmoment.sim file
  if (length(K) != dim) stop(paste("K must be of length ", dim, sep=""))
  l.Gspl <- 2*K + 1
  KK.max <- l.Gspl[1]*l.Gspl[2]
  
  if (missing(grid1)){
    stop("grid1 must be given")
  }    
  else{
    nx1 <- length(grid1)
    if (nx1 <= 0) stop("grid1 must have a positive length")
  }    
  if (missing(grid2)){
    stop("grid2 must be given")
  }    
  else{
    nx2 <- length(grid2)
    if (nx2 <= 0) stop("grid2 must have a positive length")    
  }
  
  ## Check whether needed files are available
  ## * further,  whether at least first row has correct number of elements
  ## * and determine the MC sample size
  ## =======================================================================
  filesindir <- dir(dir)                                                      ## character vector with available files
  if (!length(filesindir)) stop("Empty directory with simulated values?")  
  if (sum(!is.na(match(filesindir, paste("mweight", extens, ".sim", sep=""))))){
    mix <- read.table(paste(dir, "/mweight", extens, ".sim", sep = ""), nrows = 1)
    k.max <- length(mix)
    if (k.max != KK.max) stop("Strange number of columns in mweight, probably wrong 'K' argument")    
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
    mix <- read.table(paste(dir, "/gspline", extens, ".sim", sep = ""), nrows = 2, header=TRUE)
    lmix <- ncol(mix)
    if (lmix != 5*dim) stop("This function is implemented only for bivariate case, so that gspline.sim must have 10 columns")
  }
  else{
    stop("File with simulated values of gamma/sigma/delta/intercept/scale not found.")
  }  

  if (sum(!is.na(match(filesindir, paste("mixmoment", extens, ".sim", sep=""))))){
    name.mix <- scan(paste(dir, "/mixmoment", extens, ".sim", sep = ""), nlines=1, what=character(0), quiet=TRUE)
    if (length(name.mix) != n.in.mixmoment) stop(paste("File mixmoment.sim should have ", n.in.mixmoment, " columns", sep=""))
    mix <- as.data.frame(matrix(scan(paste(dir, "/mixmoment", extens, ".sim", sep = ""), skip=1, quiet=TRUE), byrow=TRUE, ncol=n.in.mixmoment))
    colnames(mix) <- name.mix
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
    
  if (missing(skip)) skip <- 0
  else{
    if (skip > M) stop("You ask to skip more rows from the file than available.")
    if (skip < 0) skip <- 0
  }
  if (missing(by)) by <- 1
  else{
    if (by <= 0) by <- 1
  }    

  lvalue1 <- ifelse(only.aver, nx1, nx1* (1 + (M-skip-1) %/% by))
  lvalue2 <- ifelse(only.aver, nx2, nx2* (1 + (M-skip-1) %/% by))  
  if (missing(nwrite)) nwrite <- M

  mcmc <- .C("marginal_bayesGspline",
                average1 = double(nx1),    average2 = double(nx2),
                value1 = double(lvalue1),  value2 = double(lvalue2),
                M.now = integer(1),        as.integer(only.aver),
                as.character(dir),         as.character(extens),
                as.double(grid1),          as.double(grid2),             
                as.integer(K),
                as.integer(M),             as.integer(skip),         as.integer(by),   as.integer(nwrite),
                nx1 = as.integer(nx1),     nx2 = as.integer(nx2),
                err = integer(1),
               PACKAGE = thispackage)
  
  if (mcmc$err) stop("No results produced, something is wrong.")  

  if (!only.aver){
    val1 <- matrix(mcmc$value1, nrow = mcmc$nx1, ncol = mcmc$M.now)
    rownames(val1) <- paste(grid1)
    colnames(val1) <- 1:mcmc$M.now

    val2 <- matrix(mcmc$value2, nrow = mcmc$nx2, ncol = mcmc$M.now)
    rownames(val2) <- paste(grid2)
    colnames(val2) <- 1:mcmc$M.now    
  }
  aver1 <- mcmc$average1
  to.return1 <- data.frame(grid = grid1, average = aver1)
  rownames(to.return1) <- 1:mcmc$nx1

  aver2 <- mcmc$average2
  to.return2 <- data.frame(grid = grid2, average = aver2)
  rownames(to.return2) <- 1:mcmc$nx2

  RET <- list(margin1=to.return1, margin2=to.return2)
  attr(RET, "sample.size") <- mcmc$M.now
  if (!only.aver){
    attr(RET, "sample1") <- val1
    attr(RET, "sample2") <- val2
  }  
  
  class(RET) <- "marginal.bayesGspline"
  return(RET)  
}


### ======================================
### plot.marginal.bayesGspline
### ======================================
plot.marginal.bayesGspline <- function(x, type = "l", lty=1, bty = "n", xlab1, ylab1, main1, xlab2, ylab2, main2, sub, ...)
{
  dim <- 2

  if (missing(xlab1)) xlab1 <- "y1"
  if (missing(ylab1)) ylab1 <- "g(y1)"
  if (missing(xlab2)) xlab2 <- "y2"
  if (missing(ylab2)) ylab2 <- "g(y2)"  
  if (missing(main1)) main1 <- "McMC average - margin 1"
  if (missing(main2)) main2 <- "McMC average - margin 2"  
  if (missing(sub)) sub <- paste("M = ", attr(x, "sample.size"), sep="")
  oldpar <- par(mfrow=c(2, 1), bty=bty)
  on.exit(par(oldpar))
  plot(x$margin1$grid, x$margin1$average, type=type, lty=lty, xlab=xlab1, ylab=ylab1, ...)
  title(main=main1, sub=sub)
  plot(x$margin2$grid, x$margin2$average, type=type, lty=lty, xlab=xlab2, ylab=ylab2, ...)
  title(main=main2, sub=sub)

  return(invisible(x))
}  
                              




