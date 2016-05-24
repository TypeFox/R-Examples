##############################################
#### AUTHOR:     Arnost Komarek           ####
####             (2005)                   ####
####                                      ####
#### FILE:       sampled.kendall.tau.R    ####
####                                      ####
#### FUNCTIONS:  sampled.kendall.tau      ####
####                                      ####
##############################################

### ======================================
### sampled.kendall.tau
### ======================================
##
## 26/11/2005
sampled.kendall.tau <- function(dir = getwd(), extens = "", K, skip = 0, by = 1, last.iter, nwrite)
{
  thispackage = "bayesSurv"
  #thispackage = NULL
  
  dim <- 2                      ## only for bivariate G-splines
  n.in.mixmoment <- 1 + 2 + 3   ## number of columns in mixmoment.sim file
  if (length(K) != dim) stop(paste("K must be of length ", dim, sep=""))
  KK.max <- (2*K[1]+1) * (2*K[2]+1)
  
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
    if (mix$gamma1[1] != mix$gamma1[2]) stop("gamma1 in the file gspline.sim must be constant")
    if (mix$gamma2[1] != mix$gamma2[2]) stop("gamma2 in the file gspline.sim must be constant")
    if (mix$sigma1[1] != mix$sigma1[2]) stop("sigma1 in the file gspline.sim must be constant")
    if (mix$sigma2[1] != mix$sigma2[2]) stop("sigma2 in the file gspline.sim must be constant")
    if (mix$delta1[1] != mix$delta1[2]) stop("delta1 in the file gspline.sim must be constant")
    if (mix$delta2[1] != mix$delta2[2]) stop("delta2 in the file gspline.sim must be constant")
    GAMMA1 <- mix$gamma1[1]
    GAMMA2 <- mix$gamma2[1]    
    SIGMA1 <- mix$sigma1[1]
    SIGMA2 <- mix$sigma2[1]    
    DELTA1 <- mix$delta1[1]
    DELTA2 <- mix$delta2[1]
    IND1 <- (-K[1]):K[1]
    IND2 <- (-K[2]):K[2]        
    MU1 <- GAMMA1 + DELTA1*IND1
    MU2 <- GAMMA2 + DELTA2*IND2

    ### Basis for the computation of kendall's tau
    MU1diff <- (outer(MU1, MU1, "-"))/(sqrt(2)*SIGMA1)    
    MU2diff <- (outer(MU2, MU2, "-"))/(sqrt(2)*SIGMA2)
    pn.MU1diff <- pnorm(MU1diff)
    pn.MU2diff <- pnorm(MU2diff)
    l.MU1diff <- length(MU1diff)
    l.MU2diff <- length(MU2diff)    
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

  lvalue <- 1 + (M-skip-1) %/% by
  if (missing(nwrite)) nwrite <- M

  mcmc <- .C("sampledKendallTau", tau = double(lvalue),   M.now = integer(1),
                                  as.character(dir),      as.character(extens),
                                  as.integer(K),
                                  as.double(pn.MU1diff),  as.double(pn.MU2diff),
                                  as.integer(M),          as.integer(skip),       as.integer(by),  as.integer(nwrite),
                                  err = integer(1),
              PACKAGE = thispackage)             

  if (mcmc$err) stop("No results produced, something is wrong.")
  return(mcmc$tau[1:mcmc$M.now])
}

