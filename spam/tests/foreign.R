# This is file ../spam/tests/foreign.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     








# A few rudimentary tests to check transformations.
# This is for illustration and will not be run.

if (FALSE) {

  options( echo=FALSE)
  library( spam, warn.conflict=FALSE)

  test.for.zero <- function( xtest, xtrue, tol= 1.0e-6, tag=NULL){

    if( !is.null(tag)){
      cat( "testing: ", tag, fill=TRUE)}
    
    
    test.value <- sum( abs(c(xtest) - c( xtrue) ) )
    if(   test.value < tol ){
      cat("** PASSED test at tolerance ", tol, fill=TRUE)}
    else{ cat( "## FAILED test value = ", test.value, " at tolerance ", tol,
              fill=TRUE)}
    
  }


  xn <- 3
  xm <- 2
  set.seed(14)

  X <- array(runif(xn*xm), c( xn,xm))
  
  S <- as.spam(X)

  R <- as.matrix.csr.spam(S)
  test.for.zero(as.matrix(R),X)
  
  Q <- as.spam.matrix.csr(R)
  test.for.zero(Q,X)


  U <- as.dgRMatrix.spam(S)
  test.for.zero(as.matrix(U),S)
  
   
  V <- as.spam.dgRMatrix(U)
  test.for.zero(V,X)

  W <- as.dgCMatrix.spam(S)
  test.for.zero(as.matrix(W),S)
  
  Z <- as.spam.dgCMatrix(W)
  test.for.zero(Z,X)
  
  
  

  lundAs <- read.HB(system.file("external/lund_a.rsa",package = "Matrix"))
  lundAm <- readHB(system.file("external/lund_a.rsa",package = "Matrix"))
  test.for.zero(lundAs, as.matrix(lundAm))
  




  lundAs <- read.MM(system.file("external/lund_a.mtx",package = "Matrix"))
  lundAm <- readMM(system.file("external/lund_a.mtx",package = "Matrix"))
  test.for.zero(lundAs, as.matrix(lundAm))
  
  tmp <- read.MM(gzcon(url("ftp://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/bcspwr/bcspwr01.mtx.gz")))

  image(tmp <- read.MM(gzcon(url("ftp://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/acoust/young1c.mtx.gz"))))
  
  tmp <- read.MM(gzcon(url("ftp://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/bcspwr/bcspwr01.mtx.gz")))

  
  tmp <-read.MM(gzcon(url("ftp://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/platz/plskz362.mtx.gz")))

  test.for.zero(sum(upper.tri.spam(tmp) - t( lower.tri.spam(tmp))),0)
  test.for.zero(norm(tmp,'F'),8.152348)
  test.for.zero(dim(tmp),rep(362,2))

  
  
# use
  # file=gzcon(url("ftp://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/platz/plskz362.mtx.gz"))
  #open(file)



  
}

