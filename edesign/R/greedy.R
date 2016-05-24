greedy <- function(A,nf,ns,etol=0,mattest=TRUE)
  {
    A <- as.matrix(A)
    na <- dim(A)[1]
    nf <- as.integer(nf)
    ns <- as.integer(ns)
    ne <- na-nf

    if (nf==0) stop("algorithm fails for nf=0")
    if (ns==0) stop("number of points have to be added is 0")
    if (ne==0) stop("number of eligible points is 0")

    if(mattest)
      {
        if(!is.double(eigen(A)$values))
          stop("Matrix isn't symmetric")
        if(sum(eigen(A)$values>etol)<na)
          stop("Matrix isn't positive definite")
      }


    ans<-.Fortran("grd",
                  A=as.double(A),
                  lda=as.integer(na),
                  na=as.integer(na),
                  nf=as.integer(nf),
                  ne=as.integer(ne),
                  ns=as.integer(ns),
                  S=integer(ne),
                  det=double(1),
                  integer(na), # ind
                  double(na*na), #As
                  ldas=as.integer(na),
                  double(na*na), #Bs
                  ldbs=as.integer(na),
                  double(na*na), #INV
                  ldinv=as.integer(na),
                  double(na), #v
                  double(na), #w
                  integer(1), #ierr
                  PACKAGE="edesign"
                  )

    ret <- list(S=ans$S,det=ans$det)
    ret$A <- A
    ret$na <- na
    ret$nf <- nf
    ret$ns <- ns
    ret$ne <- ne
    ret$method <- "greedy"
    class(ret) <- "monet"
    ret
  }
