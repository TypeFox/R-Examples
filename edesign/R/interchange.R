interchange <- function(A,nf,ns,S.start,etol=0,mattest=TRUE)
  {
    A <- as.matrix(A)
    na <- dim(A)[1]
    nf <- as.integer(nf)
    ns <- as.integer(ns)
    ne <- na-nf

    if ((nf==0)&(ns==1)) stop("algorithm fails for nf=0 and ns=1")
    if (ns==0) stop("number of points have to be added is 0")
    if (ne==0) stop("number of eligible points is 0")

    if(mattest)
      {
        if(!is.double(eigen(A)$values))
          stop("Matrix isn't symmetric")
        if(sum(eigen(A)$values>etol)<na)
          stop("Matrix isn't positive definite")
      }

    
    if(!is.null(S.start))
      {
        if(length(S.start)!=ne)
          {
            stop("Error in dimension of S.start")
          }
        if(length(S.start[S.start!=0])!=ns)
          {
            stop("S.start is infeasible")
          }
      } else {
        stop("S.start missing")
      }
    
    ans <- .Fortran("change",
                    A=as.double(A),
                    lda=as.integer(na),
                    na=as.integer(na),
                    nf=as.integer(nf),
                    ne=as.integer(ne),
                    ns=as.integer(ns),
                    S=as.integer(S.start),                    
                    det=double(1),
                    integer(na),
                    double(na*na),
                    ldas=as.integer(na),
                    double(na*na),
                    ldbs=as.integer(na),
                    double(na*na),
                    ldinv=as.integer(na),
                    double(na),
                    double(na),
                    integer(1),
                    PACKAGE="edesign"
                    )
    ret <- list(S.start=S.start,S=ans$S,det=ans$det)
    ret$A <- A
    ret$na <- na
    ret$nf <- nf
    ret$ns <- ns
    ret$ne <- ne
    ret$method <- "interchange"
    class(ret) <- "monet"
    ret
  }

