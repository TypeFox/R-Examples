# This is file ../spam/R/precmat.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     


# construct various precision matrices



precmat <- function(n, season=12, m=n, A=NULL, order=1, ... , type="RW1") {
  avtype <- c("rw1", "rw2", "rwn", "season","igmrfreglat","igmrfirreglat","gmrfreglat")
  method <- pmatch(tolower(type), avtype)
  if (is.na(method)) 
     stop("Precision matrix type not implemented yet. Please ask for.")
  switch(method,
         return(precmat.RW1(n)),
         return(precmat.RW2(n)),
         return(precmat.RWn(n, order)),
         return(precmat.season(n,season)),
         return(precmat.IGMRFreglat(n,m,...)),
         return(precmat.IGMRFirreglat(A,...)),
         return(precmat.GMRFreglat(n,m,...)))
}

precmat.IGMRFreglat <- function (n, m, order=1, anisotropy = 1)
{
  if((n<2)|(m<2))
    stop("n and m need to be >1")
  if(order==1) {
    if (anisotropy < 0 | anisotropy > 2) 
        stop("anisotropy parameter needs to be in [0,2]")
    return(kronecker(precmat.RW1(m), diag.spam(2-anisotropy, n)) + 
        kronecker(diag.spam(anisotropy, m), precmat.RW1(n)))
  } else if (order==2) {
    return(kronecker(precmat.RW2(m), diag.spam(n)) + 
        kronecker(diag.spam( m), precmat.RW2(n)))
  }
  if( (order<1)|(order>=min(n,m)))
        stop("order needs to be between 1 and min(n,m)-1")
  Dn <- diff.spam(diag.spam(n), lag = 1, differences = order) 
  Dm <- diff.spam(diag.spam(m), lag = 1, differences = order) 
  return(kronecker(t(Dm)%*%Dm, diag.spam(n)) +             
        kronecker(diag.spam( m), t(Dn)%*%Dn))
}

precmat.IGMRFirreglat <- function(A, eps= .Spam$eps) {
  if(!is.spam(A)) A <- as.spam(A, eps)
  A@entries <- rep.int(1, length(A))
  test <- isSymmetric.spam(A, tol = eps * 100)
  if (!isTRUE(test)) 
    stop("Input matrix  not symmetric (up to 100*eps)",call. = FALSE)
  
  return(diag.spam( diff(A@rowpointers)) - A)
}


precmat.RW1 <- function(n) {
  if(n<2)
    stop("Dimension 'n' should be larger than two")
  Q <- spam(0,n,n)
  Q@entries <- rep(-1,n-1)
  Q@colindices <- as.integer( seq.int(2, to=n,by=1))
  Q@rowpointers <- as.integer(c(seq.int(1,to=n,by=1),n))

  return(Q + t(Q) + diag.spam(c(1, rep.int(2, n-2), 1)))
}

precmat.RW2<- function(n) {
  if(n<4)
    stop("Dimension 'n' should be larger than three")
  Q <- spam(0,n,n)
  Q@entries <- c(-2,1,rep(c(-4,1),n-3),-2)
  Q@colindices <- as.integer( c(rep(2:(n-1),each=2)+c(0,1),n))
  Q@rowpointers <- as.integer(c(seq.int(1,to=2*(n-1),by=2),2*(n-1),2*(n-1)))


  return(Q + t(Q) + diag.spam(c(1,5, rep.int(6, n-4), 5,1)))
}

precmat.RWn <- function(n, order=3)
{
  if( (order<1)|(order>=n))
	stop("order needs to be between 1 and n-1")
  D <- diff.spam(diag.spam(n), lag=1, differences=order)	
  return( t(D)%*%D )
}


precmat.season <- function(n, season=12) {
  if(n<2*season)
    stop("Dimension 'n' should be larger than twice the season")
  Q <- spam(0,n,n)

  m <- season
  first <- rev( outer(1:(m-1), (m-1):1, "pmin"))
  mid <- rep.int((m-1):1, n-2*m+2)
  tmp <- outer(1:(m-1), 1:(m-1), "pmin")
  last <- rev( tmp[upper.tri(tmp)])
  
  Q@entries <- c( first, mid, last)
  Q@rowpointers <- as.integer( c( seq.int( 1, length.out=n-m+1, by=m-1),    # first and mid
                                 (n-m)*(m-1)+1+cumsum((m-1):1),
                                 (m-1)*(n-m/2)+1)  # last
                              )
  Q@colindices <- as.integer( c( rep.int(1:(m-1), n-m+1) +
                                rep.int(1:(n-m+1), rep.int(m-1, n-m+1)), # first and mid
                               n-last+1 )        # last
                             )
  return(Q + t(Q) + diag.spam(c(1:m, rep.int(m, n-2*m), m:1)))

}

