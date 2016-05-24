"CF" <- function(a, finite = FALSE, tol=0){
  ii <- is.infinite(a)
  if(any(ii)){
    n <- min(which(ii))
    a <- a[seq_len(n-1)]
    finite <- TRUE
  }
  b0 <- a[1]
  a <- a[-1]
  return(GCF(a=rep(1,length(a)), b=a, b0=b0, finite=finite, tol=tol))
}
 
"GCF" <- function(a, b, b0=0, finite=FALSE, tol=0){
  stopifnot(length(b) == length(a))
  stopifnot(length(b0) == 1)
  stopifnot(length(a) > 0)
  if(tol <= 0){tol <- .Machine$double.eps}
  n <- length(a)
  if(is.complex(a) | is.complex(b)){
    a <- as.complex(a)
    b <- as.complex(b)
    jj <- .C("contfrac_complex",
             as.double(Re(a)),
             as.double(Im(a)),
             as.double(Re(b)),
             as.double(Im(b)),
             as.integer(n),
             r = double(1),
             i = double(1),
             tol = double(1),
             PACKAGE = "contfrac")
    if(abs(jj$tol) <= tol | finite){
      return(b0 + jj$r + 1i*jj$i)
    } else {
      warning("Continued fraction (complex) not converged")
      return(NA)
    }
  } else {
    stopifnot(is.numeric(a))
    stopifnot(is.numeric(b))
    jj <-  .C("contfrac",
              a = as.double(as.vector(a)),
              b = as.double(as.vector(b)),
              n = as.integer(n),
              f = double(1),
              tol = double(1),
              PACKAGE = "contfrac"
             )
    if(abs(jj$tol) <= tol | finite){
      return(b0 + jj$f)
    } else {
      warning("Continued fraction not converged")
      return(NA)
    }
  }
}

"gconvergents" <- function(a,b,b0=0){
  n <- length(a)
  stopifnot(length(a) == length(b))
  stopifnot(length(a) > 0)
  stopifnot(length(b0) == 1)
  if(is.complex(a) | is.complex(b) | is.complex(b0)){
    a <- as.complex(a)
    b <- as.complex(b)
    jj <- .C("convergents_complex",
             as.double(as.vector(Re(a))),
             as.double(as.vector(Im(a))),
             as.double(as.vector(Re(b))),
             as.double(as.vector(Im(b))),
             as.double(Re(b0)),
             as.double(Im(b0)),
             as.integer(n),
             Ar = double(n+1),
             Ai = double(n+1),
             Br = double(n+1),
             Bi = double(n+1),
             PACKAGE = "contfrac"
             )
    return(list(A = jj$Ar + 1i*jj$Ai , B = jj$Br + 1i*jj$Bi))
  } else {
    stopifnot(is.numeric(a))
    stopifnot(is.numeric(b))
    jj <-  .C("convergents",
              as.double(as.vector(a)),
              as.double(as.vector(b)),
              as.double(b0),
              as.integer(n),
              A=double(n+1),
              B=double(n+1),
              PACKAGE = "contfrac"
              )
    return(list(A=jj$A , B=jj$B))
  }
}

"convergents" <- function(a){
  gconvergents(a=rep(1,length(a)-1),b=a[-1],b0=a[1])
}

"as_cf" <- function(x, n=10){
  stopifnot(length(x)==1)
  stopifnot(is.double(x))
  out <- double(n)

  for(i in seq_len(n)){
    jj <- floor(x)
    out[i] <- jj
    x <- 1/(x-jj)
  }
  out
}
