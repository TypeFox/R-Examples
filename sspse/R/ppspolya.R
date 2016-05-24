#' @keywords internal
ppspolya <- function(s, n=length(s), N=2*n, c=1, wts=1/s, 
                     y=1:n, seed=NULL) {
      n <- length(s)
      wts <- c(c*wts,rep(mean(wts), N-n))
      if (n > length(y))
            stop("length(y) != n")
      if (n > N)
            stop("length(s) > N")
      if (length(wts) != N)
            stop("length(wts) < N")
      if (any(wts < 0) || all(wts <= 0))
            stop("wts are foobar")
      N <- as.integer(N)
      if (N <= 0)
            stop("N must be positive integer")
      if(!is.null(seed))  set.seed(as.integer(seed))
      out<-.C("ppspolya",
               y=as.double(c(y, rep(0, N-n))),
               size=as.double(c(s, rep(0, N-n))),
               w=as.double(cumsum(wts)),
               nin=as.integer(n),
               Nin=as.integer(N),
               PACKAGE="sspse")
       return(out)
}
