oclSampling <- function(kernel, smat, alpha, beta=beta, D, By, k)
{
   N1 <- length(By)
   F1 <- ncol(smat)
   nsamp1 <- nrow(smat)
   logpostdens <- newbetaret <- rep(0,nsamp1 )
   out <- oclRun(kernel, 2*nsamp1,
      as.vector((t(smat))),
      as.vector(t(D)),
      By,
      alpha,
      beta,
      nsamp1,
      N1,
      F1,
      k)
   return( matrix(out,ncol=2) )
}


