
normgibbs<-function (N, n, a, b, cc, d, xbar, ssquared) 
{
  mat = matrix(ncol = 2, nrow = N)
  mu = cc
  tau = a/b
  mat[1, ] = c(mu, tau)
  for (i in 2:N) {
    muprec = n * tau + d
    mumean = (d * cc + n * tau * xbar)/muprec
    mu = rnorm(1, mumean, sqrt(1/muprec))
    taub = b + 0.5 * ((n - 1) * ssquared + n * (xbar - mu)^2)
    tau = rgamma(1, a + n/2, taub)
    mat[i, ] = c(mu, tau)
  }
  colnames(mat)=c("mu","tau")
  mat
}


# eof

