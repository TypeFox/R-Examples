TwoCop <-
function(x, y, Nsim=100, paired=FALSE, alpha=0.95)
{

  dimx  = dim(x)
  dimy =  dim(y)

  n1 = dimx[1]
  n2 = dimy[1]
  d  = dimx[2]
  d2 = dimy[2]

  if (d != d2) stop("Samples x and y must have the same number of dimensions!!")
  if (n1!=n2 & paired==TRUE) stop("Paired samples must have the same sample size!!")

  out0 = .C("pvalue2",
                  as.double(x),
                  as.double(y),
                  as.integer(n1),
                  as.integer(n2),
                  as.integer(d),
                  as.integer(Nsim),
                  as.integer(paired),
                  cvm    = double(1),
                  cvmhat = double(Nsim),
                  pvalue = double(1),
		  PACKAGE="TwoCop"
                  )

  VaR   = quantile(out0$cvmhat,alpha)
  list(pvalue = out0$pvalue, cvm = out0$cvm,VaR = VaR, cvmsim = out0$cvmhat)
}
