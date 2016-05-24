binom.blaker.VHadj.acc <- function(x,n,p,type=c("orig","unimod"),acc.tol=1e-10,nmax=n+1000,int.eps=1e-12,...) {
  if (n < 1 || x < 0 || x > n) stop("Parameters n = ",n,", x = ",x, " wrong!")
  if (acc.tol <= 0) stop("Numerical tolerance ",acc.tol," nonpositive!")
  type <- match.arg(type)
  if (type != "orig" && type != "unimod") stop("Unknown type ",type,"!")
  m <- length(p)
  if (m < 1) stop("Empty vector p!")
  acc <- binom.blaker.acc(x,n,p,type=type,acc.tol=acc.tol,...=...)
  nn <- n+1
  ind.susp <- 1:length(p)
  ind.susp1 <- ind.susp[p[ind.susp] <= x/n]
  ind.susp2 <- ind.susp[p[ind.susp] >  x/n]
  repeat {
    if (nn > nmax) {
      warning("n = ", n, ", x = ", x,": Upper limit of ", nmax, 
              " too low, Vos-Hudson adjustment may be incomplete.") 
      break
    }
    xx <- x*nn/n
    aa1 <- pmin(2*pbeta(p[ind.susp1],xx,nn-xx+1),1)
    aa2 <- pmin(2*pbeta(1-p[ind.susp2],nn-xx,xx+1),1)
    ind.susp1 <- ind.susp1[aa1 > acc[ind.susp1]]    
    ind.susp2 <- ind.susp2[aa2 > acc[ind.susp2]]    
    if (length(ind.susp1) + length(ind.susp2) < 1) break
    acc.nn.L <- numeric(0)
    acc.nn.R <- numeric(0)
    if (length(ind.susp1) > 0) {
      x.star <- ceiling(xx-int.eps)
      acc.nn.L <- binom.blaker.acc(x.star,nn,p[ind.susp1],type=type,acc.tol=acc.tol,...=...)
    }
    if (length(ind.susp2) > 0) {
      x.star <- floor(xx+int.eps)
      acc.nn.R <- binom.blaker.acc(x.star,nn,p[ind.susp2],type=type,acc.tol=acc.tol,...=...)
    }
    acc.nn <- c(acc.nn.L,acc.nn.R)
    ind.susp <- c(ind.susp1,ind.susp2)
    acc[ind.susp] <- pmax(acc[ind.susp],acc.nn)
    nn <- nn+1
  }

  acc
}
