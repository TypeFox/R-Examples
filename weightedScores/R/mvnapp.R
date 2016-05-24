# approximation to MVN rectangle probabilities using methods
# in Joe (1995), JASA
mvnapp<- function(lb, ub, mu, sigma, type=1, eps = 1.e-05, nsim=0)
{
  #if(!is.loaded(symbol.C("mvnapp"))) dyn.load("./mvn.so")
  m <- length(ub)
  if(m>=8 && nsim==0) nsim<-10000
  if(m!=length(lb)) stop("lengths of w and x must be the same")
  tem <- sqrt(diag(sigma))
  w <- (lb - mu)/tem
  x <- (ub - mu)/tem
  tem <- diag(1/tem)
  corr <- tem %*% sigma %*% tem
  corr<-c(corr)
  if(eps<1.e-06) eps<- 1.e-06
  
  # print(w)
  # print(x)
  # print(corr)
  # print(inf)
  perr <- 0.
  # nsim<-0
  ifail<-0
  out <- .C("mvnapp",
    as.integer(type), as.integer(m), as.double(w),
    as.double(x), as.double(corr), as.integer(nsim),
    as.double(eps), prob = as.double(perr), perr = as.double(perr),
    ifail = as.integer(ifail))
  out <- list(pr = out$prob, esterr=out$perr, ifail = out$ifail)
  out
}
