#' Supporting function
#'
#' @importFrom stats toeplitz

acfpacf.pacf <-
function(x,m){

    nx=length(x)
    r<-acfpacf.acf(x,1)
    rho=r$acf
    rho=as.matrix(rho)

   pr<-matrix(0,1,m)
    pr[1]=rho[2]

   for (k in 2:m)
     {pmat=toeplitz(rho[1:k])
     rhovec=rho[2:(k+1)]
     phi=qr.solve(pmat)%*%rhovec
      pr[k]=phi[k] }
       pacf=pr

    result = list(pacf=pacf)
    class(result) = "acfpacf.pacf"
    result
}

