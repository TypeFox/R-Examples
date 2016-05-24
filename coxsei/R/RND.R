RND0 <-
function(n,int,tol=.Machine$double.eps^0.5,...)Quant(runif(n),int=int,tolerance=tol,...)
RND <- function(n,int,tol=.Machine$double.eps^0.5,
                epsabs=1e-10,epsrel=1e-10,limit=1000){
  .External("rndhaz",int,as.integer(n),new.env(),
            tol,epsabs,epsrel,as.integer(limit))
}


