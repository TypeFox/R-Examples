nbmle <-
  function(x,mu,isize,w=rep(1,length(x)),n=sum(w),xsum=x){
    
    fun <- function(logsize,x,mu,w,n){
      size <- exp(logsize)
      fval <- (-1/n)*(sum(lgamma(size+x)) - n*lgamma(size) + sum(w*size*log(size/(size+mu))) + sum(xsum*log(mu/(size+mu))))
      return(fval)
    }
    size <- tryCatch(suppressWarnings(exp(nlm(f=fun,p=log(isize),x=x,mu=mu,w=w,n=n)$estimate)),error = function(e) isize)
    return(size)
    
  }