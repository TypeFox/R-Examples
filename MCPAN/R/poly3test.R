`poly3test` <-
function(time, status, f, type="Dunnett", cmat=NULL, method="BP", alternative="two.sided", dist="MVN", k=3, ...)

{

aargs<-list(...)

est <- poly3estf(time=time, status=status, f=f, tmax=max(time), method=method, k=k)

n <- est$n
varnames<-est$names
ngroups <- length(n)


# check the contrast matrix

    if (is.null(cmat)) {
      if(type=="Dunnett") {
        if(is.null(aargs$base)){base<-1}
        else{base<-aargs$base}
        cmat <- contrMat(n=n, type=type, base=base)
       }
       else{cmat <- contrMat(n = n, type = type)}
    }
    else {
        if (!is.matrix(cmat) || ncol(cmat) != ngroups)
         {stop("cmat must be a matrix with number of columns = number of groups")}
    }

# Check margin

out <- Waldtest(estp=est$estp, varp=est$varcor, cmat=cmat, alternative=alternative, dist=dist)

# estimates using the MLE

out$estimate <- cmat %*% est$estimate


colnames(cmat)<-varnames

out$time<-time
out$status<-status
out$f<-f
out$method<-method
out$cmat<-cmat
out$k<-k

out$sample.estimate<-est

class(out)<-"poly3test"

return(out)

}

