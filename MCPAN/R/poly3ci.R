`poly3ci` <-
function(time, status, f, type="Dunnett", cmat=NULL, method="BP", alternative="two.sided", conf.level=0.95, dist="MVN", k=3, ...)

{

paargs<-list(...)

est <- poly3estf(time=time, status=status, f=f, tmax=max(time), method=method, k=k)

n <- est$n
varnames<-est$names
ngroups <- length(n)

# check the contrast matrix

    if (is.null(cmat)) {
      if(type=="Dunnett") {
        if(is.null(paargs$base)){base<-1}
        else{base<-paargs$base}
        cmat <- contrMat(n=n, type=type, base=base)
       }
       else{cmat <- contrMat(n = n, type = type)}
    }
    else {
        if (!is.matrix(cmat) || ncol(cmat) != ngroups)
         {stop("cmat must be a matrix with number of columns = number of groups")}
    }

# Check margin

out <- Waldci(estp=est$estp, varp=est$varp, varcor=est$varcor, cmat=cmat, alternative=alternative, conf.level=conf.level, dist=dist)


out$estimate <- cmat %*% est$estimate


colnames(cmat)<-varnames

out$time<-time
out$status<-status
out$f<-f
out$method<-method
out$cmat<-cmat
out$k<-k

out$sample.estimate<-est

class(out)<-c("poly3ci", "sci")

return(out)

}

