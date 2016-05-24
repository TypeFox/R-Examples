sword <- function(h, k, GSD, x=GSD$b[k]) {
    if(k == 1) 1 - pnorm(x-h*sqrt(GSD$t[k]*GSD$Imax))
    else {
        seqmon(a=GSD$a[1:k],
               b=c(GSD$b[1:(k-1)],x)-h*sqrt(GSD$t[1:k]*GSD$Imax),
               t=(GSD$t[1:k]*GSD$Imax)/GSD$t[1]*GSD$Imax,int=500*array(c(1),k)
               )[2*k]
    }
}

