`BackcastResidualsAR` <-
function(y,phi,Q=100,demean=TRUE){
if (demean)
    z <- y-mean(y)
else
    z <- y
p <- length(phi)
if (p==0)
    a<-z
else {
    n<-length(z)
    p<-length(phi)
    nQ<-n+Q
    a<-e<-zF<-zR<-numeric(nQ)
    r<-p+1
    zR[1:n]<-rev(z)
    zF[Q+(1:n)]<-z
    for (i in r:n)
        e[i]<-zR[i]-crossprod(phi,zR[i-(1:p)])
    for (i in 1:Q)
        zR[n+i]<-crossprod(phi,zR[(n+i)-(1:p)])
    zF[1:Q]<-rev(zR[n+(1:Q)])
    for (i in r:nQ)
        a[i]<-zF[i]-crossprod(phi,zF[i-(1:p)])
    zF[Q+(1:n)]<-z
    a<-a[Q+(1:n)]
    }
a
}

