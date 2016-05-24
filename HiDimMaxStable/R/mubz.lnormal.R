
mubz.lnormal<-function(b,z,
                       params=NULL,spatial=NULL,
                       cov.matrix=spatial.cor.matrix(params,spatial),
                       details=FALSE) {
    if(all(!b)) { # case B=emptyset
        r<-1
    } else {
        sigma.b<-cov.matrix[b,b,drop=FALSE]
        sigma.b.inv<-solve(sigma.b)
 
        sigma.bc<-cov.matrix[!b,!b,drop=FALSE]
        sigma.bcb<-cov.matrix[!b,b,drop=FALSE]
        sigma.cond<-sigma.bc - sigma.bcb %*% sigma.b.inv %*% t(sigma.bcb)
 
        zb<-z[b]
 
        e<-rep(1,length(z))
        eb<-e[b]
        ebc<-e[!b]
        nu<-diag(cov.matrix)/2
 
        ebc.proj<- ebc - sigma.bcb %*% sigma.b.inv %*% eb
 
        normeb2<-as.numeric(t(eb) %*% sigma.b.inv %*% eb)
        normznub2<-as.numeric(t(log(zb)+nu[b]) %*% sigma.b.inv %*% (log(zb)+nu[b]))
        prod.ebzbnub<-as.numeric(t(eb) %*% sigma.b.inv %*% (log(zb)+nu[b]))
 
        lztb<-log(z[!b])+nu[!b]-sigma.bcb %*% sigma.b.inv %*% (log(zb)+nu[b]) -
            ebc.proj*(prod.ebzbnub-1)/normeb2

        sigma.yi<-
            sigma.cond +
                (as.numeric(ebc.proj) %o% as.numeric(ebc.proj))/normeb2
 
        1/( ((2*pi)^((sum(b)-1)/2)) * sqrt( normeb2 * det(sigma.b) )
           * prod(zb) ) *
               exp((prod.ebzbnub-1)^2/(2*normeb2) - normznub2/2) *
                   ifelse(all(b),1,pmnorm(lztb,mean=rep(0,sum(!b)),varcov=sigma.yi))
    }
}
