
mubz.normal<-function(b,z,
                      params=NULL,spatial=NULL,
                      cor.matrix=spatial.cor.matrix(c(params,1),spatial),
                      details=FALSE) {
    if(all(!b)) { # case B=emptyset
        r<-1
    } else {
        cov.matrix<-2*pi*cor.matrix
        sigma.b<-cov.matrix[b,b,drop=FALSE]
        sigma.b.inv<-solve(sigma.b)

        zb<-z[b]
        normzb<-as.numeric(sqrt(t(zb) %*% sigma.b.inv %*% zb))

        if(all(b)) { # case B=all
            r<-2^(length(b)/2)*gamma( (length(b)+1)/2 ) /
                sqrt( 2*(2*pi)^(length(b))*det(sigma.b) ) /
                    (normzb^(length(b)+1))
        } else {
            sigma.bc<-cov.matrix[!b,!b,drop=FALSE]
            sigma.bcb<-cov.matrix[!b,b,drop=FALSE]
            sigma.cond<-sigma.bc - sigma.bcb %*% sigma.b.inv %*% t(sigma.bcb)
            
            ztilde<-as.numeric(z[!b] - sigma.bcb %*% sigma.b.inv %*% zb)/normzb

            d2<-sum(!b)+1
            sigma.y<-matrix(nrow=d2, ncol=d2)
            sigma.y[1:(d2-1),1:(d2-1)]<-
                sigma.cond +
                    (ztilde %o% ztilde)
            sigma.y[1:(d2-1),d2] <- ztilde
            sigma.y[d2,1:(d2-1)] <- ztilde
            sigma.y[d2,d2]<-1

            r<-abs(pmnormpow(x=rep(0,d2),varcov=sigma.y,ipuiss=d2,puiss=sum(b))) /
                ( sqrt( (2*pi)^(sum(b)-1)*det(sigma.b) ) * normzb^(sum(b)+1) )
        }
    }
    if(details) {
        r
    } else {
        as.numeric(r)
    }
}
