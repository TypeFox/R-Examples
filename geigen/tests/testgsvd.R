
testgsvd <- function(z,A,B) {
    tol <- 100*.Machine$double.eps

    if( inherits(z,"xzgsvd") ) {
        chkA <- Conj(t(z$U)) %*% A %*% z$Q
        chkB <- Conj(t(z$V)) %*% B %*% z$Q
    } else if( inherits(z,"xdgsvd") ) {
        chkA <- t(z$U) %*% A %*% z$Q
        chkB <- t(z$V) %*% B %*% z$Q
    } else stop("argument z not result of gsvd")

    oR <- gsvd.oR(z)
    D1 <- gsvd.D1(z)
    D2 <- gsvd.D2(z)

#    print(chkA-D1 %*% oR)
#    print(chkB-D2 %*% oR)

    ret <- logical(4)
#    print(abs(chkA-D1 %*% oR))
#    print(abs(chkB-D2 %*% oR))
    ret[1] <- all(abs(chkA-D1 %*% oR) <= tol)
    ret[2] <- all(abs(chkB-D2 %*% oR) <= tol)
    ret[3] <- all(abs(t(D1) %*% D1 + t(D2) %*% D2 - diag(1,nrow=z$k+z$l)) <= tol)

#   check alpha^2 + beta^2 == 1
    salpha <- 0
    sbeta  <- 0
    m <- z$m
    k <- z$k
    l <- z$l
    if(m-k-l >= 0 ) {
        if(l > 0) {
            salpha <- z$alpha[(k+1):(k+l)]^2
            sbeta  <- z$beta[(k+1):(k+l)]^2
        }
    } else {
        if(m > k) {
            salpha <- z$alpha[(k+1):m]^2
            sbeta  <- z$beta[(k+1):m]^2
        }
    }
#    print(salpha)
#    print(sbeta)
#    print(salpha+sbeta)
    ret[4] <- all(salpha+sbeta-1<=tol)

    ret
}
