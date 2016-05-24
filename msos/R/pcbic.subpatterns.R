pcbic.subpatterns <-
function(eigenvals,n,pattern0) {
     b <- NULL
     pts <- NULL
     k <- length(pattern0)
     if(k==1) return(F)
     for(i in 1:(k-1)) {
          p1 <- pcbic.unite(pattern0,i)
          b2 <- pcbic(eigenvals,n,p1)
          b <- c(b,b2$BIC)
          pts <- cbind(pts,p1)
     }
     list(bic=b,pattern=pts)
}
