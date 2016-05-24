"aprox.Qrbar" <-
function(x)
 {
 vr <- varr(x)
 N <- sum(x$n,na.rm=TRUE)
 rb <- rbar(x)
 chi <- (N/(1-rb^2)^2)*vr
 k <- length (x$Rxy[!(is.na(x$Rxy))])
 pval <- 1 - pchisq(chi, k-1)
 mat <- matrix(c(chi,k-1,pval),ncol=3)
 colnames(mat) <- c("CHISQ", "df", "p-val")
 return(mat)
 }

