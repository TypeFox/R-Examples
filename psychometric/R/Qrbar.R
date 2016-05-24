"Qrbar" <-
function(x)
 {
 r <- x$Rxy
 n <- x$n
 rb <- rbar(x)
 chi <- sum((((n-1)*(r-rb)^2)/(1-rb^2)^2),na.rm=TRUE)
 k <- length (x$Rxy[!(is.na(x$Rxy))])
 pval <- 1 - pchisq(chi, k-1)
 mat <- matrix(c(chi,k-1,pval),ncol=3)
 colnames(mat) <- c("CHISQ", "df", "p-val")
 return(mat)
}

