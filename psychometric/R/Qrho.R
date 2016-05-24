"Qrho" <-
function(x, aproxe=FALSE)
 {
 if(!aproxe) { ve <- vare(x)}
 else {ve <- aprox.vare(x)}
 k <- length (x$Rxy[!(is.na(x$Rxy))])
 vr <- varr(x)
 vav <- varAV(x)
q <- (k*vr)/(vav+ve)
pval <- 1 - pchisq(q, k-1)
 mat <- matrix(c(q,k-1,pval),ncol=3)
 colnames(mat) <- c("CHISQ", "df", "p-val")
 return(mat)
 }

