"bRyy" <-
function(x)
 {
 Ryy <- x$Ryy
 n <- length (x$Ryy[!(is.na(x$Ryy))])
 b <- mean(sqrt(Ryy),na.rm=TRUE)
 vb <- var(sqrt(Ryy),na.rm=TRUE)*(n-1)/n
 out <- list(b,vb)
 return(out)
 }

