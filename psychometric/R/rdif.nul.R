"rdif.nul" <-
function (r1, r2, n1, n2)
 {
 z1 <- r2z(r1)
 z2 <- r2z(r2)
 z <- (z1 - z2)/sqrt(1/(n1-3)+1/(n2-3))
 p <- pnorm(z)
 return(data.frame(zDIF = z, p = 1-p))
 }

