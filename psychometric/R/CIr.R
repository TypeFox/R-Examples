"CIr" <-
function (r, n, level=.95)
 {
 
 z <- r2z(r)
 uciz <- CIz(z, n, level)[2]
 lciz <- CIz(z, n, level)[1]
 ur <- z2r(uciz)
 lr <- z2r(lciz)
 mat <- list(lr,ur)
 return(as.numeric(mat))
 }

