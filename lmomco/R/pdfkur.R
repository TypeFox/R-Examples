"pdfkur" <-
function(x,para) {
   if(! are.parkur.valid(para)) return()
   A <- para$para[1]
   B <- para$para[2]

   lo <- quakur(0, para) # just in case a location/scale version
   hi <- quakur(1, para) # is ever implemented

   f <- A*B*x^(A-1)*(1-x^A)^(B-1)
   f[x < lo | x > hi] <- NA

   names(f) <- NULL
   f[! is.finite(f)] <- NA
   f[is.na(f)] <- 0 # decision Dec. 2015
   return(f)
}

