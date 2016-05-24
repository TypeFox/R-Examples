"pdfray" <-
function(x,para) {
   if(! are.parray.valid(para)) return()
   U <- para$para[1]
   A <- para$para[2]
   AA <- A^2

   Y <- x - U
   f <- Y/AA * exp(-(Y^2/(2*AA)))

   names(f) <- NULL
   f[Y <= 0] <- NA
   f[! is.finite(f)] <- NA
   f[is.na(f)] <- 0 # decision Dec. 2015
   return(f)
}

