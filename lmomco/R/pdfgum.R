"pdfgum" <-
function(x,para) {
   if(! are.pargum.valid(para)) return()
   U <- para$para[1]
   A <- para$para[2]

   Y <- -(x - U)/A
   f <- A^(-1)*exp(Y)*exp(-exp(Y))

   names(f) <- NULL
   f[! is.finite(f)] <- NA
   f[is.na(f)] <- 0 # decision Dec. 2015
   return(f)
}

