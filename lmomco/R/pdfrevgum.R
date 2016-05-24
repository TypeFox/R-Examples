"pdfrevgum" <-
function(x,para) {
   if(! are.parrevgum.valid(para)) return()
   U <- para$para[1]
   A <- para$para[2]

   Y <- (x-U)/A
   f <- (1/A) * exp(Y) * ( exp( -exp(Y) ) )
   names(f) <- NULL
   f[! is.finite(f)] <- NA
   f[is.na(f)] <- 0 # decision Dec. 2015
   return(f)
}

