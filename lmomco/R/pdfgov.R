"pdfgov" <-
function(x,para) {
   if(! are.pargov.valid(para)) return()
   U <- para$para[1]
   A <- para$para[2]
   B <- para$para[3]
   ARG <- 1 / (A*B*(B+1))

   lo <- U
   hi <- U+A

   F <- cdfgov(x, para)
   F[x < lo] <- NA
   F[x > hi] <- NA
   f <- ARG * F^(1-B) / (1 - F)

   names(f) <- NULL
   f[! is.finite(f)] <- NA
   f[is.na(f)] <- 0 # decision Dec. 2015
   return(f)
}


