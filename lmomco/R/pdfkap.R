"pdfkap" <-
function(x,para) {
   if(! are.parkap.valid(para)) return()
   # Function based on written communication of FORTRAN source
   # from J.R.M. Hosking in late October 2007.
   XI <- para$para[1]
   A  <- para$para[2]
   K  <- para$para[3]
   H  <- para$para[4]

   Fs <- cdfkap(x,para)
   Y  <- (x - XI)/A
   if(K != 0) {
      Y <- 1 - K*Y
      ops <- options(warn=-1)
      Y <- (1-1/K)*log(Y)
      options(ops)
   }
   Y <- exp(-Y)
   f <- Y/A * Fs^(1-H)

   names(f) <- NULL
   f[! is.finite(f)] <- NA
   f[is.na(f)] <- 0 # decision Dec. 2015
   return(f)
}

