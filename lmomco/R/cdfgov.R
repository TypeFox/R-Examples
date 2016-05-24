"cdfgov" <-
function(x, para) {
   if(! are.pargov.valid(para)) return()
   U <- para$para[1]
   A <- para$para[2]
   B <- para$para[3]

   lo <- U
   hi <- U+A

   "afunc" <- function(x, f) return(x - quagov(f, para))
   f <- sapply(1:length(x), function(i) {
                   tmp <- NULL
               try(tmp <- uniroot(afunc, lower=0, upper=1, x=x[i]), silent=TRUE)
               ifelse(is.null(tmp), return(NA), return(tmp$root)) } )
   f[x < lo] <- NA
   f[x > hi] <- NA

   names(f) <- NULL
   f[! is.finite(f)] <- NA
   return(f)
}

