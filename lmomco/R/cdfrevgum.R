"cdfrevgum" <-
function(x,para) {
   if(! are.parrevgum.valid(para)) return()
   U <- para$para[1]
   A <- para$para[2]

   Y <- (x-U)/A
   f <- (1 - exp(-exp(Y)))
   names(f) <- NULL
   return(f)
}
