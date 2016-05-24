"cdfkur" <-
function(x,para) {
   if(! are.parkur.valid(para)) return()
   A <- para$para[1]
   B <- para$para[2]

   f <- 1 - (1 - x^A)^B
   f[x < 0] <- 0
   f[x > 1] <- 1

   names(f) <- NULL
   return(f)
}
