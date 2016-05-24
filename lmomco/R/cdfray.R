"cdfray" <-
function(x,para) {
   if(! are.parray.valid(para)) return()
   U <- para$para[1]
   A <- para$para[2]

   Y <- x - U
   f <- (1 - exp(-1*(Y)^2/(2*A^2)))

   names(f) <- NULL
   f[Y <= 0] <- 0
   return(f)
}

