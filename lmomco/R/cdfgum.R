"cdfgum" <-
function(x,para) {
   if(! are.pargum.valid(para)) return()
   U <- para$para[1]
   A <- para$para[2]

   Y <- -(x-U)/A
   f <- exp(-exp(Y)) # Condition GEV(K=0), HW1997, p.195
   names(f) <- NULL
   return(f)
}

