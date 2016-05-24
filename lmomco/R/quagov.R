"quagov" <-
function(f, para, paracheck=TRUE) {
   if(! check.fs(f)) return()
   if(paracheck == TRUE) {
      if(! are.pargov.valid(para)) return()
   }
   U <- para$para[1]
   A <- para$para[2]
   B <- para$para[3]
   Bp1 <- B + 1

   x <- U + A*(Bp1*f^B - B*f^Bp1)

   names(x) <- NULL
   return(x)
}

