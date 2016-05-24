"quagum" <-
function(f,para,paracheck=TRUE) {
   if(! check.fs(f)) return()
   if(paracheck == TRUE) {
     if(! are.pargum.valid(para)) return()
   }
   U <- para$para[1]
   A <- para$para[2]

   x <- U-A*log(-log(f))
   names(x) <- NULL
   return(x)
}

