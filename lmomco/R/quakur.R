"quakur" <-
function(f,para,paracheck=TRUE) {
   if(! check.fs(f)) return()
   if(paracheck == TRUE) {
     if(! are.parkur.valid(para)) return()
   }
   A <- para$para[1]
   B <- para$para[2]

   x <- (1 - (1-f)^(1/B))^(1/A)
   names(x) <- NULL
   return(x)
}
