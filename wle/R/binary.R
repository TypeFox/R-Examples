binary <- function(x, dim) {

   if (x==0) {
       pos <- 1
   } else {
       pos <- floor(log(x, 2))+1
   }

   if (!missing(dim)) {
       if (pos<=dim) {
           pos <- dim
       } else {
           warning("the value of `dim` is too small")
       }  
   }

   bin <- rep(0, pos)
   dicotomy <- rep(FALSE, pos)
   for (i in pos:1) {
        bin[i] <- floor(x/2^(i-1))
        dicotomy[i] <- bin[i]==1
        x <- x-((2^(i-1))*bin[i])
   }
   return(list(binary=bin, dicotomy=dicotomy))
}


