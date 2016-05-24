`as.snp` <-
function (x, ...) 
 {
   if (is.snp(x)) x else snp(x, ...)
 }

