hctest <-
function(x, res)
{
   n <- length(x)
   stopifnot(n==length(res) && n>2)
   indjx <- order(x)
   resx <- res[indjx]
   cor.test(resx[-1], resx[-n], method="kendall")$p.value  
 }

