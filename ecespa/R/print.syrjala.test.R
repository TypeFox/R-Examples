print.syrjala.test <-
function (x, ...) 
{
   cat("Cramer-von Misses test for the difference between\nthe spatial distributions of ",x$datanames[1], "and", x$datanames[2], "\nbased on",x$nperm, "permutations.\n\n")
   cat("   psi:    ", x$cvm.obs, "\n")
   cat("   p-value:", sum(x$cvm.sim >= x$cvm.obs)/(x$nperm + 
   1), "\n\n")
   cat("Kolmogorov-Smirnov test for the difference between\nthe spatial distributions of ",x$datanames[1], "and", x$datanames[2], "\nbased on",x$nperm, "permutations.\n\n")
   cat("   psi:    ", x$ks.obs, "\n")
   cat("   p-value:", sum(x$ks.sim >= x$ks.obs)/(x$nperm + 
   1), "\n\n")
}

