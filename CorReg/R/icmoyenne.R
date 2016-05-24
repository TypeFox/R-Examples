icmoyenne <- function(emoy=emoy, evar=evar, n=n, alpha=0.05) {
   inter = alpha/2#seuil de confiance
   lower = emoy + qt(inter, n - 1) * sqrt(evar/n)
   upper = emoy + qt(1 - inter, n - 1) * sqrt(evar/n)
   return(c(lower, upper))
}

