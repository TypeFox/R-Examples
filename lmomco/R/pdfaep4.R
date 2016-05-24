"pdfaep4" <-
function(x, para, paracheck=TRUE) {
   if(paracheck == TRUE) {
     if(! are.paraep4.valid(para)) return()
   }

   U <- para$para[1]
   A <- para$para[2]
   K <- para$para[3]
   H <- para$para[4]

   # The following appears unnecessary as the AEP4 seemingly sweeps through
   # either the Normal or Laplace distributions as the K or H pass near the
   # critical points. Hence, all of this code is commented out.
   #SMALL <- 1E-6
   #if(abs(K - 1) < SMALL & abs(H - 2) < SMALL) {
   #   warning("Normal distribution being used for k=", K, " and h=", H)
   #   SIGMA <- 0.39894228 * sqrt(pi) * A
   #   MU    <- U
   #   return(pdfnor(x, vec2par(c(MU, SIGMA), type="nor")))
   #}
   #if(abs(K - 1) < SMALL & abs(H - 1) < SMALL) {
   #    warning("Laplace distribution being used for k=", K, " and h=", H)
   #    A  <- A # L2lap = 0.75 * Alap;  L2aep = 0.75 * Aaep
   #    XI <- U
   #    return(pdflap(x, vec2par(c(XI, A), type="lap")))
   #}

   Z <- H*K / ( A * (1 + K*K) * gamma(1/H) )

   Y <- abs(x - U) / A
   f <- Z * exp(-1 * (  K^sign(x - U) * Y )^H )

   names(f) <- NULL
   f[! is.finite(f)] <- NA
   f[is.na(f)] <- 0 # decision Dec. 2015
   return(f)
}
