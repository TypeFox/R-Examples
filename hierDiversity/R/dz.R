dz <-
function(abundances, lev = "beta", q = 1) {
   wts <- FALSE
   if (wts[1] == FALSE) {
      if (is.matrix(abundances)) {
         wts <- rep(1, dim(abundances)[1])
      }
      if (is.vector(abundances)) {
         wts <- 1
      }
   }
   w <- normalizeRows(wts)
   p <- normalizeRows(abundances)
   if (lev == "alpha") {
      if (q != 1) {
         if (is.matrix(abundances)) {
            numerator <- sum(as.matrix(w^q) * 
              as.matrix(apply(p, 1, pqSum, 
               q = q)))
         }
         if (is.vector(abundances)) {
            numerator <- pqSum(p, q = q)
         }
         D.VALUE <- ((numerator/(sum(w^q)))^(1/(1 - 
            q)))
      }
      if (q == 1) {
         if (is.matrix(abundances)) {
            numerator <- sum(as.matrix(-1 * 
              w) * as.matrix(apply(p, 1, pqSum, 
              q = q)))
         }
         if (is.vector(abundances)) {
            numerator <- -1 * pqSum(p, q = q)
         }
         D.VALUE <- exp(numerator)
      }
   }
   if (lev == "beta") {
      D.VALUE <- dz(abundances, lev = "gamma", q = q)/dz(abundances, lev = "alpha", q = q)
   }
   if (lev == "gamma") {
      D.VALUE <- dz(t(as.matrix(colSums(p * 
         w))), lev = "alpha", q = q)
   }
   return(round(D.VALUE, digits = 3))
}
