checkIMRcollinearity <- function(X, tol=1e6) {
   ## This small utility checks whether inverse Mills ratio is (virtually) collinear to the other explanatory
   ## variables.  IMR is in the last columns.
   ## In case of collinearity it returns TRUE, otherwise FALSE
   X <- X[!apply(X, 1, function(row) any(is.na(row))),]
   if(kappa(X) < tol)
      return(FALSE)
   if(kappa(X[,-ncol(X)]) > tol)
      return(FALSE)
                           # it is multicollinear,
                           # but not (just) due to IMR
   return(TRUE)
}
