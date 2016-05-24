"print.dirichlet" <-
function(x,...){
  obj <- x
  if (obj$error==1) {
    cat("ERROR! nstar is too small! (nstar=",obj$nstar,"), Sum of Pn is (should be 1)",
        sum(sapply(0:obj$nstar,obj$Pn)), "\n")
  }
  cat("Number of Brands in the Category =", obj$nbrand, "\n")
  cat("Brand List",obj$brand.name, sep=" : ")
  cat("\nBrands' Market Shares:", round(obj$brand.share,3),"\n")
  cat("Brands' Penetration:  ", round(obj$brand.pen.obs,3),"\n")
  obj$period.print()
#   cat("\nCategory Penetration =", round(obj$cat.pen,2),
#       ", with Buying Rate =", round(obj$cat.buyrate,2),"\n")
  cat("\nCategory Penetration =", round(obj$cat.pen,2),
      ", with Buying Rate =", round(obj$cat.buyrate,2),"\n")
  cat("Estimated Dirichlet Model Parameters:\n")
  cat("NBD: M =", round(obj$M,2), ",  K =", round(obj$K,2),";  ")
  cat("Dirichlet: S =", round(obj$S,2), "\n\n")
}

