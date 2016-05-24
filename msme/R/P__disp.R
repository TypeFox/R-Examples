P__disp <- function(x) {
   pr <- sum(residuals(x, type="pearson")^2)
   dispersion <- pr / x$df.residual
   return(c(pearson.chi2 = pr, dispersion = dispersion))
}
