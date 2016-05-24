kfolds2coeff <- function(pls_kfolds) {
if (is.null(pls_kfolds$coeffs_kfolds)) {cat("No coefficients found\n"); return(NULL)}
if (length(pls_kfolds$coeffs_kfolds) != 1) {cat("Only if NK=1 and Jackknife (Loo) computations\n"); return(NULL)}
coef.all <- matrix(unlist(pls_kfolds$coeffs_kfolds[[1]]),nrow=length(pls_kfolds$coeffs_kfolds[[1]]),byrow=T)
if (is.null(pls_kfolds$folds)) {attr(coef.all,"folds") <- pls_kfolds$folds}
return(coef.all)}
