summary.joint <- function (object, variance = TRUE, ...) {
cat("Random effects joint model\n")
model <- object$model
cat(" Data:", deparse( object$call$data ), "\n")
cat(" Log-likelihood:", object$loglik$jointlhood, "\n")
cat("\n")
cat("Longitudinal sub-model fixed effects:",deparse(object$formulae$lformula))
lfixed <- object$coefficients$fixed$longitudinal
names(lfixed) <- ""
print(lfixed)
cat("\n")
cat("Survival sub-model fixed effects:",deparse(object$formulae$sformula))
sfixed <- data.frame(object$coefficients$fixed$survival)
if (sum(dim(sfixed)) == 0) {
cat("\n", "No survival baseline covariates specified", "\n")
}
else {
names(sfixed) <- ""
print(sfixed)
}
cat("\n")
cat("Latent association:")
lat <- data.frame(object$coefficients$latent)
names(lat) <- ""
print(lat)
cat("\n")
cat("Variance components:\n")
sigu <- diag(object$sigma.u)
names(sigu) <- rownames(object$sigma.u)
sigz <- object$sigma.z
names(sigz) <- "Residual"
vars <- c(sigu, sigz)
names(vars) <- c(names(sigu), names(sigz))
if(!variance) {vars = sqrt(vars)}
print(vars)
if(!variance) {cat(" Note: the above are standard deviations\n")}
cat("\n")
if(object$convergence == TRUE) {
    cat("Convergence at iteration:", object$numIter,"\n")
}
else {
    cat("Convergence not achieved\n")
}
cat("\n")
cat("Number of observations:", dim(object$data$longitudinal)[1],"\n")
cat("Number of groups:", dim(object$data$survival)[1],"\n")
object$nobs <- dim(object$data$longitudinal)[1]
object$ngrps <- dim(object$data$survival)[1]
class(object) <- c("summary.joint", class(object))
invisible(object)
}
