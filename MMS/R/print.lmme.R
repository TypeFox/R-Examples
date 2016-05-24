print.lmme <-
function(x,...) {
    # 0. part
    print(x$call)
cat("\n")
    # 1. part
cat("Random effects:","\n")
cat("Variance of the random effects:\n")
print(x$Psi)
cat("\n")
cat("Variance of the residuals", x$sigma_e,"\n\n")

    # 2. part
cat("Fixed effects:", "\n")
cat(x$beta,"\n")

}
