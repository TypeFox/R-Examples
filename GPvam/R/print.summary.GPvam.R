print.summary.GPvam <-
function (x, ...) 
{
    cat("Number of observations:", x$num.obs, "\n")
    cat("Number of years:", x$num.year, "\n")
    cat("Number of students:", x$num.student, "\n")
    for (i in 1:x$num.year) {
        cat("Number of teachers in year", x$key[i,1], ":", x$num.teach[i], 
            "\n")
    }
    cat("\n")
    cat("Number of EM iterations: ", x$iter, "\n")
    cat("-2 log-likelihood", -2 * x$loglik, "\n")
    cat("AIC", x$AIC, "\n")
    cat("AICc", x$AICc, "\n")
    cat("\n")
    for (i in 1:x$num.year) {
        cat("Covariance matrix for current and future year\neffects of year", 
            x$key[i,1], "teachers.\n")
        print(as.matrix(x$teach.cov[[i]]))
        cat("\n")
        cat("with correlation matrix\n")
        print(round(cov2cor(as.matrix(x$teach.cov[[i]])), 4))
        cat("\n")
    }
    if (!any(is.na(x$R_i))) {
        cat("Block of error covariance matrix (R):\n")
        print(as.matrix(x$R_i))
        cat("with correlation matrix\n")
        print(round(cov2cor(as.matrix(x$R_i)), 4))
        cat("\n")
    }
    if (!is.na(x$stu.cov)) {
        cat("Student variance component:\n")
        print(x$stu.cov)
        cat("\n")
    }
    cat("Parameter estimates:\n")
    print(x$parameters)
    cat("\n")
    cat("Distribution of marginal residuals\n")
    print(summary(x$mresid))
    cat("\n")
    cat("Distribution of raw conditional residuals\n")
    print(summary(x$cresid))
    cat("\n")
    cat("Distribution of scaled conditional residuals\n")
    print(summary(x$sresid))
    cat("\n")

   
}
