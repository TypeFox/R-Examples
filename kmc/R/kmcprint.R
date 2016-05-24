###########################################################################
#      KM EL with constraint, S3 methods!!!
###########################################################################

print.kmcS3 <- function (x,digits = max(3, getOption("digits") - 3), type = "plain",...)
{
    if (type=='md') {
        cat("|*A Recursive Formula for the Kaplan-Meier Estimator with Constraint*|\n")
        cat(":--------------------------------------------------------------------:\n")
        cat("|Information:|\n")
        cat("|Number of Constraints:\t*", length(x$g), "*|\n\n")
        cat("\n-------------------------------------------------------------------------------\n")
        cat("|Log-likelihood(Ha) | Log-likelihood(H0) | -2LLR | p-Value(df=",
        length(x$g), ")|\n")
        cat("|:---------------:|:--------------------:|:-------------------:|:--------------------:|:--------------------:|\n")
        cat("|", x[[1]], " | ", x[[2]], " | ", x[[3]], " | ",
        1 - pchisq(x[[3]], df = length(x$g)),'|',x$lambda,"|\n")
    }else {
        cat("\n---------------------------------------------------------------------------------\n")
        cat("A Recursive Formula for the Kaplan-Meier Estimator with Constraint\n")
        cat("Information:\n")
        cat("Number of Constraints:\t", length(x$g), "\n")
        cat("lamda(s):\t",x$lambda,'\n');
        cat("\n---------------------------------------------------------------------------------\n")
        names <- c("Log-likelihood(Ha)", "Log-likelihood(H0)",
        "-2LLR", paste("p-Value(df=", length(x$g), ")",sep = ""))
        re <- matrix(c(x[[1]], x[[2]], x[[3]], 1 - pchisq(x[[3]],
        length(x$g))), nrow = 1)
        colnames(re) <- names
        rownames(re) <- "Est"
        print.default(format(re, digits = digits), print.gap = 2,
        quote = FALSE, df = length(x$g))
        cat("---------------------------------------------------------------------------------\n")
        
    }
}


#summary.kmcS()
