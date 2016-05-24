`Bonferroni.sig` <-
function (x, model = "codominant", alpha = 0.05, include.all.SNPs = FALSE) 
{
    if (!inherits(x, "WGassociation")) 
        stop("x must be a 'WGassociation' object")
    x <- attr(x, "pvalues")
    model.type <- names(x)
    m <- charmatch(model, model.type, nomatch = 0)
    if (m == 0) 
        stop("this model is was not fitted")
    temp1 <- grep("no", as.character(x[, m ]))
    temp2 <- c(1:nrow(x))[is.na(x[, m ])]
    temp <- c(temp1, temp2)
    if (!(include.all.SNPs) & length(temp)>=1) {
        x <- x[-temp, c(1, (m ))]
        cut.p <- alpha/nrow(x)
    }
    else cut.p <- alpha/nrow(x)
    cat("number of tests: ", nrow(x), "\n")
    cat("alpha:", alpha, "\n")
    cat("corrected alpha:", cut.p, "\n")
    significant <- x[as.numeric(x[, 2]) <= cut.p, ]
    if (all(is.na(significant))) {
        cat("   No significant SNPs after Bonferroni correction \n")
        ans <- NULL
    }
    else {
        ans <- significant
        print(as.matrix(ans), na.print = "-", quote = FALSE)
    }
    invisible(ans)
}

