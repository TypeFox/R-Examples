mantelTest <- 
function(m1, m2, nperm = 999,
	alternative = "greater", graph = TRUE, 
	main = "Mantel's test",
	xlab = "Correlation", ...)
{
    if(!inherits(m1, c("matrix", "dist")))
	  stop("'m1' must be an object of class 'matrix' or 'dist'")
    if(!inherits(m2, c("matrix", "dist")))
	  stop("'m2' must be an object of class 'matrix' or 'dist'")
    m1 <- as.matrix(m1)
    m2 <- as.matrix(m2)
    if(diff(dim(m1)) != 0 || diff(dim(m2)) != 0)
	  stop("'m1' and 'm2' must be square matrices!")
    n <- nrow(m1)
    if (n != nrow(m2)) 
        stop("Incompatible dimensions!")
    alternative <- match.arg(alternative,
       c("two.sided", "less", "greater"))

    # Observed correlation
    obs <- cor(as.dist(m1), as.dist(m2))

    # Matrix permutation
    perm <- function(m)
    {
       o <- sample(nrow(m))
       return(m[o, o])
    }
    nullcor <- NULL
    for(i in 1:nperm) {
       nullcor[i] <- cor( as.dist(m1),
	    as.dist(perm(m2)) )
    }

    # Plot density
    if(graph) {
	  plot(density(nullcor, na.rm = TRUE),
           main = main, xlab = xlab, ...)
        points(x = c(obs, obs), 
           y = c(0, 0.3 * par("usr")[4]), 
           type = "l", col = 2)
        points(x = obs, y = 0.3 * par("usr")[4], 
           col = 2, pch = 7)
    }

    # p-value
    pval <- simpval(nullcor, obs, alternative)

    # Percentile CI
    #qu <- quantile(nullcor, p = c(alpha/2, 1 - alpha/2))
    #CI <- c(qu[1], qu[2])
    #attributes(CI) <- list(conf.level = 1 - alpha)

    # Output
    out <- list(correlation = obs,
       #CI = CI,
       p.value = pval,
       alternative = alternative,
       nullcor = nullcor)
    class(out) <- "mantelTest"
    return(out)
}

# -------------------------------------
# print method
print.mantelTest <- function(x, digits = 4, ...)
{
    cat("\n            Mantel's permutation test",
	 "\n\nCorrelation: ", x$correlation,
       #"\nPercentile CI (", 100 * attr(x$CI, "conf.level"), "%): ", 
       #   "[", x$CI[1], ", ", x$CI[2], "]",
	 "\np-value: ", x$p.value,
       ", based on ", length(x$nullcor), 
          " matrix permutations", sep = "")
    if (x$alternative == "greater") {
       mes <- "greater than 0"
    } else if (x$alternative == "two.sided") {
       mes <- "not equal to 0"
    } else { 
       mes <- "less than 0"
    }
    cat("\nAlternative hypothesis: true correlation is", mes, "\n")
    invisible(x)
}
