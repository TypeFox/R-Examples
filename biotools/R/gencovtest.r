# S3 method for "gencovtest"
gencovtest <-
function(obj, 
	geneticFactor, 
	gcov = NULL,
	residualFactor = NULL, 
	adjNrep = 1, 
	test = c("MCPR", "Wilks", "Pillai"),
	nsim = 9999,
	alternative = c("two.sided", "less", "greater"))
		UseMethod("gencovtest")

# -----------------------------------
# manova method
gencovtest.manova <- 
function(obj, 
	geneticFactor, 
	gcov = NULL,
	residualFactor = NULL, 
	adjNrep = 1, 
	test = c("MCPR", "Wilks", "Pillai"),
	nsim = 9999,
	alternative = c("two.sided", "less", "greater"))
{
    manov <- anova(obj)
    model <- obj$model
    SS <- summary(obj)$SS
    #stopifnot(geneticFactor %in% names(model))
    dfg <- manov[geneticFactor, "Df"]
    Mg <- SS[[geneticFactor]] / dfg
    if (is.null(residualFactor)) {
        dfe <- df.residual(obj)
        Me <- SS[["Residuals"]] / dfe
    } else {
        #stopifnot(residualFactor %in% names(model))
        dfe <- manov[residualFactor, "Df"]
        Me <- SS[[residualFactor]] / dfe
    }
    test <- match.arg(test)
    alternative <- match.arg(alternative)
    nvar <- nrow(Mg)

    # genetic covariance matrix input/estimation
    if (is.null(gcov)) {
        nrep <- nrow(model) / nlevels(model[[geneticFactor]])
        gcov <- (Mg - Me) / (nrep * adjNrep)
    } else {
        if (nrow(gcov) != nvar || ncol(gcov) != nvar)
           stop("'gcov' presents incompatible dimensions")
    }
    gcor <- cov2cor(gcov)

    # collinearity diagnosis
    cn <- conditionNumber(gcov)
    if (cn > 100) { 
       mess <- paste("the genetic covariance matrix presents", 
          attr(cn, "meaning"))
       warning(mess)
    }

    # Aux. function for Wilks' Lambda and Pillai's Tn
    teststat <- function(m, var1, var2)
    {
       m2 <- m[c(var1, var2), c(var1, var2)]
       Lambda <- det(m2) / prod(diag(m2))
       Tn <- m2[1, 2]^2 / prod(diag(m2))
       out <- list(Lambda = Lambda, Tn = Tn)
       return(out)
    }

    # Wishart simulation
    ratio <- Mg / Me
    WG <- rWishart(nsim, dfg, Me) / dfg 
    WR <- rWishart(nsim, dfe, Me) / dfe
    simRatio <- windata(WG / WR, p = 0.01)
    dimnames(simRatio) = list(rownames(Mg), colnames(Mg), NULL)

    # test statistics, chi-sq approx (wilks and pillai) and p-values
    pval <- matrix(NA, nvar, nvar,
       dimnames = list(rownames(Mg), colnames(Mg)))
    stat <- X2 <- pval
    for(i in 1:nvar){
       for(j in 1:nvar) {
          if (i != j) {
             if (test == "MCPR") {
                 stat <- ratio
                 X2 = NULL
                 if (alternative == "two.sided") {
                     pval[i, j] <- mean(simRatio[i, j, ] - 
                         median(simRatio[i, j, ]) >= abs(stat[i, j]) -
                         median(simRatio[i, j, ])) +
                         mean(simRatio[i, j, ] - 
                         median(simRatio[i, j, ]) <= -abs(stat[i, j]) -
                         median(simRatio[i, j, ])) 
                     pval[i, i] <- mean(simRatio[i, i, ] >= stat[i, i])
                 } else {
                 stop("This type of alternative hypothesis is not implemented yet!")
                 }
             } else if (test == "Wilks") {
                 stat[i, j] <- teststat(gcov, i, j)$Lambda
                 X2[i, j] <- -dfg * log(stat[i, j])
                 pval[i, j] <- pchisq(X2[i, j], 1, lower.tail = FALSE)
             } else {
                 stat[i, j] <- teststat(gcov, i, j)$Tn
                 X2[i, j] <- dfg * stat[i, j]
                 pval[i, j] <- pchisq(X2[i, j], 1, lower.tail = FALSE)
             }
          }
       }
    }
    
    # output
    out <- list(gcov = gcov, gcor = gcor, 
       test = test, statistics = stat, 
       p.values = pval, alternative = alternative,
       X2 = X2, simRatio = simRatio,
       dfg = dfg, dfe = dfe)
    class(out) <- "gencovtest"
    return(out)
}

# ---------------------------
# print method
print.gencovtest <- function(x, digits = 4, ...)
{
    cat("\n          Genetic Covariance Test \n")
    cat("\nGenetic (Co)variances and Correlations (upper triangular):\n")
    GR <- lower.tri(x$gcov, diag = TRUE) * x$gcov + 
        upper.tri(x$gcor) * x$gcor 
    print(round(GR, digits))
    cat("\nTest:", x$test)
    if (x$test == "MCPR") {
        cat("\nAlternative hypothesis:", x$alternative, "\n")
        nsim <- dim(x$simRatio)[3L]
        cat("\nMean Sq and Cross-Prods Ratios and p-values (upper triangular) \nbased on", nsim, "estimates:\n")
        ratP <- lower.tri(x$statistics, diag = TRUE) * x$statistics + 
            upper.tri(x$p.values) * x$p.values
        print(round(ratP, digits))
    } else {
        cat("\n\nChi-Sq (df = 1) approx. and p-values (upper triangular):\n", sep = "")
        X2.p <- lower.tri(x$X2, diag = TRUE) * x$X2 + 
            upper.tri(x$p.values) * x$p.values 
        print(round(X2.p, digits))
    }
    invisible(x)
}

# -----------------------------
# plot method
plot.gencovtest <- function(x, var1, var2, ...)
{
    if (is.character(var1) & is.character(var2)) {
       pos <- match(c(var1, var2), colnames(x$statistics))
    } else {
       pos <- c(var1, var2)
    }
    if (x$test == "MCPR") {
        y <- x$simRatio[pos[1], pos[2], ]
        if (pos[1] != pos[2]) {
            plot(density(y), 
               main = paste("variables:", var1, "and", var2), ...)
            abline(v = median(y), col = "gray")
            lines(x = c(x$statistics[var1, var2], x$statistics[var1, var2]), 
               y = c(0, 0.3 * par("usr")[4]), col = 2)
            points(x = x$statistics[var1, var2], 
               y = 0.3 * par("usr")[4], col = 2, pch = 7)
        } else {
            plot(density(y), 
               main = paste("variable:", var1), ...)
            gl1 <- x$dfg
            gl2 <- x$dfe
            curve(df(x, gl1, gl2), add = TRUE, col = "blue", lty = 2)
            lines(x = c(x$statistics[var1, var2], x$statistics[var1, var2]), 
               y = c(0, 0.3 * par("usr")[4]), col = 2)
            points(x = x$statistics[var1, var2], 
               y = 0.3 * par("usr")[4], col = 2, pch = 7)
            legend("topright", c("Kernel density", paste("F(", gl1, ", ", 
               gl2, ")", sep = "")), lty = 1:2, col = c(1, 4), cex = 0.8)
        }
    } else {
        if (pos[1] != pos[2]) {
            curve(dchisq(x, df = 1), 
               main = paste("variables:", var1, "and", var2), ...)
            lines(x = c(x$X2[var1, var2], x$X2[var1, var2]), 
               y = c(0, 0.3 * par("usr")[4]), col = 2)
            points(x = x$X2[var1, var2], 
               y = 0.3 * par("usr")[4], col = 2, pch = 7)
        }
    }
}
