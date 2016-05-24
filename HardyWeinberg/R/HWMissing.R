HWMissing <- function (X, imputecolumn = 1, m = 50, coding = c(0,1,2), verbose = FALSE, alpha = 0.05,
    varest = "oneovern", statistic = "chisquare", alternative = "two.sided", ...)
{
    skip <- FALSE #skip computations if no missing values are found.
    if (alpha <= 0 | alpha >= 1)
        stop("HWMissing: alpha should be in the range (0,1)")
    X <- as.data.frame(X)
    if(ncol(X)>1) {
      wholerowmissing <- apply(X,1,missingentirerow)
      if(sum(wholerowmissing) > 0) {
        warning(paste("HWMissing: there are",toString(sum(wholerowmissing)),
                      "rows with only missings, these are eliminated."))
        X <- X[!wholerowmissing,]
      }
      if(sum(is.na(X[,imputecolumn]))==0) {
        warning("HWMissing: there are no missing values for this marker")
        skip <- TRUE
      }
    }
    if(!skip) {
        if(ncol(X) == 1) { # with one column sample at random assuming MCAR
           X <- data.frame(MakeFactor(X[,1],coding), rep(1, length(X)))
           imp <- mice(X, m = m, predictorMatrix = quickpred(X), printFlag = FALSE, ...)
        } else {
           X[, imputecolumn] <- MakeFactor(X[, imputecolumn],coding)
           imp <- mice(X, m = m, predictorMatrix = quickpred(X), printFlag = FALSE, ...)
       }
    Xmat <- NULL
    for (i in 1:m) {
        Ximp <- complete(imp, i)
        nAA <- sum(Ximp[, imputecolumn] == "AA")
        nAB <- sum(Ximp[, imputecolumn] == "AB")
        nBB <- sum(Ximp[, imputecolumn] == "BB")
        Ximputedcounts <- c(AA = nAA, AB = nAB, BB = nBB)
        Xmat <- rbind(Xmat, Ximputedcounts)
    }
    rownames(Xmat) <- 1:m
    cout <- switch(statistic, chisquare = CombineChisquare(Xmat, alpha = alpha, varest = varest),
                              exact = CombineExact(Xmat),
                              stop("HWMissing: unknown value for argument statistic"))
    if(statistic=="chisquare") {
      Res <- c(cout$fhatimp, cout$llf, cout$ulf, cout$pvalimp,
          cout$r, cout$gamma, cout$pvalimpd, cout$pvalimpe)
      names(Res) <- c("f", "llci", "ulci", "p-value", "r", "gamma","p-v dearth","p-val excess")
      if (verbose) {
        cat("Test for Hardy-Weinberg equilibrium in the presence of missing values\n")
        cat("Inbreeding coefficient f = ", round(cout$fhatimp,
            digits = 4), "\n")
        cat(round(100 * (1 - alpha), digits = 0), "% Confidence interval (",
            round(cout$llf, digits = 4), ",", round(cout$ulf,
                digits = 4), ")\n")
        cat("p-value = ", round(cout$pvalimp, digits = 4), "\n")
        cat("Relative increase in variance of f due to missings: r = ",
            round(cout$r, digits = 4), "\n")
        cat("Fraction of missing information about f: lambda = ",
            round(cout$gamma, digits = 4), "\n")
      }
    }
     if(statistic=="exact") {
           pval <- switch(alternative,
                      greater = cout$mipgrea,
                      less = cout$mipless,
                      two.sided = 2*min(cout$mipgrea,cout$mipless),
                      stop("invalid value for parameter alternative"))
           if(verbose) {
               stringtwosided <- paste("Two-sided Exact test for Hardy-Weinberg equilibrium in the presence of missing values\n p-value = ",
                                       format(pval,scientific = FALSE), "\n")
               stringless <- paste("Exact test for heterozygote dearth in the presence of missing values\n p-value = ",
                                   format(pval,scientific = FALSE), "\n")
               stringgreater <- paste("Exact test for heterozygote excess in the presence of missing values\n p-value = ",
                                      format(pval,scientific = FALSE), "\n")
               toprint <- switch(alternative, two.sided = stringtwosided,
                                             greater = stringgreater,
                                             less = stringless)
              cat(toprint)
           }
         Res <- c(cout$mipless,cout$mipgrea)
    } 
    } # if(!skip)
    else {Res <- NA; Xmat <- NULL}
    return(list(Res = Res, Xmat = Xmat))
}

