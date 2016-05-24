zScoreFDR <- function (esets, classes, useREM = TRUE, nperm = 1000, CombineExp = 1:length(esets)) 
{
    for (i in 1:length(classes)) {
        if (!is.factor(classes[[i]])) {
            classes[[i]] <- factor(classes[[i]])
        }
        if (nlevels(classes[[i]]) != 2) {
            stop("Error: Each list in the argument \"classes\" must contain only 2 levels.")
        }
        else {
            Ref <- levels(classes[[i]])[1]
            classes[[i]] <- sapply(classes[[i]], function(x) ifelse(x == 
                Ref, 0, 1))
        }
    }
    num.studies <- length(esets)
    num.genes <- nrow(esets[[1]])
    zscoresAll <- zScores(esets, classes, useREM = useREM, CombineExp = CombineExp)
    MuQu <- zscoresAll[, c("MUvals", "MUsds", "Qvals", "df", 
        "Qpvalues", "Chisq")]
    zscore <- zscoresAll[, c(paste("zSco_Ex_", 1:num.studies, 
        sep = ""), "zSco")]
    aperms <- replicate(nperm, zScores(esets, lapply(classes, 
        sample), useREM, CombineExp = CombineExp)[, c(paste("zSco_Ex_", 
        1:num.studies, sep = ""), "zSco")], simplify = FALSE)
    k <- c("pos", "neg", "two.sided")
    All <- vector(mode = "list", length = 3)
    for (l in 1:length(k)) {
        theFDR <- multExpFDR(zscore, aperms, type = k[l])
        n <- num.studies + 1
        i <- 1:n
        theResult <- matrix(NA, ncol = n * 2, nrow = nrow(zscore))
        rownames(theResult) <- rownames(zscore)
        stopifnot(rownames(theResult) == rownames(theFDR))
        theResult[, 2 * i - 1] <- zscore
        theResult[, 2 * i] <- theFDR
        colnames(theResult) <- 1:(2 * n)
        colnames(theResult)[2 * i - 1] <- paste("zSco_Ex_", i, 
            sep = "")
        colnames(theResult)[2 * i] <- paste("FDR_Ex_", i, sep = "")
        colnames(theResult)[(2 * n - 1):(2 * n)] <- c("zSco", 
            "FDR")
        All[[l]] <- cbind(theResult, MuQu)
    }
    names(All) <- k
    return(All)
}