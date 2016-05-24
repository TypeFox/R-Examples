zScores <- function (esets, classes, useREM = TRUE, CombineExp = 1:length(esets)) 
{

   num.studies <- length(esets)
    num.genes <- nrow(esets[[1]])
    if (num.studies != length(classes)) 
        stop("Length of classes must be equal to length of esets.")
    for (i in 1:num.studies) {
        if (!is.factor(classes[[i]])) {
            classes[[i]] <- factor(classes[[i]])
        }
        if (nlevels(classes[[i]]) != 2) {
            stop("Error: Each list in the argument \"classes\" must contain only 2 levels.")
        }        else {
            Ref <- levels(classes[[i]])[1]
            classes[[i]] <- sapply(classes[[i]], function(x) ifelse(x == 
                Ref, 0, 1))
        }
    }
    tau2 <- function(Q, num.studies, my.weights) {
        vwts <- rowSums(my.weights)
        tmp2 <- rowSums(my.weights^2)
        tau2 <- pmax(0, (Q - (num.studies - 1))/(vwts - tmp2/vwts))
        return(tau2)
    }
    theNames <- rownames(esets[[1]])
    for (i in 2:num.studies) stopifnot(identical(theNames, rownames(esets[[i]])))
    ds <- matrix(NA, ncol = num.studies, nrow = num.genes)
    vars <- matrix(NA, ncol = num.studies, nrow = num.genes)
    for (i in 1:length(esets)) {
        my.d.adj <- dstar(getdF(esets[[i]], classes[[i]]), length(classes[[i]]))
        ds[, i] <- as.numeric(my.d.adj)
        vars[, i] <- as.numeric(sigmad(my.d.adj, sum(classes[[i]] == 
            0), sum(classes[[i]] == 1)))
    }
    sepZscores <- ds/sqrt(vars)
    effects <- ds
    effectsVar <- vars
    colnames(sepZscores) <- paste("zSco_Ex_", 1:num.studies, 
        sep = "")
    colnames(effects) <- paste("Effect_Ex_", 1:num.studies, sep = "")
    colnames(effectsVar) <- paste("EffectVar_Ex_", 1:num.studies, 
        sep = "")
    ds <- ds[, CombineExp]
    vars <- vars[, CombineExp]
    num.studies <- length(CombineExp)
    df <- num.studies - 1
    Qvals <- f.Q(ds, vars)
    if (useREM) 
        vars <- vars + tau2(Qvals, num.studies, my.weights = 1/vars)
    wt <- 1/vars
    MUvals <- rowSums(ds * wt)/rowSums(wt)
    MUsds <- sqrt(1/rowSums(wt))
    zSco <- MUvals/MUsds
    Qpvalues <- 1 - pchisq(Qvals, df)
    Chisq <- 1 - pchisq(zSco^2, 1)
    theResult <- cbind(sepZscores, zSco, MUvals, MUsds, Qvals, 
        df, Qpvalues, Chisq, effects, effectsVar)
    rownames(theResult) <- theNames
    return(theResult)
}