zScorePermuted <-function (esets, classes, useREM = TRUE, CombineExp = 1:length(esets)) 
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
    zScores(esets, lapply(classes, sample), useREM, CombineExp = CombineExp)[, 
        "zSco"]
}
