hfp <- function (obj, gen, ind, Y){
    Y <- as.matrix(Y)
    genes <- row.names(Y)
    indxs <- colnames(Y)
    num1 <- length(gen)
    num2 <- length(ind)
    gen.index <- rep(0, num1)
    ind.index <- rep(0, num2)
    for (i in 1:num1) gen.index[i] <- which(genes == gen[i])
    for (j in 1:num2) ind.index[j] <- which(indxs == ind[j])
    W_imp <- obj$PLS.imp
    if (max(W_imp) != 0) 
        heatmap(t(W_imp[gen.index, ind.index]), Rowv = NA, Colv = NA, 
            labRow = as.character(ind), labCol = gen, xlab = "Specified genes", ylab = "Specified subjects")
    if (max(W_imp) == 0) 
        print("No trace of hidden subject-specific variation is found in the data")
}
