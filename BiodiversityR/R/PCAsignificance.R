`PCAsignificance` <-
function(pca,axes=8) {
    eigen <- pca$CA$eig
    tot <- sum(eigen)
    p <- length(eigen)
    if (p < axes) {axes <- p}
    varexplained <- array(dim=c(7,p))
    varexplained[1,] <- eigen[1:p]
    varexplained[2,] <- varexplained[1,]/tot*100
    varexplained[3,1] <- varexplained[2,1]
    for (i in 2:p) {varexplained[3,i] <- varexplained[3,i-1]+varexplained[2,i]}
    for (i in 1:p) {varexplained[6,i] <- 1/i}
    for (i in 1:p) {varexplained[4,i] <- sum(varexplained[6,i:p])/p*100}
    varexplained[5,1] <- varexplained[4,1]
    for (i in 2:p) {varexplained[5,i] <- varexplained[5,i-1]+varexplained[4,i]}
    for (i in 1:p) {
        if(varexplained[2,i]>varexplained[4,i]) {
            varexplained[6,i] <- TRUE
        }else{
            varexplained[6,i] <- FALSE
        }    
        if(varexplained[3,i]>varexplained[5,i]) {
            varexplained[7,i] <- TRUE
        }else{
            varexplained[7,i] <- FALSE
        }
    }
    rownames(varexplained) <- c("eigenvalue","percentage of variance","cumulative percentage of variance",
        "broken-stick percentage","broken-stick cumulative %","% > bs%","cum% > bs cum%")
    colnames(varexplained) <- c(1:p)
    return(varexplained[,1:axes])
}

