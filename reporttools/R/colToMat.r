colToMat <- function(tab, cols){

## break a n * p data.frame in a data.frame
## with ceiling(n / cols) rows and
## cols * p columns

tab <- as.matrix(tab)        
n.row <- ceiling(dim(tab)[1] / cols)
n.col <- cols * dim(tab)[2]
mat <- matrix(NA, ncol = n.col, nrow = n.row)
for (c in 1:cols){
    ind.col <- (1 + (c - 1) * dim(tab)[2]):(c * dim(tab)[2])
    ind.row <- (1 + (c - 1) * n.row):(c * n.row)
    ind.row <- unique(pmin(ind.row, dim(tab)[1]))
    up <- n.row
    if (c == cols){
        ind.row <- (1 + (c - 1) * n.row):length(tab[, 1])
        up <- max(ind.row %% n.row)
        if (min(ind.row %% n.row) == 0){up <- n.row}}
    mat[1:up, ind.col] <- tab[ind.row, ]
}

mat <- as.data.frame(mat)
dimnames(mat)[[2]] <- rep(dimnames(tab)[[2]], cols)

return(mat)
}
