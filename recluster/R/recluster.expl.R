recluster.expl<-function(mat,clust){
mat<-as.matrix(mat)
beta <- sum(mat)
cluster <- NULL
betapartial <- 0
        for (row in 1:nrow(mat)) {
            for (col in 1:ncol(mat)) {
                if (clust[row] != clust[col]) {
                  betapartial <- betapartial + mat[row, col]
                }
            }
        }
return(betapartial/beta)
}
