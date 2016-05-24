scaledata <-
function(mat, scale.Y=FALSE, col.means, col.norms)
{
 n <- NROW(mat)
  k <- NCOL(mat)
  cnms <- colnames(mat)
    col.means <- colMeans(mat)
  matcenter <- scale(mat, center=col.means, scale=FALSE)
     col.norms <- sqrt(colSums(matcenter*matcenter))
  matscale <- scale(matcenter, center=FALSE, scale=col.norms)
list(Y=matcenter,X=matscale)
}
