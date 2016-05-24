`mstlines` <-
function(mst, x, y = NULL, pts.names = NULL, ...) {
  a <- colnames(mst)
  if (ncol(as.matrix(x)) >= 2 & nrow(as.matrix(x)) >= 2 & is.null(y) == TRUE) b <- x[,1:2]
  else b <- cbind(x,y)
  if (is.null(pts.names) == FALSE) rownames(b) <- pts.names
  else rownames(b) <- colnames(mst)
  v <- colnames(mst)
    for (i in 1:nrow(mst)) {
        for (j in 1:ncol(mst)) {
            if (mst[i, j] == 1) {
                lines(c(b[v[i], 1], b[v[j],1]), c(b[v[i], 2], b[v[j],2]), ...)
            }
        }
    }
}
