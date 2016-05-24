lt2full <-
function(x){
nn <- getNumNodes(x, type="adjmatrixlt")
y <- matrix(0, nn, nn)
y[lower.tri(y)] <- x
y <- y + t(y)

return(as.vector(y))
}
