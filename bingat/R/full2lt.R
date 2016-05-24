full2lt <-
function(x){
nn <- getNumNodes(x, type="adjmatrix")
y <- matrix(x, nn, nn, byrow=FALSE)
y <- y[lower.tri(y)] 

return(as.vector(y))
}
