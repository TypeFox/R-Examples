coupler.classes.adj.i <-
function(cl, mat)
{
    return(names(which(mat[,cl] == 1)))
}
