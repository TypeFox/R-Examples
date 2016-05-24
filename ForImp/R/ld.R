ld <-
function(mat)
{
mat[-which(is.na(mat),arr.ind=TRUE)[,1],]
}

