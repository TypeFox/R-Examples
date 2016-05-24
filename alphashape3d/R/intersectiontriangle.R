intersectiontriangle <-
function (triangles, coord, points) 
{
    points = as.matrix(points, nc = 3)
    inside <- integer(dim(points)[1])
    tempdir = as.numeric(2 * runif(3) - 1)
    retour <- .C("pointinashape", as.integer(triangles), dim(triangles)[1], 
        as.numeric(coord), dim(coord)[1], as.numeric(points), 
        dim(points)[1], as.integer(inside), tempdir)
    return(retour[[7]])
}
