"matcor" <-
function (X, Y) 
{
    matcorX = cor(X, use = "pairwise")
    matcorY = cor(Y, use = "pairwise")
    matcorXY = cor(cbind(X, Y), use = "pairwise")
    return(list(Xcor = matcorX, Ycor = matcorY, XYcor = matcorXY))
}

