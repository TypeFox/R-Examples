"Make.Z" <-
function(x)
{

    x <- as.matrix(x)
    x <- sweep(x, 2, colMeans(x, na.rm=TRUE))
zx <- sweep(x,2,apply(x, 2, sd, na.rm=TRUE),"/")
return(zx)
}

