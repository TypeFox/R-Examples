randomcc <-
function(x, fraction=0.1)
{
size <- nrow(x)
newSize <- round(size * fraction)
randomRows <- sample(size, newSize)
randomObjects <- x[randomRows,]
return(randomObjects)
}
