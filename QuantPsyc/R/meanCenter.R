"meanCenter" <-
function (x)
{
mcx <- x - mean(x, na.rm=TRUE)
return(mcx)
}

