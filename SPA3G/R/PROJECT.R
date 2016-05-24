PROJECT <-
function(MM)
{
n <- nrow(MM)
PMM <- MM-(matrix(1, n, 1)%*%apply(MM, 2, sum))/n
return(PMM)
}
