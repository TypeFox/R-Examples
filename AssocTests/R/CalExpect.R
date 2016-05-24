CalExpect <- function(H)
{
    a <- apply(H, 1, sum)
    b <- apply(H, 2, sum)
    n <- sum(H)
    
    kronecker(a,b)/n
}
