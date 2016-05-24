"makematrix" <-
function (n) 
{
    x <- runif(n)
    y <- runif(n)
    Ag <- outer(x, x, "-")^2 + outer(y, y, "-")^2
    Ag <- (2 - Ag)/10
    diag(Ag) <- 0
    diag(Ag) <- 1/n + apply(Ag, 2, sum)
    Ag
}
