rST <-
function(n, df=10, skew=1)
{
zstar <- rt(n=n,df=df)
weight <- skew/(skew + 1/skew)
z <- runif(n, -weight, 1 - weight)
signz <- sign(z)
epsilon <- skew^signz
result <- -abs(zstar)/epsilon * signz
return(result)
}
