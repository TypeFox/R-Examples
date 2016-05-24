"pcochran" <-
function(q,n,k)
{

f <- (1/q-1)/(k-1)
p <- 1 - pf(f,(n-1)*(k-1),n-1)*k

p[p<0] <- 0
p[p>1] <- 1

return(p)

}

