ISIMP <-
function(a,b,f)
{

# if an even number take the first obs out

n <- (length(f)-1)/2
sum <- f[1] + 4*f[2*n] + f[2*n+1]

for (i in 2:(2*n-1))
{    
    if (i/2 - as.integer(i/2) == 0)
    c <- 4
    else
    c <- 2
sum <- sum + c*f[i]
}

h <- (b - a)/(2*n)
sum <- h*sum/3
return(sum)
}
