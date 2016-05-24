"alpha" <-
function(x)
{
x <- na.exclude(as.matrix(x))
Sx <- sum(var(x))
SumSxi <- sum(apply(x,2,var))
k <- ncol(x)
alpha <- k/(k-1)*(1-SumSxi/Sx)
return(alpha)
}

