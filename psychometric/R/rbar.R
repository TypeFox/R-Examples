"rbar" <-
function(x)
{
rxy <- x$Rxy
n <- x$n
rbar <- sum(n*rxy, na.rm=TRUE)/sum(n,na.rm=TRUE)
return(rbar)
}

