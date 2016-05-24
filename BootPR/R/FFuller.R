FFuller <-
function(x,p,h)
{
x <- as.matrix(x)
n <- nrow(x)
b <- fuller(x,p)
e <- RESID(x,b)
f <- {}
if(h > 0)
f <- AR.Fore(x,b,h)
return(list(coef=b,resid=e,forecast=f))
}
