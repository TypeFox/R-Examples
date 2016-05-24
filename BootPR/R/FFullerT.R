FFullerT <-
function(x,p,h)
{
x <- as.matrix(x)
n <- nrow(x)
b <- fullerT(x,p)
e <- RESIDT(x,b)
f <- {}
if(h > 0)
f <- ART.Fore(x,b,h)
return(list(coef=b,resid=e,forecast=f))
}
