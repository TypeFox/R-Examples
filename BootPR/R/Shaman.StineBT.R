Shaman.StineBT <-
function(x,p,h)
{
x <- as.matrix(x)
n <- nrow(x)
b <- LSMBT(x,p)$coef
bc <- c(StineT(b,n,p),b[p+1],b[p+2])
bc <- adjust(b,bc,p)

bc[(p+1):(p+2),] <- RE.LSMTB(x,p,bc)
e <- RESIDT(x,bc)
f <- {}
if(h > 0)
f <- ART.Fore(x,bc,h)
return(list(coef=bc,resid=e,forecast=f))
}
