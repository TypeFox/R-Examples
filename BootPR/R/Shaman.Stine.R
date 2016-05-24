Shaman.Stine <-
function(x,p,h)
{
x <- as.matrix(x)
n <- nrow(x)
b <- LSM(x,p)$coef
bc <- c(Stine(b,n,p),b[p+1])

    bc <- adjust(b,bc,p)
    
bc[p+1] <- mean(x)*(1-sum(bc[1:p]))
e <- RESID(x,bc)
f <- {}
if(h > 0)
f <- AR.Fore(x,bc,h)
return(list(coef=bc,resid=e,forecast=f))
}
