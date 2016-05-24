Shaman.Stine <-
function(x,p,a)
{
x <- as.matrix(x)
n <- nrow(x)
#b <- LSM(x,p)$coef
b=matrix(a)
tem = Stine(b,n,p); ahat0=tem$ahat.old
bc <- c(tem$ahat,b[p+1])
bc <- adjust(b,bc,p)
bc[p+1] <- mean(x)*(1-sum(bc[1:p]))
e <- RESID(x,bc)

return(list(coef=bc,coef0=ahat0,resid=e,mat=tem$mat))
}
