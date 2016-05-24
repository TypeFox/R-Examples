VAR.LRtest <-
function(x,p,restrict0,restrict1,type="const")
{
options(warn=-1)
if(type=="none") add <- 0
if(type=="const") add <- 1
if(type=="const+trend") add <-2
k <- ncol(x)

var1 <- VAR.est(x,p,type)
n = nrow(x)
b <- var1$coef
sigu <- var1$sigu
{
if(restrict1=="full") {

VARmat = VAR.Rmat(p,k,restrict0,type)
Rmat2 = VARmat$Rmat; rmat2 = VARmat$rvec
sigu1 = var1$sigu *((n-p)-ncol(b))/(n-p); nor=nrow(restrict0)
sigu2=VAR.Rtem(x,p,sigu,Rmat2,rmat2,type)$sigu
df3=ncol(var1$coef)}
else {

VARmat = VAR.Rmat(p,k,restrict1,type)
Rmat1 = VARmat$Rmat; rmat1 = VARmat$rvec

VARmat = VAR.Rmat(p,k,restrict0,type)
Rmat2 = VARmat$Rmat; rmat2 = VARmat$rvec

sigu2=VAR.Rtem(x,p,sigu,Rmat1,rmat1,type)$sigu 
var2=VAR.Rtem(x,p,sigu,Rmat2,rmat2,type)
sigu1=var2$sigu
nor=nrow(restrict0)-nrow(restrict1); df3=ncol(var2$coef)
      }
}
LR=(n-p) * abs( log(det(sigu2)) - log(det(sigu1)) ); LRpval=1-pchisq(LR,df=nor)
Fstat=LR/nor; Fpval = 1-pf(Fstat,df1=nor,df2=(nrow(x)-p)-df3)
return(list(LRstat=LR,LRpval=LRpval,Fstat=Fstat,Fpval=Fpval))
}
