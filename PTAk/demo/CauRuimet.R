#########
demo.CauRuimet <- function(ker=0.84,withingroup=TRUE,openX11s=FALSE)
#########
{
data(iris)
iris2 <- as.matrix(iris[,1:4])
dimnames(iris2)[[1]] <- iris[,5]
D2 <- CauRuimet(iris2,ker=ker,withingroup=withingroup)
D2 <- Powmat(D2,(-1))
iris2 <- sweep(iris2,2, apply(iris2,2,mean))
res <- SVDgen(iris2,D2=D2,D1=1)
res2 <- SVDgen(iris2)
if(openX11s)X11(width=6,height=4)
par(mfrow=c(1,2))

if(withingroup)plot(res,nb1=1,nb2=2,cex=0.5,main=expression(paste("Unknown ",Wg^{-1}) ),type="n")
else plot(res,nb1=1,nb2=2,cex=0.5,main=expression(paste(Wo^{-1})),type="n")
 mtext(paste("ker= ",ker),side=3,cex=1,line=0)
plot(res2,nb1=1,nb2=2,cex=0.5,main="Id",type="n")
summary(res,testvar=0)
cat("\n","SVD with metrics in variables space","\n")
browser()
print(demo.CauRuimet)
invisible(res)
}
######

demo.CauRuimet()
cat("\n","\n","args(demo.CauRuimet)","\n")
print(args(demo.Cauruimet))
