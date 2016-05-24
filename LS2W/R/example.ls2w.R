`example.ls2w`<-function(n=25,size=64)
{
data(textures,envir=environment())

# the lines below (and the x=Aim variable call in sample.stats) make for cleaner
# local variable manipulation (it avoids confusion if there are multiple variables
# of the same name in different environments.

Aim<-get("A",envir=environment())
Bim<-get("B",envir=environment())
Cim<-get("C",envir=environment())

A.stats<- sample.stats(x=Aim,n,size)
B.stats<- sample.stats(x=Bim,n,size)
C.stats<- sample.stats(x=Cim,n,size)
all.stats<-rbind(A.stats, B.stats, C.stats)
imlabels<-c(rep("a",n), rep("b",n),rep("c",n))
all.stats.lda<-lda(all.stats, imlabels)
all.stats.ld<-predict(all.stats.lda, dimen=2)$x
plot(all.stats.ld, type="n", xlab="First Linear Discriminant", ylab="Second Linear Discriminant")
text(all.stats.ld, imlabels)
return(all.stats.lda)
}
