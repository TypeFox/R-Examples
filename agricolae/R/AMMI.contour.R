`AMMI.contour` <-
function(model,distance,shape,...)
{
G<- subset(model$biplot,model$biplot$type=="GEN")
x<-G$PC1
y<-G$PC2
d<-sqrt(x^2+y^2)
r <-distance*max(d)
x<-seq(-r,r,length=shape)
A<-cbind(x,y=sqrt(r^2-x^2))
B<-cbind(x,y=-sqrt(r^2-x^2))
lines(A,type="l",...)
lines(B,type="l",...)
Gin <- d<=r
Gout<- d >r
GEN.in<-row.names(G)[Gin]
GEN.out<-row.names(G)[Gout]
cat("\nLimit, radio:",r)
cat("\nGenotype  in:",length(GEN.in))
cat("\nGenotype out:",length(GEN.out),"\n\n")
distance<-data.frame(row.names=row.names(G),distance=d)
return(list("GENOTYPE IN"=GEN.in, "GENOTYPE OUT"=GEN.out,Distance=distance))
}

