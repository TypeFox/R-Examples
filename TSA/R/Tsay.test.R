`Tsay.test` <-
function(x,order,...){
if(missing(order)) order=ar(x,...)$order 
x=as.vector(x)
m=order
X=NULL
for (i in 1:m) X=cbind(X,zlag(x,i))
X1=cbind(x,X)
for (j in 1:m) { for (k in 1:j) X=cbind(X,zlag(x,j)*zlag(x,k))}
X2=cbind(x,X)
X1=na.omit(X1)
X2=na.omit(X2)
y=X1[,1]
lm1=lm(y~X1[,-1])
lm2=lm(y~X2[,-1])
a1=anova(lm1,lm2)
list(test.stat=signif(a1[[5]][2],4),p.value=signif(a1[[6]][2],4),order=order)
}

