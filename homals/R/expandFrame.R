expandFrame<-function(tab, clean = TRUE, zero = TRUE, returnFrame = TRUE) 
{
# expandFrame expands a matrix or data-frame to an indicator supermatrix and 
# optionally converts this to a data-frame again. By default NA becomes zero 
# and constant rows and columns are eliminated.

n<-dim(tab)[1]
m<-dim(tab)[2]
g<-matrix(0,n,0)
l<-rep("",0)
lab1<-labels(tab)[[1]]
lab2<-labels(tab)[[2]]
for (j in 1:m) {
	y<-as.factor(tab[,j]); h<-levels(y)
	g<-cbind(g,ifelse(outer(y,h,"=="),1,0))
	l<-c(l,paste(lab2[j],"_",h,sep=""))
	}
if (zero) g<-ifelse(is.na(g),2,as.matrix(g))
if (clean) {
	g<-g[which(rowSums(g)>0),which(colSums(g)>0)]
	g<-g[,which(colSums(g)<n)]
	}
if (!returnFrame) return(g)
g<-as.data.frame(g,row.names=lab1)
names(g)<-l
return(g)
}