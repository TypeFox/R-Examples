data(iris)
Y  = iris[,-5]
lbl= as.numeric(iris[,5])
K  = 3
set.seed(2)
Tinit = t(rmultinom(150,1,c(rep(1/K,K))))
U  = fem(Y,3,init='user',Tinit=Tinit,model='AkB',maxit=1,method='gs',eps=1e-3)$U
cls1 = cls = max.col(Tinit)
for (i in 1:15){
	  # x11()
	  # def.par <- par(no.readonly = TRUE)
	  x = as.matrix(Y)%*%U
	  min1= round(min(x[,1]),1)-1
	  max1= round(max(x[,1],1))+1
	  min2= round(min(x[,2]),1)-1
	  max2= round(max(x[,2]),1)+1
	  cls[cls==1]='lightgreen'
	  cls[cls==2]='gray'
	  cls[cls==3]='darkblue'
	  x1 = split(x[,1],cls)
	  x2 = split(x[,2],cls)

	  xhist1 <- hist(x1$'lightgreen',breaks=seq(min1,max1,0.1), plot=FALSE) 
	  xhist2 <- hist(x1$'gray',breaks=seq(min1,max1,0.1), plot=FALSE) 
	  xhist3 <- hist(x1$'darkblue',breaks=seq(min1,max1,0.1), plot=FALSE) 
	  yhist1 <- hist(x2$'lightgreen',breaks=seq(min2,max2,0.1), plot=FALSE) 
	  yhist2 <- hist(x2$'gray',breaks=seq(min2,max2,0.1), plot=FALSE) 
	  yhist3 <- hist(x2$'darkblue',breaks=seq(min2,max2,0.1), plot=FALSE) 
	  topX <- max(c(xhist1$counts,xhist2$counts,xhist3$counts))
	  topY <- max(c(yhist3$counts,yhist3$counts,yhist3$counts))
	  nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE) 
	  xrange <- c(min1,max1)
	  yrange <- c(min2,max2)

	  par(mar=c(3,3,1,1)) 
	  plot(x,col=cls,pch=cls1,xlim=xrange, ylim=yrange,cex=1.5)
	  # , xlim=xrange, ylim=yrange, xlab="", ylab="") 
	  par(mar=c(0,3,1,1)) 
	  barplot(xhist1$counts, axes=FALSE,col='lightgreen',ylim=c(0,topX), space=0) 
	  barplot(xhist2$counts, axes=FALSE,col='gray', ylim=c(0,topX), space=0,add=T)
	  barplot(xhist3$counts, axes=FALSE,col='darkblue', ylim=c(0,topX), space=0,add=T)
	  par(mar=c(3,0,1,1)) 
	  barplot(yhist1$counts, axes=FALSE,col='lightgreen', xlim=c(0,topY), space=0, horiz=TRUE) 
	  barplot(yhist2$counts, axes=FALSE,col='gray', xlim=c(0,topY), space=0, horiz=TRUE,add=T) 
	  barplot(yhist3$counts, axes=FALSE,col='darkblue', xlim=c(0,topY), space=0, horiz=TRUE,add=T) 
	  res = fem(Y,3,init='user',Tinit=Tinit,model='AkB',maxit=1,method='gs',eps=1e-3)
	  Tinit = res$P
	  cls1 = cls = res$cls
	  U = res$U
}
