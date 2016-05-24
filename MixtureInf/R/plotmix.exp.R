plotmix.exp <-
function(x,theta,hist=1,comp=TRUE,k=20,h=1,main="",xlab="Observations",ylab="")
#x: 		data, can be either a vector or a matrix with the 1st column being the observed values 
#        		and the 2nd column being the corresponding frequencies. 
#theta:  	output of pmle.pois or emtest.pois, or a vector of parameter values, include mixing proportions and mixing parameters.
#hist:	style of histogram. We can select to obtain two styles of histogram (hist=1 or 2). hist=0 means no histogram.
#comp:	parameter for compontent fitted density. 
#		comp=T means component fitted probability mass function are drawn, and comp=F means no component fitted probability mass function.
#k:		a single number giving the number of cells for the histogram.
#main: 	title of graph.
#xlab: 	label of x-axis.
#ylab: 	label of y-axis.
{
	if (is.data.frame(x))
	{	
		if (ncol(x)==2)
			x=as.matrix(x)
		if (ncol(x)==1 | ncol(x)>2)
			x=x[,1]
	}
	a=c()
	b=c()
	if (is.matrix(x))
	{
		x=x[order(x[,1]),]
		a=x[,1]
		b=x[,2]
	}
	if (is.vector(x))
	{
		n=length(x)
		x=sort(x)
		r=(x[n]-x[1])/k
		for (i in 1:k)
		{
			a[i]=x[1]+(i-1/2)*r
			b[i]=sum(x>=a[i]-1/2*r & x<a[i]+1/2*r)
		}
		b[k]=sum(x>=a[i]-1/2*r & x<=a[i]+1/2*r)
	}

	l=length(a)
	d=(a[l]-a[1])/(l-1)/2
	b=b/sum(b)/d/2
	b=c(0,b)

	if (is.list(theta))
	{
		if (is.vector(theta[[1]]))
			theta=c(theta[[1]],theta[[2]])
		if (is.matrix(theta[[1]]))
			theta=c(theta[[1]][1,],theta[[1]][2,])
	}

	m=length(theta)/2
	xx=seq(a[1]-d,a[l]+d,(a[l]-a[1]+2*d)/100)
	dd=matrix(0,m,length(xx))
	ee=matrix(0,1,length(xx))
	for (j in 1:m)
	{
		dd[j,]=theta[j]*dexp(xx,1/theta[m+j])
		ee=ee+dd[j,]
	}
	hh=1.05*max(c(b,ee))
	plot(c(a[1]-d,a[l]+d),c(0,h*hh),"n",main=main,xlab=xlab,ylab=ylab)
	if (hist==1)
	{
		lines(c(a[1]-d,a[l]+d),c(0,0),"l",col="black")
		for (i in 1:l)
		{
			lines(c(a[i]-d,a[i]+d),c(b[i+1],b[i+1]),"l",col="blue")
			lines(c(a[i]-d,a[i]-d),c(b[i],b[i+1]),"l",col="blue")
		}
		lines(c(a[l]+d,a[l]+d),c(b[l+1],0),"l",col="blue")
	}
	if (hist==2)
	{
		lines(c(a[1]-d,a[l]+d),c(0,0),"l",col="black")
		for (i in 1:l)
		{
			lines(c(a[i]-d,a[i]+d),c(b[i+1],b[i+1]),"l",col="blue")
			lines(c(a[i]-d,a[i]-d),c(0,b[i+1]),"l",col="blue")
			lines(c(a[i]+d,a[i]+d),c(0,b[i+1]),"l",col="blue")
		}
	}

	lines(xx, ee, lty=1)
	if (comp==T)
	{
		for (j in 1:m)
		lines(xx, dd[j,], lty=2)
	}
}
