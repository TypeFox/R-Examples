#Core logic same as 3d, modified to work better with 29c
pp <- function(r=0.8,n,data,oth,k)
{
	polar	<- gen_points(n=n,d=k)
	cart	<- angles(polar)
	stopifnot(ncol(data) == ncol(oth))
	n1	<- nrow(data)
	n2	<- nrow(oth)
	full	<- as.matrix(rbind(data,oth))
	ev	<- evaluator(n=n1+n2,p=k)
	proj	<- matrix(0,nrow=ncol(data),ncol=k)
	srank.data	<- matrix(0,nrow=n,ncol=k)
	srank.oth	<- matrix(0,nrow=n,ncol=k)

	index <- function(mat)
	{
		mat	<- matrix(mat,ncol=k)	#sometimes not a matrix!
		proj	<-  t(full %*% mat)	#First n1 columns are for data, next n2 for oth (saves time)
		spmed	<-  trust(ev,parinit=apply(proj,MARGIN=1,FUN=median),samp=t(proj),u=rep(0,k),rinit=0.5,rmax=2e5)
		tmax	<-  max(sqrt(colSums((proj - spmed$argument) ^ 2)))	#max distance of all points from spatial median
		ev.points <- cart * tmax * r
		ev.points <- t(ev.points) - spmed$argument
		for(i in 1:n)
		{
			one	<- proj - ev.points[,i]
			norms	<- sqrt(colSums(one^2))
			srank.data[i,]	<- colMeans(t(one[,1:n1]) / norms[1:n1])
			srank.oth[i,]	<- colMeans(t(one[,(n1+1):(n1+n2)]) / norms[(n1+1):(n1+n2)])
		}
		tmp <- (srank.data - srank.oth)
		vol <- ((sqrt(pi) * tmax) ^ (k/2)) / gamma(k/2 + 1)
		return(mean(sqrt(rowSums(tmp ^ 2))) * vol)
	}
	return(index)
}
