thinl<-function(coord, tol)
{
	tol2<-tol*tol

	x<-coord[[1]]
	y<-coord[[2]]

	npoints<-length(x)
	index<-c(1)
	j<-1
	
	for(i in 2:(npoints-1))
	{
		xx<-x[j]-x[i]
		yy<-y[j]-y[i]
		if((xx*xx+yy*yy)>tol2)
		{
			j<-i
			index<-c(index,i)
		}

	}

	index<-c(index,npoints) #Suggested by Erich Neuwirth	

	list(x[index], y[index])
}

thinlines<-function(arc, tol)
{
	newarc<-list()

	narcs<-length(arc[[1]][[1]])

	newarc<-lapply(arc[[2]],thinl, tol=tol)

	newtable<-arc[[1]]

	newtable$NVertices<-as.numeric( lapply(newarc, function(X){length(X[[1]])})  )
	
	list(newtable, newarc)
}
