rsrings <- function(pt, numofrings=c(5), clr=FALSE, ...)
{
	if(!is.matrix(pt)) stop ("first argument pt must be a 2 dimensional matrix");
	n=length(pt)/2
	m=n


 plot(pt[,1],pt[,2],pch=20,xlab='x',ylab='y');
       m1 = .C("rs_depthrings",
                   xpoints=as.double(pt[,1]), 
  		   ypoints=as.double(pt[,2]),
  		   outx=double(m),
  		   outy=double(m),
                   sz=as.integer(m),
                    PACKAGE = "rsdepth");

##don't return the convex hull... start from 2.
numofrings=numofrings+1
	for(i in 2:numofrings)
	{
	
		x=m1$outx[(1+( (i-1)*m/numofrings)):m]
		y=m1$outy[(1+( (i-1)*m/numofrings)):m]
#		print(m1$outx)
#		print(m1$outy)
	
		k=chull(x,y)
		polygon(x[k],y[k])
	}
	points(m1$outx[m],m1$outy[m],pch=20,col="red")
	
	


}
