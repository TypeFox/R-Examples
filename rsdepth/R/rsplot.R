rsplot <- function(x, y=NULL, factorsecondbag=2, mring=T,...)
{
	len = length(x)	
	if( is.null(y) )	{
    		if(!is.matrix(x)) stop ("first argument pt must be a 2 dimensional matrix");
                   x1=as.double(x[,1])
		   y=as.double(x[,2])
		   x=x1;
		   len=len/2;
	}
	pt=cbind(x,y);
	m=len;

 plot(x,y,pch='.',...)
       m1 = .C("rs_depthrings",
                   xpoints=as.double(pt[,1]), 
  		   ypoints=as.double(pt[,2]),
  		   outx=double(len),
  		   outy=double(len),
                   sz=as.integer(len),
                    PACKAGE = "rsdepth");

	x=m1$outx[(1+m/2):m]
	y=m1$outy[(1+m/2):m]

	cbag=convexhull(x,y)
	
	#y=matrix(y,nc=2);
	inflatedcbag=inflate(cbag,factor=factorsecondbag);

	x=pt[,1]
	y=pt[,2]
	k = chull(pt[,1],pt[,2]);


chbag = .C("polygonintersection",
            p1X=as.double(x[k]), p1Y=as.double(y[k]), 
            p2X=as.double(inflatedcbag[,1]), p2Y=as.double(inflatedcbag[,2]), 
            p3X=double( 2*(length(k)+length(inflatedcbag)) ), p3Y=double( 2*(length(k)+length(inflatedcbag)) ), 
            n=as.integer(length(k) ),
            m=as.integer(length(inflatedcbag[,1]) ),
            p=integer(c(1)),
            PACKAGE = "rsdepth");

k= chull( chbag$p3X[1:chbag$p], chbag$p3Y[1:chbag$p] );
#print(chbag$p3X)
#print("two")
#print(chbag$p3X[k])
#print("three")
#print(k)
polygon(chbag$p3X[k],chbag$p3Y[k],col="#737CA1")
polygon(cbag[,1],cbag[,2],col="#43C6DB")

if(mring==T)
{

	rtvalue = .C("rs_getcenter",
                    ptX=as.double(x), ptY=as.double(y), 
                    num=as.integer(len),
                    center=double( (len)^4 ), cSize=integer(c(1)),
                    PACKAGE = "rsdepth");

	center = matrix(rtvalue$center[1:(rtvalue$cSize)],ncol=2,byrow=TRUE)
	cent=convexhull(center[,1],center[,2]);

	#print(rtvalue$cSize)
	polygon(cent[,1],cent[,2],col="#151B8D");
	md=centroid(cent[,1],cent[,2])
	
	points(x,y,pch=20,...)	
	points(md[1],md[2],pch=20,col='red')

}
else
{

points(x,y,pch=20,...)	
points(m1$outx[m],m1$outy[m],pch=20,col="red")
}

}

