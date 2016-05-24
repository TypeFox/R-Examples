getbag <- function(x, y=NULL, factorsecondbag=2,...)
{
	len = length(x)	
	if( is.null(y) )	{
    		if(!is.matrix(x)) stop ("first argument pt must be a 2 dimensional matrix");
                   x1=as.double(x[,1])
		   y=as.double(x[,2])
		   x=x1;
		   len=len/2;
	}
	pt=cbind(x,y)
	bagi=double( (len)^4 );
	numi=(len);

#x1=rbind(x,c(-1))
#y1=rbind(y,c(-1))
plot(x,y,pch='.',...)

#drawcompletegraph(x,y,start=F)        

	rtvalue = .C("rs_getbag",
                    ptX=as.double(x), ptY=as.double(y), 
                    bag=double( (len)^4 ), num=as.integer(numi),
                    bagsz=integer(c(1)),
                    center=double( (len)^4 ), cSize=integer(c(1)),
                    PACKAGE = "rsdepth");

##truncate bag (2d matrix byrow) to up till num/2
bag = matrix(rtvalue$bag[1:(rtvalue$bagsz)],ncol=2,byrow=TRUE)
center = matrix(rtvalue$center[1:(rtvalue$cSize)],ncol=2,byrow=TRUE)
##points(rtvalue$bagx,rtvalue$bagy,pch='x')
##return (bag);
##print(bag)
##find convex hull
cbag=convexhull(bag[,1],bag[,2]);

#y=matrix(y,ncol=2);
icbag=inflate(cbag,factor=factorsecondbag);

k = chull(x,y);

##polygon(icbag[,1],icbag[,2],col="red")
##polygon(x[k],y[k],col="green")
##take intersection of bag with convex hull

chbag = .C("polygonintersection",
            p1X=as.double(x[k]), p1Y=as.double(y[k]), 
            p2X=as.double(icbag[,1]), p2Y=as.double(icbag[,2]), 
            p3X=double( 2*(length(k)+length(icbag)) ), p3Y=double( 2*(length(k)+length(icbag)) ), 
            n=as.integer(length(k) ),
            m=as.integer(length(icbag[,1]) ),
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

cent=convexhull(center[,1],center[,2]);

#print(rtvalue$cSize)
polygon(cent[,1],cent[,2],col="#151B8D")

points(pt[,1],pt[,2],pch=20)

#md=rsmed(pt)
md=centroid(cent[,1],cent[,2])

points(md[1],md[2],pch='*',col='red')

	return (y);

}

