rstinterval <- function(pt, beta=c(0.90), sampleSize=c(250), M=c(50),clr=FALSE, ...)
{
	if(!is.matrix(pt)) stop ("first argument pt must be a 2 dimensional matrix");
	n=length(pt)/2
	m=n
	if(m<sampleSize) stop ("we expect a matrix of size more than sampleSize or 250");



 plot(pt[,1],pt[,2],pch=20,xlab='x',ylab='y',col="#747170");
 
 #pick 
 x=pt[,1]
 y=pt[,2]

 index = floor(runif(sampleSize,1,m))
 
 sampleX=x[index]
 sampleY=y[index]
 
       m1 = .C("rs_depthrings",
                   xpoints=as.double(sampleX), 
  		   ypoints=as.double(sampleY),
  		   outx=double(sampleSize),
  		   outy=double(sampleSize),
                   sz=as.integer(sampleSize),
                    PACKAGE = "rsdepth");

betaPoint = cbind( m1$outx[floor(sampleSize-sampleSize*beta)],m1$outy[floor(sampleSize-sampleSize*beta)] )
depth = rsdepth(pt,betaPoint)

count = 1;
outX=rnorm(M*sampleSize)
outY=rnorm(M*sampleSize)
for(i in 1:M)
{

	index = floor(runif(sampleSize,1,m))
	 
	sampleX=x[index]
	sampleY=y[index]
	
	for(j in 1:sampleSize)
	{
		iPoint = cbind( sampleX[j],sampleY[j] )
		if( rsdepth(pt,iPoint)>=depth )
		{
			outX[count] = sampleX[j]
			outY[count] = sampleY[j]
			count=count+1;
		}
		
	} 
#	print(i);
#	print(count);
}

outx=outX[1:(count-1)]
outy=outY[1:(count-1)]
k=chull(outx,outy);
#print(outx)
#print(outy)
polygon(outx[k],outy[k],...);

return (0);   

}

