rsmed <- function(pt, eps=c(0), ...)
{
	if(!is.matrix(pt)) stop ("first argument pt must be a 2 dimensional matrix");
	if(eps>1) stop ("approximation argument eps must be greater than 0 and less than 1");
	if(eps<0) stop ("approximation argument eps must be >= 0 and <= 1");
	n=length(pt)/2
	m=n
	if(eps>0)
	{	
		m = ceiling(1/eps^2*log(1/eps))
		if(m>n)
		{
			m=n	
		}
		else
		{
	
			for(i in 1:m)
			{
				##generate a random integer
				id = floor(runif(1,i,n+1))
				##swap point at id with point at 
				temp = pt[id,1]
				pt[id,1]=pt[i,1]
				pt[i,1]=temp
	
				temp = pt[id,2]
				pt[id,2]=pt[i,2]
				pt[i,2]=temp	
	
			}	
			print(c("approximating median by using a sample of size :"))
			print(m[1])
	
		}
	
	
	}

        median = .C("rs_med",
                   xpoints=as.double(pt[,1]), 
  		   ypoints=as.double(pt[,2]), 
                    md=double(2), 
                    sz=as.integer(m),
                    PACKAGE = "rsdepth")$md;
                    
	return (median);

}

