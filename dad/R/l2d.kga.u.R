l2d.kga.u <-
function(x1,x2, check = FALSE)
{
	var1<-var(x1);
	var2<-var(x2);
  if(check)
    {if(var1<.Machine$double.eps |var2<.Machine$double.eps) 
      {stop("At least one variance is zero") 
      }
    }
	n1<-length(x1);
	n2<-length(x2);
	expo<-matrix(0,ncol=n2,nrow=n1);
	w1<-bandwidth.parameter(1,n1);
	w2<-bandwidth.parameter(1,n2);
	vars<-w1^2*var1+w2^2*var2;
	x1<-as.vector(x1);
	x2<-as.vector(x2);
	for (i1 in 1:n1)
		{for (i2 in 1:n2)
			{expo[i1,i2]=exp((-1/2)*(x1[i1]-x2[i2])^2/vars) }};
	sum(expo)/(sqrt(2*pi*vars)*n1*n2)
}
