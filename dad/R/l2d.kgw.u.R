l2d.kgw.u <-
function(x1,varw1,x2,varw2, check = FALSE)
{
  if(check)
    {if(abs(det(varw1))<.Machine$double.eps | abs(det(varw2))<.Machine$double.eps )
      {stop("One of the bandwidths is zero")
      }
    }  
  n1=length(x1)
  n2=length(x2)
	p<-1;
	expo<-matrix(0,ncol=n2,nrow=n1);
	varsom<-varw1+varw2;
	varinv<-1/varsom;
	x1<-as.vector(x1);
	x2<-as.vector(x2);
	for (i1 in 1:n1)
		{for (i2 in 1:n2)
			{expo[i1,i2]=exp((-1/2)*(x1[i1]-x2[i2])*varinv*(x1[i1]-x2[i2])) }};
	(1/(n1*n2))*(1/(2*pi)^(p/2))*(1/varsom^(1/2))*sum(expo)
}
