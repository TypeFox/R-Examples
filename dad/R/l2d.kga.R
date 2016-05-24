l2d.kga <-
function(x1,x2, check = FALSE)
{
	var1<-var(x1);
	var2<-var(x2);
  if(check)
   {if(abs(det(var1))<.Machine$double.eps | abs(det(var2))<.Machine$double.eps)
    {stop("One of the sample variances is degenerate")
    }
   }
	p<-ncol(x1);
	n1<-nrow(x1);
	n2<-nrow(x2);                                             
	expo<-matrix(0,ncol=n2,nrow=n1);
  w1<-bandwidth.parameter(p,n1);
	w2<-bandwidth.parameter(p,n2);
	vars<-w1^2*var1+w2^2*var2;
	varinv<-solve(vars);
	x1<-as.matrix(x1);
	x2<-as.matrix(x2);
	for (i1 in 1:n1)
		{for (i2 in 1:n2)
			{expo[i1,i2]=exp((-1/2)*(x1[i1,]-x2[i2,])%*%varinv%*%(x1[i1,]-x2[i2,])) 
      }
    }
	(1/(n1*n2))*(1/(2*pi)^(p/2))*(1/det(vars)^(1/2))*sum(expo)
}

