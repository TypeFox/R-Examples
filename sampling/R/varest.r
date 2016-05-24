varest<-function(Ys,Xs=NULL,pik,w=NULL)
{
 if (any(is.na(pik))) 
        stop("there are missing values in pik")
 if (any(is.na(Ys))) 
        stop("there are missing values in y")
 if (length(Ys) != length(pik)) 
        stop("y and pik have different sizes")
if(!is.null(Xs))
 {if(is.data.frame(Xs)) Xs=as.matrix(Xs)
  if(is.vector(Xs) & (length(Ys)!= length(Xs)))  
        stop("x and y have different sizes")
  if(is.matrix(Xs) & (length(Ys) != nrow(Xs)))  
        stop("x and y have different sizes")
}
 a=(1-pik)/sum(1-pik)
 if(is.null(Xs))
	{A=sum(a*Ys/pik)
	 var=sum((1-pik)*(Ys/pik-A)^2)/(1-sum(a^2))
	}
 else 
   {B=t(Xs*w) 
   beta=ginv(B%*%Xs)%*%B%*%Ys
   e=Ys-Xs%*%beta
   A=sum(a*e/pik)
   var=sum((1-pik)*(e/pik-A)^2)/(1-sum(a^2))
}
var
}