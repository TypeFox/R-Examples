"UPpivotal" <-
function(pik,eps=1e-6)
{
if(any(is.na(pik))) stop("there are missing values in the pik vector")
N<-length(pik)
s<-rep(0,times=N)
a<-pik[1]
b<-pik[2]
i<-1
j<-2
k<-3
while(k<=N)
{
u<-runif(1)
if(a>=eps & a<= 1-eps & b>=eps & b<= 1-eps)
if(a+b>1)
	{
		if(u<(1-b)/(2-a-b)) 
			{b<-a+b-1;a<-1} 
		else {a<-a+b-1;b<-1} 
      }
	else{ if(u< b/(a+b)) 
			{b<- a+b;a<-0} 
            else {a<- a+b;b<-0} 
           }
      if( (a<eps | a > 1-eps)& (k<=N)) 
                {s[i]=a;a=pik[k];i=k;k=k+1;} 
      if( (b<eps | b > 1-eps)& (k<=N) ) 
                  {s[j]=b;b=pik[k];j=k;k=k+1;} 
}
u<-runif(1)
if(a>=eps & a<= 1-eps & b>=eps & b<= 1-eps)
if(a+b>1)
          {
            if(u<(1-b)/(2-a-b)) {b<-a+b-1;a<-1} 
            else {a<-a+b-1;b<-1} 
          }
 else{ if(u< b/(a+b)) 
                 		{b<-a+b;a<-0} 
       else {a<- a+b;b<-0} 
       } 
s[i]=a; s[j]=b;
s
}

