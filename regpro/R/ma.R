ma<-function(x,h=1,kernel="exp",k=length(x))
{
if (kernel=="exp") 
   ker<-function(xx){ return( exp(-xx) ) }
if (kernel=="bart") 
   ker<-function(xx){ return( 1-xx^2 ) }
if (kernel=="gauss") 
   ker<-function(xx){ return( exp(-xx^2/2) ) }
if (kernel=="uniform") 
   ker<-function(xx){ ans<-(xx^2 <= 1) 
                      return( ans ) }

t<-seq(0,k-1,1)
w<-ker(t/h)
w<-w/sum(w)
w<-w[length(w):1]
ans<-sum(w*x)

return(ans)
}


