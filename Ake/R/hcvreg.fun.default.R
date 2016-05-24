hcvreg.fun.default <-
function(Vec,y,type_data=c("discrete","continuous"),
			ker=c("bino","triang","dirDU","BE","GA","LN","RIG"),
			h=NULL,
			a0=0,
			a1=1,
			a=1,
			c=2,
			...)
{
    	 if(is.null(h)) 
	
            {
		h=seq(0.001,(max(Vec)-min(Vec))/2, length.out=1000)
	    }
	
    x<-Vec
    n <- length(x)
    m=matrix(0,n,length(Vec))
    m2=matrix(0,n,length(Vec))
    A=rep(0,length(h))
    for(k in 1:length(h)){
    for(i in 1:n){
         m[i,]= kef(x[i],Vec,h[k],type_data,ker,a0,a1,a,c)
         m2[i,]= m[i,]*y
	        }
	diag(m)<-0
	diag(m2)<-0
	G<-apply(m2,1, sum) 
        E<-apply(m,1, sum)
   	A[k]<-(1/length(Vec))*sum((y-(G/E))^2)
    }
    index<-which.min(A)
   hcv<-h[index]
 structure(list(kernel = ker,hcv=hcv,CV=A,seq_bws=h),class="hcvreg.fun")
}
