hcvc.fun.default <-
function(Vec,bw=NULL,type_data,ker,a0=0,a1=1,...)

{

if(is.null(bw)) {
	
		bw=seq((max(Vec)-min(Vec))/200,(max(Vec)-min(Vec))/2, length.out=100)
	    }


        result1<-0
	result2<-0
	x=seq(min(Vec),max(Vec),length.out=100)
	n1 <- length(x)
	n2 <- length(Vec)
	m1=matrix(0,n1,length(Vec))
	m2=matrix(0,n2,n2)
	Dak<-Vectorize(kef,vectorize.args=c('x','t'))
      for(k in 1:length(bw)){
              m1 <-outer(x,Vec,Dak,bw[k],"continuous",ker,a0,a1)
	      m2 <-outer(Vec,Vec,Dak,bw[k],"continuous",ker,a0,a1) 
		
              res1<-apply(m1,1,mean)
	      diag(m2)<-0
	      res2<-apply(m2,1,sum)
              result1[k]=simp_int(x,res1^2)
	      result2[k]=(2/((n2-1)*n2))*sum(res2)
    }
  	CV=result1-result2    
	index<-which.min(CV)  #to compute the optimal bandwidth
	hcv<-bw[index]

structure(list(hcv=hcv,seq_h=bw,CV=CV),class="hcvc.fun")

}
