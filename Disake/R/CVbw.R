CVbw <-
function(Vec,seq_bws=NULL,ker,a=1,c=2){
###########################################################################################################
# INPUTS:
#   "Vec" 	:The sample of data 
#   "seq_bws"   :The sequence of bandwidths where 
#		     the CV function will be evaluated (and minimized)
#   "ker" 	: The kernel:  "dirDU" DiracDU,"bino" Binomial,"triang" Discrete Triangular.
#   "a" 	:The arm of the Dicrete Triangular kernel. Default value is 1.
#   "c" 	:The  number of categories in the dirDU kernel. Default value is 2.     
#    
# OUTPUT:Returns a list containing:
#   "hcv" 	: The bandwidth obtained by cross-validation
#   "CV"        : The cross validation values 
#   "seq_bws"	: The sequence of the bandwiths used
###########################################################################################################

if(is.null(seq_bws)) 
	if (ker=="bino")
             {
	        seq_bws=seq((max(Vec)-min(Vec))/500,1, length.out=100) 
	     }
	else
            {
		seq_bws=seq((max(Vec)-min(Vec))/200,(max(Vec)-min(Vec))/2, length.out=100)
	    }


        result1<-0
	result2<-0
	if(ker=="dirDU"){x=0:(max(Vec))}# The values on the support must be up to max
  	else {x=0:(max(Vec)+2)}# and up to two points after the max for other kernels.
	n1 <- length(x)
	n2 <- length(Vec)
	m1=matrix(0,n1,length(Vec))
	m2=matrix(0,n2,n2)
        Dak<-Vectorize(kf,vectorize.args=c('x','t')) 
      for(k in 1:length(seq_bws)){
              m1 <-outer(x,Vec,Dak,seq_bws[k],ker,a,c)
	      m2 <-outer(Vec,Vec,Dak,seq_bws[k],ker,a,c) 
		
              res1<-apply(m1,1,mean)
	      diag(m2)<-0
	      res2<-apply(m2,1,sum)
              result1[k]=sum(res1^2)
	      result2[k]=(2/((n2-1)*n2))*sum(res2)
    }
  	CV=result1-result2    
	index<-which.min(CV)  #to compute the optimal bandwidth
	hcv<-seq_bws[index]

	 return(list(hcv=hcv,CV=CV,seq_bws=seq_bws))

}
