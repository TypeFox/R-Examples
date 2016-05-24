kpmfe <-
function(Vec,h,ker,a=1,c=2){
###########################################################################################################
# INPUTS:
#   "Vec" 		: Sample of data 
#   "h" 		: Bandwidth.
#   "ker" 		: The kernel function:  "dirDU" DiracDU,"bino" Binomial,"triang" discrete Triangular. 		   
#   "a" 		: The arm is used only for the Discrete Triangular kernel. The default value is 1.
#   "c" 		: The number of categories in the Aitchison and Aitken kernel  is used only for DiracaDU.The default value is 2.
# OUTPUT: Returns a list containing:
#    "n" 		: The number of observations.
#    "support" 		: The support of fn.
#    "C_n" 		: The normalizant constant.
#    "ISE_0" 		: The integrated squared error when using the naive distribution instead of fn.
#    "f_0" 		: The couples (x,f_0(x)).
#    "f_n" 		: The couples (x,f_n(x)).
#    "f0" 		: The empirical p.m.f.
#    "fn" 		: The estimated p.m.f. containing estimated values after normalization.
###########################################################################################################
  		 V=data.frame(table(Vec),row.names=NULL)
   		 N=V$Freq 
  		      if(ker=="dirDU"){x=0:(max(Vec))}
   		      else {x=0:(max(Vec)+2)}
 			 t1=rep(0,length(x))
 			 t2=rep(0,length(x))
   			 n <- length(x)
   			 f0=c(N/sum(N),rep(0,length(x)-length(N)))
   			 m=matrix(0,n,length(Vec))
    		     for(i in 1:n){
       			m[i,]= kf(x[i],Vec,h,ker,a,c)
               	   }
    		res<-apply(m,1,mean)

    		result<-res/sum(res)
   		E0=sum((result-f0)^2)
     		for (i in 1:n){
                 t1[i]=paste(x[i],";",f0[i])
		 t2[i]=paste(x[i],";",result[i])

           	 }
   		return(list(n=length(Vec),support= x,C_n=sum(res),ISE_0=E0,f_0=t1,f_n=t2,f0=f0,fn=result))

}
