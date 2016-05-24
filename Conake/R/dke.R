dke <-
function(vec_data,ker,bw,x=NULL,a=0,b=1)
# INPUTS:
#   "ker" kernel function: "GA" Gamma,	"BE" extended beta, 
#                          "LN" lognormal, RIG "reciprocal inverse Gaussian"         
#   "vec_data" sample of data 
#   "bw" bandwidth
#   "x" single value or grid where the kernel estimation is computed
#   "a"  left bound of the support
#   "b"  right bound of the support
# OUTPUT:Returns a list containing:
#    "C_n" the normalizing constant
#    "f_n" vector containing the estimated function
#    in the grid values
{
		n <- length(vec_data)
		if(is.null(x)) {
		x=seq(min(vec_data),max(vec_data),length.out=100)}
		aux <- matrix(data=vec_data,nrow=length(vec_data),ncol=length(vec_data),byrow=TRUE)
		#for(i in 1:n){
       		#	aux[i,]= kef(x[i],vec_data,bw,ker,a,b)
                # 	   }
	 	 aux <- kef(x,aux,bw,ker,a,b)
		res<- apply(aux,1,mean)	 # density without normalization
     		C<-simp_int(x,res)	 # Normalizant constant 
		result<-res/C 		 # density normalized
		 result<-res/C
		
return(list(C_n=C,f_n=result))
}
