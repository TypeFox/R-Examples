dke.fun.default <-
function(Vec,h,type_data=c("discrete","continuous"),ker=c("BE","GA","LN","RIG"),x=NULL,a0=0,a1=1,...)
              
# INPUTS:
#   "ker" kernel function: "GA" Gamma,	"BE" extended beta, 
#                          "LN" lognormal, RIG "reciprocal inverse Gaussian"         
#   "Vec" sample of data 
#   "h" bandwidth
#   "x" single value or grid where the kernel estimation is computed
#   "a0"  left bound of the support
#   "a1"  right bound of the support
# OUTPUT:Returns a list containing:
#    "C_n" the normalizing constant
#    "f_n" vector containing the estimated function
#    in the grid values
{
		n <- length(Vec)
		if(is.null(x)) {
		x=seq(min(Vec),max(Vec),length.out=100)}
		aux <- matrix(data=Vec,nrow=length(x),ncol=length(Vec),byrow=TRUE)
		#for(i in 1:n){
       		#	aux[i,]= kef(x[i],vec_data,bw,ker,a,b)
                # 	   }
	 	 aux <- kef(x,aux,h,"continuous",ker,a0,a1)
		res<- apply(aux,1,mean)	 # density without normalization
     		C<-simp_int(x,res)	 # Normalizant constant 
		result<-res/C 		 # density normalized
structure(list(data=Vec,n=length(Vec),hist=hist(Vec,prob=TRUE),eval.points= x,h=h, kernel=ker,C_n=C,est.fn=result),class="dke.fun") 
}
