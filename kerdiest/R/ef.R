ef <-
function(type_kernel="n",vec_data,
c, bw=PBbw(type_kernel="n", vec_data, 2), Dmin=0, Dmax=15,
size_grid=50, lambda)
# INPUTS:
#   "type_kernel" kernel function: "e" Epanechnikov,	"n" Normal, 
#                                  "b" Biweight, "t" Triweight         
#   "vec_data" data sample
#   "c"  an specific value on which to calculate the exceedance function
#   "Dmin" min value for D time units  (years, days... )
#   "Dmax" max value for D time units  (years, days... )
#    "size_grid" length of a grid where to compute the exceedance function
#    "landa" mean activity rate
#   "bw" bandwidth
# OUTPUT:Returns a list containing:
#    "Estimated_values" vector containing the estimated function
#    in the grid values
#    "grid" the used grid 
#    "bw" value of the  bandwidth
{
		
		D<-seq(Dmin, Dmax,length.out=size_grid)
		F_c<- kde(type_kernel, vec_data=vec_data,
		y=c,bw=bw)$Estimated_values
		exc<-1-exp(-lambda*D*(1-F_c))
  	return(list(Estimated_values = exc, grid=D, bw = bw))
}
