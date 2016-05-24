rl <-
function(type_kernel="n",vec_data,
T, bw=PBbw(type_kernel="n", vec_data, 2))
# INPUTS:
#   "type_kernel" kernel function: "e" Epanechnikov,	"n" Normal, 
#                                  "b" Biweight, "t" Triweight         
#   "vec_data" data sample
#   T a sequence of times (years, days... )
#   "bw" bandwidth
# OUTPUT:Returns a list containing:
#    "Estimated_values" vector containing the estimated function
#    in the values of T
#    "grid" T
{
		
	s0<-min(vec_data)
	s1<-max(vec_data)
	p<-1-1/T
	nw<-length(p)
	rl<-0
	for (i in 1:nw)
	{
	rl[i]<-dichotomy_fun(type_kernel="n", vec_data,s0,s1,bw,p[i])
	}
	return(rl)
}
