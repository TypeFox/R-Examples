CVbw <-
function(type_kernel="n",vec_data,n_pts=100,
seq_bws=NULL)
# INPUTS:
#   "type_kernel" kernel function: "e" Epanechnikov,	"n" Normal, "b" Biweight
#   "vec_data" data sample
#   "n_pts" the number of points to be used in the integration
#   "seq_bws" sequence of bandwidths where the CV function will be evaluated (and minimized)
 
# OUTPUT:
# Returns the CV bandwidth according to the criterion by Bowman, Hall and Prvan (1998)  
{
    prob_quantile <- 0
    ss <- quantile(vec_data, c(prob_quantile,1-prob_quantile))
    y <- seq(ss[1],ss[2],length.out=n_pts)
	if(is.null(seq_bws)) 
	seq_bws=seq((max(vec_data)-min(vec_data))/200,(max(vec_data)-min(vec_data))/2, length.out=50)
    n_bws <- length(seq_bws)
    CVfunction <- numeric(length=n_bws)
    for (i in 1:n_bws)
		{
      integrand <- apply(((outer(y,vec_data,"-")>=0)-t(kernel_distribution_without_i(type_kernel,y,vec_data,seq_bws[i])))^2,1,mean)
      CVfunction[i] <- simp_int(y,integrand)$value
    }
    i0 <- which.min(CVfunction)
    CVbw_val <- seq_bws[i0]
    return(list(seq_bws=seq_bws, CVfunction=CVfunction, bw=CVbw_val))
}
