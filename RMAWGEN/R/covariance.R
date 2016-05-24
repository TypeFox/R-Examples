NULL

#'
#'  Calculates the covariance matrix of the normally standardized variables obtained from the columns of \code{x} 

#' @param x variable
#' @param data a sample of data on which a non-parametric pghjjrobability distribution is estimated
#' @param cpf cumulative probability distribution. If \code{NULL} (default) is calculated as \code{\link{ecdf}(data)}
#' @param mean mean (expected value) of the normalized random variable. Default is 0.
#' @param sd standard deviation of the normalized random variable. Default is 1.
#' @param step vector of values in which step discontinuities of the cumulative probability function occur. Default is \code{NULL}
#' @param prec amplitude of the neighbourhood of the step discontinuities where cumulative probability function is treated as non continuous.
#' @param type see \code{\link{quantile}}
#' @param extremes logical variable. 
#'  If \code{TRUE} (default) the probability or frequency is multiplied by \deqn{\frac{N}{N+1}} where \eqn{N} is the length of \code{data}
#' @param sample information about sample or probability distribution. Default is \code{NULL}
#' @param origin_x date corresponding to the first row of \code{x}
#' @param origin_data date corresponding to the first row of \code{data}
# @param set.mcovcorrection logical value. If \code{TRUE} \code{\link{mcovcorrection}} is enabled. Warning, \code{step=0} is required.
# @param cov_min minimum admitted value of covariance. See \code{\link{mcovcorrection}}
# @param cov_max maximum admitted value of covariance. See \code{\link{mcovcorrection}}
#' 
#' @export 
#'  
#' @param use see \code{\link{cov}} 
#' 
#' @author Emanuele Cordano, Emanuele Eccel
#' @return a matrix with the normalized variable or its inverse   
#' 
#' @seealso   \code{\link{normalizeGaussian_severalstations}},\code{\link{normalizeGaussian}}
#'   
#'     
#'  @note It applies \code{\link{normalizeGaussian_severalstations}} to \code{x} and \code{data} and then calculates the covariances among the column.
#' See the R code for further details







covariance <-
function(x,
		data=x,
		cpf=NULL,
		mean=0,
		sd=1,
		step=NULL,
		prec=10^-4,
		use="pairwise.complete.obs",
		type=3,extremes=TRUE,
		sample=NULL,
		origin_x=NULL,
		origin_data=origin_x
# mcovcorrectionparameter
#		set.mcovcorrection=FALSE,		
	#	cov_min=0.001,
	#	cov_max=1.0,
	#	use.mcov=FALSE

) {

	
	
	
	out <- array(NA,c(ncol(x),ncol(x)))
	
	if (is.null(sample)) sample="NULL"
	
	if (sample=="gaussian") {
		y <- x
	} else if (sample=="gaussian_check"){
		
		x_inv <- normalizeGaussian_severalstations(x=x,data=data,cpf=cpf,mean=mean,sd=sd,step=step,prec=prec,type=type,extremes=extremes,sample="monthly",origin_data=origin_data,origin_x=origin_x,inverse=TRUE)
		y <- normalizeGaussian_severalstations(x=x_inv,data=data,cpf=cpf,mean=mean,sd=sd,step=step,prec=prec,type=type,extremes=extremes,sample="monthly",origin_data=origin_data,origin_x=origin_x)
		
	} else {
		if (sample=="NULL") sample=NULL
		y <- normalizeGaussian_severalstations(x=x,data=data,cpf=cpf,mean=mean,sd=sd,step=step,prec=prec,type=type,extremes=extremes,sample=sample,origin_data=origin_data,origin_x=origin_x)
	}
	
#	if (use.mcov) {
#		step=0
#		s_y <- array(NA,ncol(y)) 
#		for (j in 1:length(s_y)) {
#			
#			s_y[j] <- normalizeGaussian(x=step+prec/2,data=data[,j],cpf=cpf,mean=mean,sd=sd,type=type,extremes=extremes,sample=NULL) 	
#		}
	
#		out <- mcov(s_x=s_y,x=y,cov_min=cov_min,cov_max=cov_max)
		
		
	# TO DO	
		
#	} else {
		
	for (r in 1:nrow(out)){
			
		for (c in 1:ncol(out)) {
				
			out[r,c] <- cov(y[,r],y[,c],use=use)
				
		}
	}
		
#}
	
	
	
	
	return(out)
}

