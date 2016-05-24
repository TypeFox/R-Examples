NULL

#'
#'  DEPRECATED Converts several samples \code{x} random variable (daily precipitation values) extracted by populations represented by the columns of \code{data} respectively or \code{sample}
#'  to a normally-distributed samples with assinged mean and standard deviation or vice versa in case \code{inverse} is \code{TRUE} using the function \code{\link{normalizeGaussian_prec}}
#'  
#' 
  
#' @param x value to be converted 
#' @param data a sample of data on which a non-parametric probability distribution is estimated
#' @param cpf cumulative probability distribution. If \code{NULL} (default) is calculated as \code{\link{ecdf}(data)}
#' @param mean mean (expected value) of the normalized random variable. Default is 0.
#' @param sd standard deviation of the normalized random variable. Default is 1.
#' @param inverse  logical value. If \code{TRUE} the function works inversely (the opposite way). Default is \code{FALSE}.
#' @param qnull probability of no precipitation occurence. (It can be a matrix in case \code{sample="monthly"}
#' @param valmin minimum value of precipitation to consider a wet day


# @param step vector of values in which step discontinuities of the cumulative probability function occur. Default is \code{NULL}
# @param prec amplitude of the neighbourhood of the step discontinuities where cumulative probability function is treated as non-continuous.
#' @param type see \code{\link{quantile}}
#' 
#' @param extremes logical variable. 
#'  If \code{TRUE} (default) the probability or frequency is multiplied by \deqn{\frac{N}{N+1}} where \eqn{N} is the length of \code{data}
#' @param sample information about sample or probability distribution. Default is \code{NULL}
#' @param origin_x date corresponding to the first row of \code{x}
#' @param origin_data date corresponding to the first row of \code{data}

#' 
#' @export 
#' @author Emanuele Cordano, Emanuele Eccel
#' @return a matrix or a data.frame with the normalized variable or its inverse   
#' 
#' @seealso   \code{\link{normalizeGaussian_prec}}
#'   
#'     
#' @note In the version 1.2.5 of \pkg{RMAWGEN} This function is deprecated and not used.
# It applies \code{\link{normalizeGaussian}} for each column of \code{x} and \code{data}.
# See the R code for further details
#' 




normalizeGaussian_severalstations_prec <-
function(x,
		data=x,
		cpf=NULL,mean=0,
		sd=1,
		inverse=FALSE,
		qnull=NULL,
		valmin=0.5,
	#	step=NULL,
	#	prec=10^-4,
		type=3,
		extremes=TRUE,
		sample=NULL,
		origin_x=NULL,
		origin_data=NULL
		


) {
	
	

	
	out <- x*NA
	

	
	
	if (is.null(sample)) {
		
		
		
		for (i in 1:ncol(x)) {
			out[,i] <- normalizeGaussian_prec(x=x[,i],data=data[,i],cpf=cpf,mean=mean,sd=sd,inverse=inverse,type=type,extremes=extremes,sample=sample,valmin=valmin,qnull=qnull) 
			
		}
		

	} else if (sample=="monthly") {
		
		months <- months((0.5:11.5)*365/12,abbreviate=TRUE)
		
		for (m in 1:length(months)) {
			
			i_months_x <- extractmonths(data=1:nrow(x),when=months[m],origin=origin_x)
		
			i_months_data <- extractmonths(data=1:nrow(data),when=months[m],origin=origin_data)
		
			for (i in 1:ncol(x)) {
				
				out[i_months_x,i] <- normalizeGaussian_prec(x=x[i_months_x,i],data=data[i_months_data,i],cpf=cpf,mean=mean,sd=sd,inverse=inverse,valmin=valmin,qnull=qnull,type=type,extremes=extremes,sample=NULL) 	
					
			}
			

			
			
		
		
			
		}	
		
		
		
		
	
	
	} else if (sample=="monthly_year"){
	 
	# TO DO 	year <-
		# TO DO 
		
	} else {
		
		print("Error in normalizaGaussian_sevaralStation_prec: sample option not yet implemented!!")
	
	}
	

	names(out) <- names(x)

	return(out)
	
}

