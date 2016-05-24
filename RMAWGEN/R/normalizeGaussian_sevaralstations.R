NULL

#'
#'  Converts several samples \code{x} random variable  extracted by populations represented by the columns of \code{data} respectively or \code{sample}
#'  to a normally-distributed samples with assinged mean and standard deviation or vice versa in case \code{inverse} is \code{TRUE}
#'    
#' @param x value to be converted 
#' @param data a sample of data on which a non-parametric probability distribution is estimated
#' @param cpf cumulative probability distribution. If \code{NULL} (default) is calculated as \code{\link{ecdf}(data)}
#' @param mean mean (expected value) of the normalized random variable. Default is 0.
#' @param sd standard deviation of the normalized random variable. Default is 1.
#' @param inverse  logical value. If \code{TRUE} the function works inversely (the opposite way). Default is \code{FALSE}.
#' @param step vector of values in which step discontinuities of the cumulative probability function occur. Default is \code{NULL}
#' @param prec amplitude of the neighbourhood of the step discontinuities where cumulative probability function is treated as non-continuous.
#' @param type see \code{\link{quantile}}
#' @param extremes logical variable. 
#'  If \code{TRUE} (default) the probability or frequency is multiplied by \deqn{\frac{N}{N+1}} where \eqn{N} is the length of \code{data}
#' @param sample information on how to sample \code{x} and \code{data}. Default is \code{NULL}, this means that the values of each column of \code{x} and \code{data} belong to the same sample. If \code{x} and \code{data} are sampled for each month seperately, it is set to \code{monthly}.
#' @param origin_x date corresponding to the first row of \code{x}
#' @param origin_data date corresponding to the first row of \code{data}

#' 
#' @export 
#' @author Emanuele Cordano, Emanuele Eccel
#' @return a matrix with the normalized variable or its inverse   
#' 
#' @seealso   \code{\link{normalizeGaussian}}
#'   
#'     
#' @note It applies \code{\link{normalizeGaussian}} for each column of \code{x} and \code{data}.
#' See the R code for further details

#'
#' @examples 
#' 
#' library(RMAWGEN) 
#' N <- 30
#' x <- rexp(N)
#' y <- x+rnorm(N)
#' df <- data.frame(x=x,y=y)
#' 
#' dfg <- normalizeGaussian_severalstations(df,data=df,extremes=TRUE,inverse=FALSE)
#' 
#' dfi <- normalizeGaussian_severalstations(dfg,data=df,extremes=TRUE,inverse=TRUE)
#' 
#' 
#'  
#' 




normalizeGaussian_severalstations <-
function(x,
		data=x,
		cpf=NULL,mean=0,
		sd=1,
		inverse=FALSE,
		step=NULL,
		prec=10^-4,
		type=3,
		extremes=TRUE,
		sample=NULL,
		origin_x=NULL,
		origin_data=NULL
		


) {
	
	

	
	out <- x*NA
	

	
	
	if (is.null(sample)) {
		
		
		
		for (i in 1:ncol(x)) {
			out[,i] <- normalizeGaussian(x=x[,i],data=data[,i],cpf=cpf,mean=mean,sd=sd,inverse=inverse,step=step,prec=prec,type=type,extremes=extremes,sample=sample) 
			
		}
		

	} else if (sample=="monthly") {
		
		months <- months((0.5:11.5)*365/12,abbreviate=TRUE)
		
		for (m in 1:length(months)) {
			
			i_months_x <- extractmonths(data=1:nrow(x),when=months[m],origin=origin_x)
		
			i_months_data <- extractmonths(data=1:nrow(data),when=months[m],origin=origin_data)
		
			for (i in 1:ncol(x)) {
				
				out[i_months_x,i] <- normalizeGaussian(x=x[i_months_x,i],data=data[i_months_data,i],cpf=cpf,mean=mean,sd=sd,inverse=inverse,step=step,prec=prec,type=type,extremes=extremes,sample=NULL) 	
					
			}
			

			
			
		
		
			
		}	
		
		
		
		
	
	
	} else if (sample=="monthly_year"){
	 
	# TO DO 	year <-
		# TO DO 
		
	} else {
		
		print("Error in normalizaGaussian_sevaralStation: sample option not yet implemented!!")
	
	}
	

	names(out) <- names(x)

	return(out)
	
}

