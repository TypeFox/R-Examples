#' Pattern and Level Reliability via Profile Analysis
#' 
#' The \code{pr} function uses subscores from two parallel test forms and computes profile reliability coefficients as described in Bulut (2013).
#' 
#' Profile pattern and level reliability coefficients are based on the profile analysis approach described in Davison and Davenport (2002) and Bulut (2013). Using the parallel test forms or multiple administration of the same test form, pattern and level reliability coefficients are computed.  Pattern reliability is an indicator of variability between the subscores of an examinee and the level reliability is an indicator of  the average subscore variation among all examinees. For details, see \href{http://conservancy.umn.edu/bitstream/handle/11299/155592/Bulut_umn_0130E_13879.pdf}{Bulut(2013)}
#' 
#' 
#' @export
#' @importFrom stats cov
#' @param form1,form2 Two data matrices or data frames; rows represent individuals, columns represent subscores. Both forms should have the same individuals and subscores in the same order. Missing subscores have to be inserted as NA.
#' @return An object of class prof is returned, listing the following components:
#' \itemize{
#' \item \code{reliability} - Within-in person, between-person, and overall subscore reliability
#' \item \code{pattern.level} - A matrix of all pattern and level values obtained from the subscores
#' }
#' 
#' @examples 
#' \dontrun{
#' data(EEGS)
#' result <- pr(EEGS[,c(1,3,5)],EEGS[,c(2,4,6)])
#' print(result)
#' plot(result)}
#' 
#' @references Bulut, O. (2013). \emph{Between-person and within-person subscore reliability: Comparison of unidimensional and multidimensional IRT models}. (Doctoral dissertation). University of Minnesota. University of Minnesota, Minneapolis, MN. (AAT 3589000).
#' @references Davison, M. L., & Davenport, E. C. (2002). Identifying criterion-related patterns of predictor scores using multiple regression. \emph{Psychological Methods, 7}(4), 468-484. DOI: 10.1037/1082-989X.7.4.468
#' @seealso \code{\link{plot.prof}}
#' @keywords methods

pr <- function(form1,form2) {
	
	n <- ncol(form1)
	k <- nrow(form1)
	
	#Average score for each person (level)
	f1 <- as.matrix(rowMeans(form1,na.rm = TRUE),ncol=1,nrow=k)
	f2 <- as.matrix(rowMeans(form2,na.rm = TRUE),ncol=1,nrow=k)
	
	pattern1 <- matrix(,ncol=n, nrow=k)
	pattern2 <- matrix(,ncol=n, nrow=k)
	
	#Creating pattern scores
	for (i in 1:n) {
		pattern1[,i] <- form1[,i]-f1
		pattern2[,i] <- form2[,i]-f2
	}
	
	#Overall profile reliability
	covar1 <- matrix(,ncol=n, nrow=1)
	
	for (i in 1:n) {
		covar1[,i]=cov(form1[,i],form2[,i])
	}
	num1=rowSums(covar1,na.rm = TRUE)
	
	variance.form1 <- sum(apply(form1,2,var),na.rm = TRUE)
	variance.form2 <- sum(apply(form2,2,var),na.rm = TRUE)
	
	denum1 <- sqrt(variance.form1*variance.form2)
	
	overall <- num1/denum1
	
	
	#Level reliability
	num2 <- n*(cov(f1,f2))
	denum2 <- sqrt((n*var(f1))*(n*var(f2)))
	level <- num2/denum2
	
	
	#Pattern reliability
	covar2 <- matrix(,ncol=n, nrow=1)
	
	for (i in 1:n) {
		covar2[,i]=cov(pattern1[,i],pattern2[,i])
	}
	num3=rowSums(covar2,na.rm = TRUE)
	
	var.form1 <- sum(apply(pattern1,2,var),na.rm = TRUE)
	var.form2 <- sum(apply(pattern2,2,var),na.rm = TRUE)
	
	denum3 <- sqrt(var.form1*var.form2)
	
	pattern <- num3/denum3
	
	result1 <- rbind(level,pattern,overall)
	rownames(result1) <- c("Level","Pattern","Overall")
	colnames(result1) <- c("Estimate")
	
	colnames(pattern1) <- c(paste("Form1_pattern",c(1:n),sep=""))
	colnames(pattern2) <- c(paste("Form2_pattern",c(1:n),sep=""))
	colnames(f1) <- c("Level1")	
	colnames(f2) <- c("Level2")	
	
	result2 <- cbind(pattern1,f1,pattern2,f2)
        
        call<- match.call()
        
	output <- list(call=call,reliability=result1,pattern.level=result2)
	class(output) <- "prof"
	output
}





