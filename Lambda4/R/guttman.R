#' Guttman's 6 Lambda Coefficients
#' 
#' @description Calculates all 6 of Guttman's lambda coefficients.
#' 
#' @return
#' \item{Lambda1}{Guttman's Lambda1 estimate of reliability.}
#' \item{Lambda2}{Guttman's Lambda2 estimate of reliability.}
#' \item{Lambda3}{Guttman's Lambda3 estimate of reliability. Also known as Cronbach's alpha or coefficient alpha.}
#' \item{Lambda4}{Guttman's maximimal Lambda4 estimate of reliability.}
#' \item{Lambda5}{Guttman's Lambda5 estimate of reliability.}
#' \item{Lambda6}{Guttman's Lambda6 estimate of reliability.}
#' 
#' @param x Can be either a data matrix or a covariance matrix
#' @param missing How to handle missing values.
#' @param standardize When TRUE Results are standardized by using the correlation matrix instead of the covariance matrix for computation.
#' 
#' @note The estimate for Lambda4 is maximized.
#' 
#' @references
#' Guttman L (1945). "A Basis for Analyzing Test-Retest Reliability." Psychometrika, 10, 255-282.
#' @author Tyler Hunt \email{tyler@@psychoanalytix.com}
#' @examples
#' guttman(Rosenberg)
#'
#' @export

guttman <- function(x, missing="complete", standardize=FALSE){
	
	sigma <- impute.cov(x, missing)
	
	if(standardize==TRUE){
		sigma <- cov2cor(sigma)
	}

  Lambda1 <- as.numeric(lambda1(sigma))
  Lambda2 <- as.numeric(lambda2(sigma))
  Lambda3 <- as.numeric(lambda3(sigma)$lambda3[1])
  Lambda4 <- as.numeric(cov.lambda4(sigma, show.splits=FALSE)$lambda4[2])
  Lambda5 <- as.numeric(lambda5(sigma))
  Lambda6 <- as.numeric(lambda6(sigma))


  result <- list(Lambda1=Lambda1, 
                 Lambda2=Lambda2, 
                 Lambda3=Lambda3, 
                 Lambda4=Lambda4, 
                 Lambda5=Lambda5, 
                 Lambda6=Lambda6)
  
  class(result) <- c("guttman")
  return(result)
}

