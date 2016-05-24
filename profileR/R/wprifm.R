#' Within-Person Random Intercept Factor Model
#' 
#' @importFrom lavaan sem
#' @export
#' @param data Data.frame containing the manifest variables.
#' @param scale Should the data be scaled? Default = FALSE
#' @param save_model Should the temporary lavaan model syntax be saved. Default = FALSE
#' @return an object of class \code{lavaan}
#'
#'@details This function performs the within-person random intercept factor model described in Davison, Kim, and Close (2009). For information about this model, please see this reference. This function returns an object of \code{lavaan} class and thus any generics defined for \code{lavaan} will work on this object. This function provides a simple wrapper for \code{lavaan}.
#'
#' @examples
#' data <- HolzingerSwineford1939[,7:ncol(HolzingerSwineford1939)]
#' wprifm(data, scale = TRUE)
#' @references Davison, M., Kim, S.-K., Close, C. (2009). Factor analytic modeling of within person variation in score profiles. Multivariate Behavioral Research, 44(5), 668 - 687. DOI: 10.1080/00273170903187665

wprifm <- function(data, scale = FALSE, save_model = FALSE){
	if(scale){
		data <- as.data.frame(scale(data))
	} 
	
    # Define random Intercept Factor
	f1x <- paste("1", paste("x", 1:ncol(data), sep=""), sep="*")
    f1 <- paste("f1 =~ ", paste(f1x[1:length(f1x)-1],"+ ", collapse=""), f1x[length(f1x)], sep = "") 

    # Define factor of interest
    f2params <- paste("lam",1:ncol(data),sep="")
    f2x <- paste(f2params, paste("x", 1:ncol(data), sep=""), sep="*")
 	f2 <- paste("f2 =~ NA*x1 + ", paste(f2x[1:length(f2x)-1],"+ ", collapse=""), f2x[length(f2x)], sep = "")

 	# Constraints
 	constraints <- "f1 ~~ f2\nf2 ~~ 1*f2"
 	cx <- paste("0 == ", paste(f2params[1:length(f2params)-1], " + ",sep = "", collapse = ""), f2params[length(f2params)], sep = "")

 	# Write to a temporary file
 	write(c(f1,f2,constraints,cx), file = "tmp.out")
 	model <- readLines("tmp.out")

 	# Remove file unless asked not to
 	if(save_model == FALSE){
 		system("rm tmp.out")
 	}
 	fit <- sem(model, data = data)
 	return(fit)
}


