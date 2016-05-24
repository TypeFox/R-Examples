#' Summarize matrices
#'
#' The function summarize matrices areas of the located spot matrices.
#' @param object S3 object of working TLC
#'
#' @return Returns S3 object with new values \code{object$spot_sums}.
#'
#' @author Ivan D. Pavicevic, \email{ivanp84@@gmail.com}
#' 
#' @examples
#' # This interactive example shows the most
#' # common usage of the qtlc library.
#' fname01 <- system.file("extdata", "test025to100sp.tiff", package="qtlc")
#' testTLC <- createTLC(fname01, RGB=FALSE)
#' print(testTLC)
#' 
#' # now using mouse select the spots with testTLC <- spot2D(testTLC)
#' # but, for automatic tests, we'll imitate that step...
#' testTLC$spots$x <- c(40.93354, 83.18687, 121.59899, 160.01111, 203.54485,
#'                      239.39616, 280.36909, 320.06161, 362.31494, 399.44666,
#'                      439.13919, 480.11211, 518.52423, 559.49716, 599.18969)
#' testTLC$spots$y <- c(198.3160, 198.3160, 199.2833, 198.3160, 198.3160,
#'                      198.3160, 198.3160, 198.3160, 197.3487, 198.3160,
#'                      199.2833, 198.3160, 199.2833, 199.2833, 199.2833)
#' 
#' # and now the select2D selects 30x30 pixels areas around spots
#' testTLC <- select2D(testTLC, 30, 30)
#' 
#' # forming spots matrices
#' testTLC <- matrices2D(testTLC)
#' 
#' # and finaly sumarizing spots areas
#' testTLC <- summat2D(testTLC)
#'
#' #eventually we'll examine the linear model
#' C <- rep(c(0.25, 1, 6.25, 25, 100), each=3) #imaginative concentrations
#' #now creates data frame with values
#' testTLC.df <- data.frame(C, testTLC$spot_sums)
#' names(testTLC.df) <- c("Concentration", "Signal")
#' # now the linear model
#' testTLC.lm <- with(testTLC.df, lm(Signal ~ Concentration))
#' # and finaly the plot
#' plot(testTLC.df)
#' abline(testTLC.lm)
#' summary(testTLC.lm)
#' 
#' @export
#' 
summat2D <- function(object) {
	sums <- c();
	mat <- object$spot_matrices;
	for(i in 1:dim(mat)[3]) {
		sums <- c(sums, sum(mat[,,i]));
		}
	object$spot_sums <- sums;
	return(object);
	}	

