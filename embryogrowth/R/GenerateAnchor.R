#' GenerateAnchor Generate a set of anchored parameters
#' @title Generate a set of anchored parameters
#' @author Marc Girondot
#' @return A vector with parameters
#' @param temperatures A vector with temperatures to serve as anchors
#' @param nests Formated nest data or result object obtained from searchR()
#' @param parameters A set of parameters value
#' @param number.anchors Number of anchors
#' @description Generate a set of anchored parameters.\cr
#' It is important that the anchors (i.e. the temperatures used as anchors) encompass 
#' the highest and lowest temperatures that are present in nests.
#' @examples
#' \dontrun{
#' # Example to generate anchored parameters
#' newp <- GenerateAnchor()
#' newp <- GenerateAnchor(temperatures=seq(from=20, 
#'   to=35, length.out=7))
#' newp <- GenerateAnchor(number.anchors=7)
#' data(nest)
#' formated <- FormatNests(nest, previous=NULL)
#' newp <- GenerateAnchor(nests=formated)
#' newp <- GenerateAnchor(nests=formated, number.anchors=10)
#' data(resultNest_4p)
#' newp <- GenerateAnchor(nests=resultNest_4p, number.anchors=7)
#' newp <- GenerateAnchor(nests=resultNest_4p, temperatures=seq(from=20,
#'  to=35, length.out=10))
#' newp <- GenerateAnchor(nests=resultNest_4p, number.anchors=7)
#' newp <- c(newp, Scale=1)
#' }
#' @export

GenerateAnchor <- function(temperatures=NULL, nests=NULL, parameters=NULL,
                      number.anchors=7) {

if (number.anchors<7 | (length(temperatures)<7 & !is.null(temperatures))) {
	print("At least 7 anchors must be defined")
	return()
}

if (is.null(temperatures)) {

	 if (!is.null(nests)) {
	
		t <- hist(nests, plot=FALSE)
		temperatures <- seq(from=floor(range(t$temperatures)[1]+273.15-1), 
		 to=floor(range(t$temperatures)[2]+273.15+1), length.out=number.anchors)
		 
	} else {
		 temperatures <- seq(from=20, to=35, length.out=number.anchors)+273.15		 
		}
} else {
  
  if (any(temperatures < 273.15)) temperatures <- temperatures+273.15
  
  
}

  newx <- rep(2, length(temperatures))
  
if (!is.null(parameters)) {
	newx <- .SSM(temperatures, parameters)[[1]]*1E5
} else 
  if (class(nests)=="NestsResult") 	newx <- .SSM(temperatures, nests$par)[[1]]*1E5


names(newx) <- temperatures

return(newx)
}

