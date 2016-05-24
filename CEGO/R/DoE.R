###################################################################################
#' Initialize Random Design
#' 
#' Create a random initial population or experimental design, given a specifed creation function,
#' as well as a optional set of user-specified design members and a maximum design size.
#' Also removes duplicates from the design/population.
#'
#' @param x Optional list of user specified solutions to be added to the design/population, defaults to NULL
#' @param cf Creation function, creates random new individuals
#' @param size size of the design/population
#'
#' @return Returns list with experimental design (or population) without duplicates
#'
#' @keywords internal
#' @export
###################################################################################
designRandom <- function(x=NULL,cF,size,control=list()){
	## initialization
	if(is.null(x)){
    x <- list()
    k=0
  }else{ #given start population
    k=length(x)
  }
		
	if(k>size){
		x <- x[1:size]
	}else if(k<size){
		## CREATE initial population
		x <- c(x, replicate(size-k , cF(),simplify=FALSE))
	}#else if k==size do nothing.
		
	## REPLACE duplicates from initial population with unique individuals
	x <- removeDuplicates(x, cF)
}
