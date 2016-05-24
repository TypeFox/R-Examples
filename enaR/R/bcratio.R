#' bcratio --- calculates the ratio of positive to 
#' negative elements in a flow model
#' INPUT = flow matrix
#' OUTPUT = ratio of positive to negetive elements
#' 
#' S. Borrett | July 2011
#' Singh, Hines, Borrett | Update June 2014
#' ---------------------------------------------------
bcratio <- function(x='matrix'){
	r=sum(x[x>0])/sum(abs(x[x<0]))  # creates the ratio of positive elements to negative elements.
  return(r)
}
