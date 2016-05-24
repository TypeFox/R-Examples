###############################################################################
###							DEFINITION OF THE STRUCTURE OF REPGRID CLASS							  ###
###############################################################################

# In this file the repgrid classes are defined.
# Design note: 	the objects will be defined using S4 classes.

### NOTE:
# in the current approach (e.g. function implementation, use and definition, 
# especially in the arguments) that you might find restricting or tedious as 
# an advanced R user programmer. The main user group will be newcomers to R 
# though. Thus, implementation tries to provide maximal ease of use without
# requiring deeper R knowledge. Advanced R users will easily find their way 
# around those limitations.


# Definition of  repgrid class
#
# @slot meta           A list to store meta data for the repertory grid. 
#                      This includes name of interviewer and interviewee,
#                      data, miscellaneous notes etc.
# @slot scale          The rating scale used (minimum, maximum ec.). 
# @slot elements       The elements of the grid including meta information like
#                      "ideal" etc.
# @slot constructs     The constructs of the grid, containing meta information 
#                      like pole preference or different ladders.
# @slot elicitation    Information about the elicitation procedure used.
# @slot ratings        The ratings.
# @slot coupled        If the grid is coupled (standard) or decoupled (sci:vesco)
#                      format, allowing bent constructs.
# @slot calcs          Results from calculations.
# @slot plotdata       Information for plotting the grid.
#
# @export
# @author  Mark Heckmann
#
setClass( "repgrid", 
		  representation( meta = "list",
						          scale = "list",
						          coupled = "logical",
						          elements = "list",
						          constructs = "list",
						          elicitation = "list",
						          ratings = "array",
						          calcs = "list",
						          plotdata = "data.frame"))
						

#' Constructor for repgrid class
#' 
#' @return \code{repgrid} object       
#' @export
#' @keywords internal
#' @author  Mark Heckmann
#'
makeEmptyRepgrid <- function(){
	x <- new("repgrid")
	x@ratings <- array(NA, c(0, 0, 3)) 				# ,,1 = coupled ratings; decoupled ratings: ,,2 left pole  ,,3 right pole
	dimnames(x@ratings) <- list(constructs=NULL, elements=NULL,			# set up layers for coupled and decoupled rating
	  	layer=c("coupled", "left pole decoupled", "right pole decoupled"))			
	x
}


# #' Show method for testClass
# #' @param testClass object
# setMethod("show", "repgrid", function(object){
#   cat("object of class 'repgrid'")
# })

# #' Show method for repgrid
# #' @param repgrid object
# setMethod("show", signature= "repgrid", function(object){
#   x <- object 
#   showMeta(x)
#   showScale(x)    #print scale info
# })








