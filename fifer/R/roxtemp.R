##' Generate a Roxygen Template
##'
##' @title Generate a Roxygen Template
##' @author Dustin Fife
##' @export
roxtemp = function(){
	f = paste0(
"##' Insert title here
##'
##' Write the description here
##'	
##' Write the details here
##' @param
##' @param
##' @aliases
##' @seealso \\code{\\link{MASS}}
##' @references
##' @return \\item{}{}
##' @author
##' @export
##' @examples"
)
	cat(f)
}