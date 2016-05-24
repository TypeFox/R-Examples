##' Automatically convert a vector of strings into a color for easy plotting
##'
##' @title Convert between strings to colors
##' @param string a vector of strings representing groups. 
##' @param colors a vector of colors, one for each unique element in \code{string}.
##' @export
##' @return a vector of colors, one for each element in \code{string} 
##' @author Dustin Fife
##' @examples
##' groups = sample(LETTERS[1:5], size=100, replace=TRUE)
##' plot(rnorm(100), rnorm(100), col=string.to.color(groups))
##' plot(rnorm(100), rnorm(100), col=string.to.color(groups), 
##'    pch=as.numeric(string.to.color(groups, colors=c(16:20))))
##' @note This function can also be used to specify pch values, cex values, or any other plotting values
##' the user may wish to differ across groups. See examples. 
string.to.color = function(string, colors=NULL){
	if (!is.null(colors)){
		if (length(colors)!=length(unique(string))){
			break("The number of colors must be equal to the number of unique elements.")
		}
		else {
			conv = cbind(unique(string), colors)
		}
	} else {
		conv = cbind(unique(string), rainbow(length(unique(string))))
	}
	unlist(lapply(string, FUN=function(x){conv[which(conv[,1]==x),2]}))
}
