##' A function to create useless hashes (####) for ease of commenting
##'
##' @title Create useless hashes
##' @param j How many columns of hashes
##' @param i How many rows of hashes
##' @author Dustin Fife
##' @export
##' @examples
##' hash(j=100, i=11)
hash = function(j=50, i=4){
	k = paste(rep("#", times=j),collapse="")
	for (m in 1:i){
		cat("\n")
		cat(k)
	}     
	cat("\n\n")   
	for (m in 1:i){

		cat(k)
				cat("\n")
	}     
}    
