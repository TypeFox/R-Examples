printme <- function(x, SEP=" "){
    if(length(x)==1)return(maybeQuote(x))
    if(length(x) > 1){
		X <- maybeQuote(x)
		res <- paste(paste(X, SEP, sep=""), collapse="")
		return(res)
		}
}

