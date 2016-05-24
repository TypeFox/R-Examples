mad <- function(x, red = TRUE){
	lvls <- levels(as.factor(x))
	res <- data.frame(sapply(lvls, function(y) ifelse(x==y, 1, 0)))
	names(res) <- lvls
	if(!red){res <- res[,-length(lvls)]}
	return(res)
}