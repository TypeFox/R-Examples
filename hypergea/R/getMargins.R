getMargins <-
function(x){
	if( !is.array(x) ){ warning(paste0(sQuote("x"), " must be matrix or array")) }
return(lapply(1:length(dim(x)), function(v){ apply(x,v,sum) }))
}
