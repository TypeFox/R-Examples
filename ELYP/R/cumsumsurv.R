cumsumsurv <- function(x){
 	if ( any(is.na(x)) ) stop('NaNs');
 	s=x;
 	.C('cumsumsurv', x=as.numeric(x), s=as.numeric(s), LLL=length(x))$s
 	}



