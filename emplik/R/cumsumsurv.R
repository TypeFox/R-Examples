cumsumsurv <- function(x){
 	if(any(is.na(x))) stop('NaNs');    ## if (sum(is.na(x))>0) stop('NaNs');  3/2015 MZ
 	s=x;
 	.C('cumsumsurv',x=as.numeric(x),s=as.numeric(s),LLL=length(x))$s
 	}



