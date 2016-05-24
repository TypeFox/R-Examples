revise.tpm <- function(xi,mixture) {
	if(mixture)  matrix(apply(xi,2,sum)/sum(xi),byrow=TRUE,
                            nrow=nrow(xi),ncol=ncol(xi))
	else xi/apply(xi,1,sum)
}
