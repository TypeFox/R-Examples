`checkDistMat` <-
function(dist,n.row,rn){
	if(nrow(dist)!=ncol(dist))
		stop("dist must be a squared matrix.",call.=FALSE)
	if(nrow(dist)!=n.row)
		stop("x and dist must have the same number of rows.",call.=FALSE)
	tmp<-t(dist)
	if(any(dist!=tmp))
		stop("dist must be a symmetric matrix.",call.=FALSE)
	tmp<-rownames(dist)
	if(!is.null(tmp) && !is.null(rn) && any(tmp!=rn))
		stop("x and dist must have the same row names.",call.=FALSE)
}

