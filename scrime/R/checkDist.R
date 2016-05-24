`checkDist` <-
function(x,dist){
	if(!is.matrix(dist) & !is.character(dist))
		stop("dist must be either a matrix or a character string.",call.=FALSE)
	if(is.matrix(dist))
		checkDistMat(dist,nrow(x),rownames(x))
	else{
		if(!dist%in%c("smc","cohen","pcc"))
			stop("dist must be either smc, cohen or pcc.",call.=FALSE)
		FUN<-match.fun(dist)
		dist<-FUN(x,dist=TRUE)
	}
	dist
}

