#Utility function for building a correlation matrix 
#from a list of off-diagonal elements

buildCor<-function(n, cors) {
	rv <- diag(1,n)
		
	if(!missing(cors)) {
		if(is.vector(cors)) cors<-matrix(cors, byrow=T, nrow=1)	

		if(any( abs(cors[,3])>1 ) ) 
			stop(sprintf( "Correlations outside +-1 in %s",
				deparse(substitute(cors)) ), call.=TRUE)

		for(Row in 1:nrow(cors) ) {
			rv[ cors[Row,1], cors[Row,2] ] <-
				rv[ cors[Row,2], cors[Row,1] ] <-
					cors[Row,3]
		}
	}
	
	if(det(rv)<=0) warning("Returned value is not a valid correlation matrix", call.=TRUE)
	
	return(rv)
}

updateCor <- function(cor, cors, cor.names) {
	
	if(!missing(cors)) {
		if(is.vector(cors)) cors<-matrix(cors, byrow=T, nrow=1)	

		for(Row in 1:nrow(cors) ) {
			cor[ cors[Row,1], cors[Row,2] ] <-
				cor[ cors[Row,2], cors[Row,1] ] <-
					cors[Row,3]
		}
	}
	
	if(!missing(cor.names)){
		if(length(cor.names) == nrow(cor)) {
			dimnames(cor) <- list(cor.names,cor.names)
		} else {
			stop("cor.names length not compatible with correlation matrix size", call.=TRUE)
		}
	}
	
	if(det(cor)<=0) warning("Returned value is not a valid correlation matrix", call.=TRUE)
	
	return(cor)
}

buildCov <- function(s, covs, vars=s^2, cors, cov.names){
	if(missing(s) && missing(vars)) 
		stop("At least one of s and vars must be given", call.=TRUE)
		
	if(missing(s)) s <- sqrt(vars)
	
	if( !is.null(names(s)) && missing(cov.names) ) cov.names <- names(s)
	
	rv <- diag(vars)
	
	if(!missing(cors)){
		if(is.vector(cors)) cors<-matrix(cors, byrow=T, nrow=1)	

		if(any( abs(cors[,3])>1 ) ) 
			stop(sprintf( "Correlations outside +-1 in %s",
				deparse(substitute(cors)) ), call.=TRUE)
		
		for(Row in 1:nrow(cors) ) {
			rv[ cors[Row,1], cors[Row,2] ] <-
				rv[ cors[Row,2], cors[Row,1] ] <-
					cors[Row,3]*sqrt(rv[ cors[Row,1], cors[Row,1] ] *
							rv[ cors[Row,2], cors[Row,2] ])
		}
	}

	if(!missing(covs)){
		if(is.vector(covs)) covs<-matrix(covs, byrow=T, nrow=1)	

		for(Row in 1:nrow(covs) ) {
			rv[ covs[Row,1], covs[Row,2] ] <-
				rv[ covs[Row,2], covs[Row,1] ] <-
					covs[Row,3]
		}
	}

	if(!missing(cov.names)){
		if(length(cov.names) == nrow(cor)) {
			dimnames(cor) <- list(cov.names,cov.names)
		} else {
			stop("cor.names length not compatible with correlation matrix size", call.=TRUE)
		}
	
	}

	if(det(rv)<=0) warning("Returned value is not a valid correlation matrix", call.=TRUE)

	return(rv)
}

updateCov <- function(cov, covs, cors, cov.names) {
	for(Row in 1:nrow(covs) ) {
		cov[ covs[Row,1], covs[Row,2] ] <-
			cov[ covs[Row,2], covs[Row,1] ] <-
				covs[Row,3]
	}
	
	if(!missing(cors)){
		if(is.vector(cors)) cors<-matrix(cors, byrow=T, nrow=1)	

		if(any( abs(cors[,3])>1 ) ) 
			stop(sprintf( "Correlations outside +-1 in %s",
				deparse(substitute(cors)) ), call.=TRUE)
		
		for(Row in 1:nrow(cors) ) {
			cov[ cors[Row,1], cors[Row,2] ] <-
				rv[ cors[Row,2], cors[Row,1] ] <-
					cors[Row,3]*sqrt(cov[ cors[Row,1], cors[Row,1] ] *
							cov[ cors[Row,2], cors[Row,2] ])
		}
	}

	if(!missing(cov.names)){
		if(is.vector(covs)) covs<-matrix(covs, byrow=T, nrow=1)	

		if(length(cov.names) == nrow(cov)) {
			dimnames(cov) <- list(cov.names,cov.names)
		} else {
			stop("cov.names length not compatible with covariance matrix size", call.=TRUE)
		}
	}
	
	if(det(cov)<=0) warning("Returned value is not a valid covariance matrix", call.=TRUE)
	
	return(cov)
}
