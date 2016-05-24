computeMW <-
function(DATA,masses=NULL,weights=NULL){
	diagmats <- TRUE
	DATA_dimensions <- dim(DATA)
	
	if((DATA_dimensions[1] > 1000) || (DATA_dimensions[2] > 1000)){
		#use the vector approach to masses & weights.
		diagmats <- FALSE
	}

	if(!is.null(masses)){
		if(length(masses)==dim(DATA)[1]){
			#if(diagmats){
			#	M <- diag(masses)
			#}else{
				M <- masses
			#}
			
		}else if(dim(masses)[1]==dim(masses)[2]){
			#if(diagmats){
			#	M <- masses
			#	masses <- diag(M)
			#}else{
				M <- diag(masses)
				masses <- M
			#}			
		}else{
			masses <- c(matrix(1/nrow(DATA),1,nrow(DATA)))
			#if(diagmats){
			#	M <- diag(masses)	
			#}else{
				M <- masses
			#}
		}
	}else{
		masses <- c(matrix(1/nrow(DATA),1,nrow(DATA)))
		#if(diagmats){
		#	M <- diag(masses)	
		#}else{
			M <- masses
		#}
	}
	
	if(!is.null(weights)){
		if(length(weights)==dim(DATA)[2]){	
			#if(diagmats){
			#	W <- diag(weights)
			#}else{
				W <- weights
			#}
		}else{
			weights <- c(matrix(1/ncol(DATA),1,ncol(DATA)))
			#if(diagmats){
			#	W <- diag(weights)
			#}else{
				W <- weights
			#}
		}
	}else{
		weights <- c(matrix(1/ncol(DATA),1,ncol(DATA)))
		#if(diagmats){
		#	W <- diag(weights)
		#}else{
			W <- weights
		#}
	}
	return(list(M=M,W=W))
}
