#orthbasis <- function(D){
#		V <- matrix(0,nrow=D,ncol=D-1)
#		for (i in 1:ncol(V)){
#			V[1:i,i] <- 1/i
#			V[i+1,i] <- (-1)
#			V[,i] <- V[,i]*sqrt(i/(i+1))
#		}
#		return(V)
#}



#' Orthonormal basis
#' 
#' Orthonormal basis from cenLR transformed data to isomLR transformated data.
#' 
#' For the chosen balances for \dQuote{isomLR}, this is the orthonormal basis
#' that transfers the data from centered logratio to isometric logratio.
#' 
#' @param D number of parts (variables)
#' @return the orthonormal basis.
#' @author Karel Hron, Matthias Templ
#' @seealso \code{\link{isomLR}}, \code{\link{cenLR}}
#' @keywords manip
#' @export
#' @examples
#' 
#' data(expenditures)
#' V <- orthbasis(ncol(expenditures))
#' xcen <- cenLR(expenditures)$x.clr
#' xi <- as.matrix(xcen) %*% V
#' xi2 <- isomLR(expenditures)
#' all.equal(xi, xi2)
#' 
orthbasis <- function(D){
	ilrBase <- NULL
	gsicomp <- 
			function (W = c(1, -1)) 
	{
		## this function is a copy of function gsi.buildilr* from 
		## package compositions. The authors are 
		## Raimon Tolosana-Delgado, K.Gerald v.d. Boogaart
		## and the intellecutal properties belongs to them.
		if (length(W) < 2) {
			return(ilrBase(D = 1))
		}
		if (length(dim(W)) == 0) {
			return(ilrBase(D = 2))
		}
		if (length(dim(W)) > 0) {
			W = as.matrix(W)
			nc = ncol(W)
			D = nrow(W)
			isPos = (W > 0)
			isNeg = (W < 0)
			nPos = matrix(1, D, D) %*% isPos
			nNeg = matrix(1, D, D) %*% isNeg
			W = (isPos * nNeg - isNeg * nPos)
			nn = sapply(1:nc, function(i) {
						1/sqrt(W[,i] %*% W[,i]) 
					})
			nn = matrix(nn, ncol = ncol(W), nrow = nrow(W), byrow = TRUE)
			W = W * nn
			return(W)
		}
	}
	codes=matrix(rep(0,(D-1)*D),ncol=D)
	for(i in 1:(D-1)){
		for(j in 1:D){
			codes[i,]=c(rep(0,i-1),-1,rep(1,D-i))
		}
	}   
	t(codes)	
	V <- gsicomp(t(codes))
	return(V)
}
