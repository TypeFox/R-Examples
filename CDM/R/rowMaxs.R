################################################################################
# utility method for computing intermediate information                        #
################################################################################
# rowMaxs <- function(mat){
# Call: from din()
# Input: numeric matrix
# Output: row maxima of input matrix
#    n <- nrow(mat)
#    p <- ncol(mat)
#    x <- as.vector(mat)
#    x <- matrix(x[order(rep(1:n, p), x)], p, n)
#    x[p , ]
#}
#########################################################
rowMaxs <- function(mat){
    n <- nrow(mat)
    p <- ncol(mat)
	maxval <- mat[,1]
	for ( cc in 2:p){
		maxval <- ifelse( mat[,cc] > maxval  , mat[,cc] , maxval )
				   }
	return(maxval)
}
###########################################################
rowMaxs2 <- function(mat){
    n <- nrow(mat)
    p <- ncol(mat)
	maxval <- mat[,1]
	maxind <- 1
	for ( cc in 2:p){
		ind <- ( mat[,cc] > maxval )
		maxval <- ifelse( ind , mat[,cc] , maxval )
		maxind <- ifelse( ind , cc , maxind)
				   }
	res <- list( "maxval" = maxval , "maxind" = maxind )
	return(res)
}
#############################################################
rowMaxs3 <- function(mat){
    n <- nrow(mat)
    p <- ncol(mat)
	maxval <- mat[,1]
	maxind <- 1
	for ( cc in 2:p){
		maxval <- ifelse( mat[,cc] > maxval , mat[,cc] , maxval )
		maxind <- maxind + ( cc - maxind ) * (mat[,cc] ==maxval )		
				   }
	res <- list( "maxval" = maxval , "maxind" = maxind )
	return(res)
}
# rowMaxs3 is faster than rowMaxs2!

