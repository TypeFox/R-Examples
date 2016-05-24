

#############################################################
#############################################################
# Utility Functions                                         #
# Function calculates the rowwise product of a matrix       #
# Note that the entries must be nonnegative                 #
rowProds <- function(matr){
    # INPUT:                                
    # matrix with positive entries          
    exp( rowSums( log(matr + 10^(-300) ) ) )
    }
# for nonnegative entries use this function in combination      ##
# with sign(matr)                                               ##
# rowProds( matr ) * rowProds( sign(matr) )                     ##
##################################################################
# alternative to rowProds
rowProds2 <- function(matr){
	y <- matr[,1]
	for (ii in 2:dim(matr)[2]){
		y <- y * matr[,ii] }
	return(y)
		}
#...................................................................

# columnwise product
colProds <- function(matr){ 
        exp( colSums( log(matr + 10^(-300) ) ) )
        }

#---------------------------------------------------------------------
# Function of D. Rizoupoulos                                         #
rowMedians <- function(mat){
    stopifnot(is.matrix(mat), typeof(mat) == "double")
    if (any(is.na(mat)))
        stop("'mat' should not contain NAs")
    n <- nrow(mat)
    p <- ncol(mat)
    x <- as.vector(mat)
    x <- matrix(x[order(rep(1:n, p), x)], p, n)
    if (p%%2)
        x[0.5 * (p + 1), ]
    else
        0.5 * (x[0.5 * p, ] + x[0.5 * p + 1, ])
}
colMedians <- function(mat){ rowMedians( t(mat) ) }
#-----------------------------------------------------------------------


#-------------------------------------------------------------------
rowVars <- function(mat , na.rm= F ){
    n <- rowSums( 1 - is.na(mat) ) 
    ( rowSums( mat^2 , na.rm= T) - n * rowMeans( mat , na.rm = na.rm )^2 ) / ( n - 1 )
    }
#*****
colVars <- function( mat , na.rm=F){ rowVars( t(mat) , na.rm ) }
#*****
rowSds <- function( mat , na.rm=F){ sqrt(rowVars( mat , na.rm ) ) }
#*****
colSds <- function( mat , na.rm=F){ sqrt(colVars( mat , na.rm ) ) }
#-------------------------------------------------------------------
min.vec <- function(a,b){ifelse( a >= b , b , a ) }
#*********************************
rowMins2 <- function(matr){
#	y <- matr[,1]
#	for (ii in 2:dim(matr)[2]){
#		y <- min.vec( y , matr[,ii] ) }
	y <- do.call( pmin , as.data.frame(matr) )
	return(y)
		}
#------------------------------------------------------------------------#
		
#------------------------------------------------------------------------#
rowMaxs <- function(mat){
    n <- nrow(mat)
    p <- ncol(mat)
    x <- as.vector(mat)
    x <- matrix(x[order(rep(1:n, p), x)], p, n)
    x[p , ]
}
#------------------------------------------------------------------------#
colMaxs <- function(mat){ t( rowMaxs( mat ) ) }
#------------------------------------------------------------------------#

#------------------------------------------------------------------------#
whichrowMaxs <- function(mat){
    n <- nrow(mat)
    p <- ncol(mat)
    x <- as.vector(mat)
    dfr <- data.frame( x , rep(1:n, p) , rep( 1:p ,each= n ) )
    ind <- order(rep(1:n, p), x)
    arg <- matrix(  dfr[ ind , 3]  , p , n )[p,]
    x <- matrix(x[ind], p, n)   
    val <- x[p , ]
    list( "val" =val , "arg" = arg )
}
#------------------------------------------------------------------------#
whichcolMaxs <- function(mat){ t( whichrowMaxs( mat ) ) }
#------------------------------------------------------------------------#
whichrowMins <- function(mat){  whichrowMaxs( - mat ) }
#------------------------------------------------------------------------#
whichcolMins <- function(mat){  whichcolMaxs( - mat ) }
#------------------------------------------------------------------------#


#--------------------------------------------------------------
# rowwise cumsum operation on matrices                        #
rowCumsums <- function( mat , multmat = NULL){
# multmat <- NULL
    if (is.null(multmat)){
        m <- ncol(mat)
        multmat <- matrix( 1 , ncol  = m , nrow= m ) 
        multmat[ lower.tri( multmat )] <- 0
                }    
    mat %*% multmat
    }
#--------------------------------------------------------------#
colCumsums <- function( mat  ){  rowCumsums( t(mat) ) }
#--------------------------------------------------------------#


#-------------------------------------------------------------------------------
# compute column-bundlewise rowSums of a matrix
rowCumsums.colbundles <- function( mat , ind , multmat = NULL ){
    if (is.null(multmat)){  matlist <- as.list( rep(1,length(ind)))
        for (ii in 1:( length(ind) )){
            mat1 <- matrix( 1 , nrow=ind[ii] , ncol=ind[ii] )
            mat1[ lower.tri( mat1 ) ] <- 0
            matlist[[ii]] <- mat1
                }
        multmat <- Matrix::bdiag(lapply(matlist, as.matrix))
        }
    as.matrix( mat %*% multmat )
        }
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# compute column-bundlewise rowSums of a matrix
rowSums.colbundles <- function( mat , ind , multmat = NULL ){
    if (is.null(multmat)){
        matlist <- as.list( rep(1,length(ind)))
        for (ii in 1:( length(ind) )){
            mat1 <- matrix( 1 , nrow=ind[ii] , ncol=ind[ii] )
            matlist[[ii]] <- mat1
                }
        multmat <- Matrix::bdiag(lapply(matlist, as.matrix))
        }
    as.matrix( mat %*% multmat )
        }
#-------------------------------------------------------------------------------



#..........................................................................
# function for calculation of Cronbach's Alpha
.cronbach <- function(matr){
                matr <- stats::na.omit(matr)
                p <- ncol(matr)
                p / ( p - 1 )* ( sum( matr ) - sum( diag( matr ) ) )/ sum(matr )
                    }
#..........................................................................
