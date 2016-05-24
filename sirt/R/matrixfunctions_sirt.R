
##########################################################################
# rowwise maximum and minimum function
#****************************
# rowMaxs function
rowMaxs.sirt <- function(matr){ 
	.Call("rowMaxsCPP_source", matr , PACKAGE = "sirt")
					}
#*****************************					
# rowMins function
rowMins.sirt <- function(matr){
	matr2 <- - matr
    res2 <- rowMaxs.sirt( matr2 )
    res <- list( "minval" = - res2$maxval ,
					"minind" = res2$maxind )
	return(res)
		}
##########################################################################
# rowwise cumulative sum
rowCumsums.sirt <- function(matr){ 
	.Call("rowCumsums2_source", matr , PACKAGE = "sirt")
					}
# The C code was posted by Romain Francois at
# http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2010-October/001198.html

##########################################################################
# rowwise cumulative sum
colCumsums.sirt <- function(matr){ 
	t( rowCumsums.sirt( t(matr) ) )
					}

##########################################################################
#****					
# 'interval_index' searches an index when a frequency is exceeded
# -> used in plausible value imputation
rowIntervalIndex.sirt <- function(matr,rn){ 
	.Call("interval_index_C", matr , rn , PACKAGE = "sirt")
					}	

##########################################################################
# extract k smallest elements in a row of a matrix
rowKSmallest.sirt <- function( matr , K , break.ties=TRUE){
    M1 <- matr
    N1 <- dim(M1)[1] ; N2 <- dim(M1)[2]
    # generate random number matrix
    rM1 <- matrix( round( stats::runif( N1*N2 ) ) , N1 , N2 )
    if ( ! break.ties ){ rM1 <- 0*rM1 }
    # define integer matrix
    indexmatr <- matrix( 1:N2 , N1 , N2 , byrow=TRUE )
    # apply function for extracting k smallest elements
	a1 <- .Call("rowKSmallest_C", matr , K , indexmatr , 
					rnmatr=rM1 , PACKAGE = "sirt")
	## OUTPUT:					
	## return List::create(_["smallval"]=SMALLVAL ,
	##		              _["smallind"]=SMALLIND ) ;  					
    return(a1)
        }
##########################################################################
# Ksmallest -> different implementation
rowKSmallest2.sirt <- function(matr , K ){
    Nmis <- nrow(matr)
    disty <- matr
    donors <- K
    indvec <- 1:Nmis
    M1 <- max(disty)+1
    smallval <- donor.ind <- matrix( 0 , nrow=Nmis , ncol=donors )
    res1 <- rowMins.sirt(matr=disty)
    donor.ind[,1] <- res1$minind
    smallval[,1] <- res1$minval    
    for (ii in 2:donors){
        disty[ cbind(indvec , donor.ind[,ii-1] ) ] <- M1
        res1 <- rowMins.sirt(matr=disty)
        donor.ind[,ii] <- res1$minind
        smallval[,ii] <- res1$minval    
                        }          
    res <- list( "smallval"=smallval , "smallind" = donor.ind )
    return(res)
    }
