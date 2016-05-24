###################################################################
# Function for defining different response patterns
resp.pattern3 <-
function( x ){
    n <- nrow(x)
    p <- ncol(x)
    mdp <- (x %*% (2^((1:ncol(x)) - 1))) + 1
    misspattern <- mdp[,1]
    misspattern <- list( "miss.pattern" = mdp[,1] , 
                "mp.index" = match( mdp[,1] , sort( unique(mdp[,1] ) ) ) )
    return( misspattern )
        }
####################################################################
rowcumsums <-
  function(m1){
    g1 <- 0*m1
    g1[,1] <- m1[,1]
    for (ss in seq(2,ncol(m1))){
      g1[,ss] <- g1[,ss-1] + m1[,ss] 
    }
    return(g1)
  }

###################################################################
rowCumsums.TAM <- function(matr){ 
	.Call("rowCumsums2_source", matr , PACKAGE = "TAM")
					}
###################################################################					
#****					
# 'interval_index' searches an index when a frequency is exceeded
# -> used in plausible value imputation
interval_index <- function(matr,rn){ 
	res <- .Call("interval_index_C", matr , rn , PACKAGE = "TAM")
	res <- res + 1
	return(res)
					}					
#############################################################
# search the maximum in each matrix row
rowMaxs <-
  function(mat, na.rm = FALSE){    
    # Call: from designMatrix()
    # Input: 
    # mat: numeric matrix
    # na.rm: logical. Should missing values (including NaN) be omitted from the calculations?
    # Output: row maxima of input matrix    
    n <- nrow(mat)
    p <- ncol(mat)
    x <- as.vector(mat)
    x <- matrix(x[order(rep(1:n, p), x, na.last = !na.rm)], p, n)
    x[p , ]
  }

#############################################################
# rewrite theta.sq function into Rcpp
theta.sq <-
  function(theta){
    theta2 <- array(,dim = c(nrow(theta), ncol(theta) , ncol(theta) ) )
    for( qq in 1:nrow(theta) ){
#		theta2[qq,,] <- theta[qq,] %*% t(theta[qq,])  
		theta2[qq,,] <- tcrossprod( theta[qq,] )  		
				}
    return("theta2" = theta2)
  }
#*******************************
# faster Rcpp function  
theta.sq2 <- function(theta){
	theta2 <- .Call("theta_sq_cpp", theta , PACKAGE = "TAM")
    return("theta2" = theta2)
  }  
  
#############################################################

add.lead <- function(x, width=max(nchar(x))){
  sprintf(paste('%0', width, 'i', sep=''), x) 
}

##############################################################