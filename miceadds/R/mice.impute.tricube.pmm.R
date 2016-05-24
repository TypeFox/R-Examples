mice.impute.tricube.pmm <- function (y, ry, x, tricube.pmm.scale= .2 , tricube.boot = FALSE , ...){
    x <- cbind(1, as.matrix(x))
	# print some informations
        vname <- get("vname", pos = parent.frame()) # get variable name        
    cat( "\n" , paste( vname , "  Imputation Method tricube.pmm with scaling factor" , tricube.pmm.scale , "\n"))
    parm <- .norm.draw(y, ry, x, ...)
    yhatobs <- x[ry, ] %*% parm$coef
    yhatmis <- x[!ry, ] %*% parm$beta
	ymean <- mean( y[ry] )
	#***
	# bootstrap
	if ( tricube.boot ){
		y1 <- y[ry] ; x1 <- x[ry,]
		B <- length(y1)
		ind <- sample( 1:B , replace = TRUE )
		parm2 <- .norm.draw(y = y1[ind] , ry = rep(TRUE,B) , x = x1[ind,] , ...)	
		yhatmis <- x[!ry, ] %*% parm2$beta		
				}
	#***
	# R2 calculations
		R2.orig <- 1 - sum( ( y[ry] - yhatobs )^2 ) / sum( ( y[ry] - ymean)^2 ) 
		R2.samp <- 1 - sum( ( x[ry, ] %*% parm$beta - y[ry] )^2 ) / sum( ( y[ry] - ymean)^2 ) 
		cat( paste( "  R2 (original data): " , round(R2.orig,4)  , "\n"))
		cat( paste( "  R2 (sampled coefficients): " , round(R2.samp,4)  , "\n"))	
		if ( tricube.boot ){ 
				R2.boot <- 1 - sum( ( x[ry, ] %*% parm2$beta - y[ry] )^2 ) / sum( ( y[ry] - ymean)^2 ) 		
				cat( paste( "  R2 (Bootstrap): " , round(R2.boot,4)  , "\n"))	
					}
	#***
    # extract scale parameter for tricube pmm
    vname <- get("vname", pos = parent.frame()) 
    tricube.pmm.scale <- .extract.list.arguments( micearg = tricube.pmm.scale , 
                           vname = vname , miceargdefault = .2 )
    flush.console()
    # doing tricube pmm
    x1 <- apply(as.array(yhatmis), 1, .tricube.pmm.match , yhat = yhatobs, y = y[ry], 
                            tricube.pmm.scale = tricube.pmm.scale , ... )
    return(x1)
    }

	
	



#---------------------------------------------------##
# tricube predictive mean matching                  ##
#   weighted according Tukey's tricube function     ##
.tricube.pmm.match <- function (z, yhat = yhat, y = y, donors = 3, tricube.pmm.scale = .2 , ...) 
{
    d <- abs(yhat - z)
    donorset <- which( rank(d, ties.method = "ran") <= donors )
    s.tricube <- tricube.pmm.scale * mean( d )
    prob.x <- unlist( sapply( d , FUN = function(dd){ ( 1- min( dd / s.tricube , 1 )^3  )^3 } ) )
    # prevent the case that all weights are equal to zero
    prob.x[ donorset ] <- prob.x[donorset] + .0001
    # standardize weights to probabilities
    prob.x <- prob.x / sum(prob.x )
    m <- sample( y , size = 1 , prob = prob.x )
    return(m)
}
#-----------------------------------------------------

