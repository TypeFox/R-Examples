
##############################################################
# function lambda part
align.optim.lambda <- function( lambda , psi0 , psi0b ,
            align.scale , align.pow , wgt , eps=.0001 ,
			group.combis){
    # optimization with respect to country SDs
    lambda1 <- lambda / psi0
    lambda1b <- lambda / psi0b    
    # function
    fopt <- 0
	I <- ncol(lambda)
    for (ii in 1:I){
        # ii <- 1
        fopt1 <- ( lambda1[ group.combis[,1] , ii ] - lambda1b[ group.combis[,2] , ii ] )^2
        fopt <- fopt + wgt[ group.combis[,1],ii] * wgt[ group.combis[,2],ii] * 
                        ( fopt1 / align.scale^2 + eps )^align.pow
                }
    res <- rowsum( fopt , group.combis[,1] )
	return(res[,1] )
        }

##############################################################
# function lambda part
align.optim.nu <- function( lambda , nu , psi0 , psi0b , 
			alpha0 , alpha0b ,
            align.scale , align.pow , wgt , eps=.0001 ,
			group.combis){
    # optimization with respect to country SDs
#    nu1 <- nu - alpha0 * lambda / psi0
#    nu1b <- nu - alpha0b * lambda / psi0
    nu1 <- nu - alpha0 * lambda 
    nu1b <- nu - alpha0b * lambda
    # function
    fopt <- 0
	I <- ncol(lambda)
    for (ii in 1:I){
        # ii <- 1
        fopt1 <- ( nu1[ group.combis[,1] , ii ] - nu1b[ group.combis[,2] , ii ] )^2
        fopt <- fopt + wgt[ group.combis[,1],ii] * wgt[ group.combis[,2],ii] * 
                        ( fopt1 / align.scale^2 + eps )^align.pow
                }
    res <- rowsum( fopt , group.combis[,1] )
	return(res[,1] )
        }
#################################################


#################################################
# alignment Newton Raphson step
align.newton.raphson <- function( ll0 , ll1 , ll2 , max.increment , h ){
    d1 <- ( ll1 - ll2  ) / ( 2 * h )  
    # second order derivative
    # f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
    d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
    # change in item difficulty
    d2[ abs(d2) < 10^(-10) ] <- 10^(-10)
    increment <- - d1 / d2
    increment <- ifelse( abs( increment) > max.increment , 
            max.increment*sign(increment) , increment ) 
    return(increment)
        }
###############################################


# auxiliary function for calculation of correlations
ai.calc.corr <- function(parsM){
	# parsM <- t(lambda.aligned)
	cM <- stats::cor( parsM)
	I <- ncol(cM)
	rbar <- ( sum(cM) - I )/ ( I^2 - I)
	return(rbar)
		}
