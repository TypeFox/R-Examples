
#########################################
# data preparation exchanegable raters
lc2.data.prep <- function(data){
	maxK <- max( data )
	m1 <- matrix( NA , maxK+1 , maxK+1 )
	rownames(m1) <- colnames(m1) <- paste0("Cat" , 0:maxK)
	for (kk1 in 0:maxK){
	for (kk2 in 0:maxK){
		m1[kk1+1,kk2+1] <- sum( ( data[,1] == kk1 ) * ( data[,2] == kk2 ) )
						}
				}
	m1 <- m1+t(m1)
	m1 <- m1/2
	res <- list("m1"=m1 , "maxK" = maxK)
	return(res)
	}
###############################################

###################################################
# agreement measures: input is a frequency table
lc2.agreement <- function( m1 ){
	# agreement and kappa measures
	#***
	# agreement
	res <- list( "agree0" = sum( diag( m1 ) ) / sum( m1 ) )
	K <- ncol(m1)
	qk <- outer( 1:K , 1:K  , "-")
	res$agree1 <- sum( m1 * ( abs(qk) <=1 ) ) / sum(m1)
	p1 <- m1 / sum(m1)
	#***
	# kappa
	pk <- colSums( m1 )
	pk <- pk/sum(pk)
	pek <- sum( pk^2 )
	res$kappa <- ( sum( diag(p1) ) - pek ) / ( 1 - pek )

	#***
	# weighted kappa

	# define weights
	wgt <- 1 - abs( qk ) / ( K - 1 )
	pa <- sum( wgt * p1 )
	pek <- sum( outer( pk , pk )*wgt )
	res$wtd.kappa.linear <- ( pa - pek) / ( 1 - pek )

	wgt <- 1 - qk^2 / ( K - 1 )^2
	pa <- sum( wgt * p1 )
	pek <- sum( outer( pk , pk )*wgt )
	res$wtd.kappa.quadratic <- ( pa - pek) / ( 1 - pek )

	#***
	# Gwet's AC(1) statistic
	pek <- sum( pk*(1-pk) ) / ( K - 1 )
	pa <- sum( diag(p1))
	res$AC1 <- ( pa - pek) / ( 1 - pek )
	res$alpha.aickin <- alpha.aickin(m1)	
	res <- unlist(res)
	}
################################################################

######################################################
# Aickin's alpha
# Aickin (1990) Biometrics
alpha.aickin <- function(m1){
    # Aickin's alpha statistic
    p1 <- m1 / sum(m1)
    pk <- colSums(p1)
    pek <- sum( pk^2 )
    pa <- sum( diag(p1) )
    ckappa <- ( pa - pek ) / ( 1 - pek )
    Q <- ncol(m1)    
    alpha <- ckappa
    pkH <- pk
    pet <- sum( pkH^2 )
    conv <- .001
    iter <- 0
    maxit <- 100
    par.change <- 1    
    while( ( par.change > conv ) & ( iter < maxit ) ){
        alpha0 <- alpha
        pkH0 <- pkH
        # update hard-to-classify probabilities
        pkH <- pk / ( ( 1 - alpha ) + alpha * pkH / pet )
        # update alpha
        pet <- sum( pkH^2 )
        alpha <- ( pa - pet ) / ( 1 - pet )
        par.change <- max( abs( c( alpha - alpha0 , pkH - pkH0 ) ) )
        iter <- iter + 1
#        print( par.change )
                }
    return(alpha)
        }
############################################################