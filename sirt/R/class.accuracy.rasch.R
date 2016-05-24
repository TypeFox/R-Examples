############################################################
# classification accuracy in the Rasch model
class.accuracy.rasch <-
function( cutscores , b , meantheta , sdtheta , theta.l , n.sims=0 , seed=988){
    dat0 <- matrix( 1 , ncol=length(b) , nrow=length(theta.l))
    set.seed(seed)
    cat("Cut Scores \n")
    print( cutscores , digits=3 )

	CC <- length(cutscores)+1 # number of classes
    L <- length(theta.l)    # number of theta points
    cutscores.long <- c(-Inf , cutscores , Inf )
    # MLE standard error
    semle <- rasch.info.mle(dat=dat0, theta=theta.l , b=b )
    wgttheta <- stats::dnorm( theta.l , mean = meantheta , sd = sdtheta )
    wgttheta <- wgttheta / sum( wgttheta )
    
    # distribution of estimated theta values
    m1 <- matrix( NA , nrow = L , ncol=L )
    for (tt in 1:L ){ 
        # tt <- 1
        w.tt <- stats::dnorm( theta.l , mean= theta.l[tt]  , sd = semle[tt] )
        w.tt <- w.tt / sum( w.tt )
        m1[tt,] <- w.tt
                    }
    # calculate probability classification
    m2 <- matrix( NA , CC , CC )
    classes <- paste("Class" , 1:CC , sep="")
    rownames(m2) <- paste("True" , classes ,sep="_" )
    colnames(m2) <- paste("Est" , classes ,sep="_")
    m0 <- matrix( wgttheta  , L , L , byrow=FALSE  )
    for (cc in 1:CC){
        # cc <- 1
        for (ee in 1:CC){
            # ee <- 1
            ind.cc <- which( ( theta.l >= cutscores.long[cc] ) & ( theta.l < cutscores.long[cc+1] ) )
            ind.ee <- which( ( theta.l >= cutscores.long[ee] ) & ( theta.l < cutscores.long[ee+1] ) )
            m2[cc,ee] <- sum( m0[ind.cc , ind.ee ] * m1[ ind.cc , ind.ee ] )
                    }
            }
    stats <- data.frame( "agree0"= sum( diag(m2)) )
    stats$agree1 <- sum( m2 * ( abs( outer( 1:CC  , 1:CC , "-" ) ) <= 1 ) )
    
    # calculate kappa calculation
    chance.prob <- rep(NA, CC)
    for (cc in 1:CC){ 
        ind.cc <- which( ( theta.l >= cutscores.long[cc] ) & ( theta.l < cutscores.long[cc+1] ) )
        chance.prob[cc] <- sum( m0[ ind.cc , 1 ] )
                }
    phi.c <- sum( chance.prob^2 )
    stats$kappa <- ( stats$agree0 - phi.c ) / ( 1 - phi.c )
	rownames(stats) <- "analytical"	
	#########################################
	# accuracy by simulation
	if ( n.sims > 0 ){
		theta0 <- stats::rnorm(n.sims, mean=meantheta , sd=sdtheta )
		theta.sim <- data.frame( "theta.true" = theta0 )
		theta.sim$true.class <- as.numeric( paste( cut( theta0 , breaks= cutscores.long , labels=1:CC) ))
		theta.sim$random.class <- sample( theta.sim$true.class )
		dat.sim <- sim.raschtype( theta=theta0 , b=b )
		dat.sim2 <- sim.raschtype( theta=theta0 , b=b )		
		wle1 <- wle.rasch( dat.sim , b = b )
		theta.sim$WLE <- wle1$theta
		theta.sim$WLE.class <- as.numeric( paste( cut( theta.sim$WLE , breaks= cutscores.long , labels=1:CC )) )
		wle2 <- wle.rasch( dat.sim2 , b = b )
		theta.sim$WLE2 <- wle2$theta
		theta.sim$WLE2.class <- as.numeric( paste( cut( theta.sim$WLE2 , breaks= cutscores.long , labels=1:CC )) )
		# calculate WLE reliability
		v2 <- mean( ( theta.sim$WLE - theta0)^2 )
		wle.rel <- stats::var( theta0) / ( var( theta0) + v2 )
		stats2 <- stats
		rownames(stats2) <- "simulated"
		stats2$agree0 <- mean( theta.sim$true.class == theta.sim$WLE.class )
		stats2$agree1 <- mean( abs( theta.sim$true.class - theta.sim$WLE.class ) <=1 )
		phi.c <- mean( theta.sim$true.class == theta.sim$random.class )
		stats2$kappa <- ( stats$agree0 - phi.c ) / ( 1 - phi.c )
		stats <- round( rbind( stats , stats2 ) , 3 )
		cat("\nWLE reliability (by simulation) =" , round( wle.rel , 3) , "\n" )
		stats[2,"consistency"] <- mean( theta.sim$WLE2.class == theta.sim$WLE.class )	
		cat("WLE consistency (correlation between two parallel forms) = " )
		cat( round( stats::cor( theta.sim$WLE , theta.sim$WLE2 ) , 3 ), "\n")
					}
	cat("\nClassification accuracy and consistency\n")
    print( stats , digits=3 )
    cat("\nProbability classification table \n")
    print( round(m2,3) , digits=3 )
    res <- list( "class.stats" = stats , "class.prob" = m2)
    invisible(res)
        }
