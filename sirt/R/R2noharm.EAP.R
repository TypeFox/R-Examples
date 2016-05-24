R2noharm.EAP <-
function( noharmobj , theta.k =seq(-6,6,len=21) , print.output=TRUE ){
    mod <- noharmobj
    mod$guess <- guess <- mod$guesses
    f0 <- mod$final.constants
    FL <- mod$loadings.theta
    D <- mod$dimensions
    P <- mod$factor.cor
    I <- length(mod$guess)
    dat <- mod$dat
    N <- nrow(dat)
    theta.k <- as.matrix( expand.grid( as.data.frame( matrix( rep( theta.k , D) , ncol = D ) ) ) )
    TT <- nrow(theta.k)		
	if ( TT > 1000 ){
		theta.k <- qmc.nodes( snodes = 1000 , ndim=D )
						}
    TT <- nrow(theta.k) 
    #***
    # calculate probabilities
    probs <- matrix( f0 , nrow= I  , ncol=TT , byrow=F)
    for (dd in 1:D){
        # dd <- 1
        probs <- probs + outer( FL[,dd] , theta.k[,dd] )
                }
    guessM <- matrix( mod$guess , nrow=I , ncol=TT )
    upperM <- matrix( mod$upper , nrow=I , ncol=TT )	
    probs1 <- guessM + ( upperM - guessM) * stats::pnorm( probs )
    probs <- array( 0 , dim=c(I , TT , 2 ) )
    probs[,,1] <- 1-probs1
    probs[,,2] <- probs1
    
    #****
    # evaluate posterior distribution
    prior.density <- mvtnorm::dmvnorm( as.matrix(theta.k) , mean = rep(0,D) ,sigma = as.matrix(P) )
    prior.density <- prior.density / sum( prior.density )
    prior.density <- matrix( prior.density , nrow=N , ncol=TT , byrow=T )
    like <- matrix( 1 , nrow=N , ncol=TT )
    for (ii in 1:I){
        # ii <- 1
        ind.ii <- which( ! is.na( dat[,ii] ) )
        like[ind.ii,] <- like[ ind.ii , ] * t( probs[ ii , , dat[ ind.ii ,ii ] + 1 ] )
            }
    posterior <- like * prior.density
    # posterior <- like
    posterior <- posterior / rowSums( posterior )
        
    #***
    # calculate EAP of all dimensions and reliabilities
    person <- data.frame( "case" = 1:N)
    nstudl <- rep(1,N)
    EAP.rel <- rep(0,D)
    hwtE <- posterior
    pweights <- nstudl
        for ( dd in 1:D ){
        # dd <- 1  # dimension
        person$EAP <- rowSums( hwtE * outer( nstudl , theta.k[,dd] ) )
        person$SD.EAP <- sqrt(rowSums( hwtE * outer( nstudl , theta.k[,dd]^2 ) ) - person$EAP^2)    
        #***
        # calculate EAP reliability
        # EAP variance
        EAP.variance <- stats::weighted.mean( person$EAP^2 , pweights ) - 
							( stats::weighted.mean( person$EAP , pweights ) )^2
        EAP.error <- stats::weighted.mean( person$SD.EAP^2 , pweights )
        EAP.rel[dd] <- EAP.variance / ( EAP.variance + EAP.error )    
        colnames(person)[ which( colnames(person) == "EAP" ) ] <- paste("EAP.Dim" , dd , sep="")
        colnames(person)[ which( colnames(person) == "SD.EAP" ) ] <- paste("SD.EAP.Dim" , dd , sep="")                
        }
    if ( print.output ){
		cat("EAP Reliabilities:\n")
		print( round (EAP.rel,3) )
					}
    res <- list( "person" = person , "theta"=theta.k , 
			"posterior" = posterior , "like"= like , 
            "EAP.rel" = EAP.rel , "probs" = probs )
    return(res)
    }
