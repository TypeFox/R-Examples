			
######################################################
# compute deviance
.mcmc.deviance.3pno.testlet <- function( aM , bM , a.testletM , theta , 
		guess ,	gamma.testlet , testletgroups , dat , dat.resp , 
				weights , eps , param){
    guessM <- matrix( guess , nrow=nrow(aM) , ncol=length(guess) )
	gamma.testletM <- gamma.testlet[ , testletgroups ]
	if (param==1){ tau.ni <- aM * theta + gamma.testletM + bM }
	if (param==2){ tau.ni <- aM * theta + aM*gamma.testletM + bM }	
	if (param==3){ tau.ni <- aM * theta + a.testletM*gamma.testletM + bM }		
	pij <- guessM + ( 1 - guessM )* stats::pnorm( tau.ni )
	llij <- log( dat.resp * ( dat*pij + ( 1-dat )*(1-pij) ) + eps )
	if ( is.null( weights ) ){ deviance <- -2*sum( llij ) }
	if ( ! is.null( weights ) ){ 
		deviance <- -2*sum( rowSums(llij) * weights ) 
					}
	return(deviance)
				}
				
######################################################				
# subfunction for calculating the DIC
.mcmc.ic.3pno.testlet <- function( a.chain , b.chain , 
		a.testlet.chain , N , I , 
		theta.chain , c.chain , gamma.testlet.chain , TT , 
		testletgroups , dat , dat.resp , weights , eps , param ,
		deviance.chain ){
	#************
	aM <- matrix( colMeans( a.chain ) , nrow=N , ncol=I , byrow=TRUE )
	bM <- matrix( colMeans( b.chain ) , nrow=N , ncol=I , byrow=TRUE )
	if (param==3){	
			a.testletM <- matrix( colMeans( a.testlet.chain ) , nrow=N , ncol=I , byrow=TRUE )
				}
	theta <- colMeans( theta.chain )
	guess <- colMeans( c.chain )	
	gamma.testlet <- colMeans( gamma.testlet.chain )
	gamma.testlet <- matrix( gamma.testlet , nrow=N , ncol=TT+1 )
	Dhat <- .mcmc.deviance.3pno.testlet( aM , bM , a.testletM , theta , 
					guess ,	gamma.testlet , testletgroups , dat , dat.resp , 
					weights , eps , param )
	Dbar <- mean( deviance.chain )
	pD <- Dbar - Dhat
	ic <- list( "Dhat"=Dhat , "Dbar"=Dbar , "pD"=pD , "DIC" = Dhat + 2*pD )
	return(ic)
		}
#################################################################
# output mcmc.list
.mcmc.list.3pno.testlet <- function( a.chain , b.chain , a.testlet.chain , 
	I , deviance.chain , est.slope , c.chain , sigma.chain ,
		est.guess , theta.chain , sigma.testlet.chain , TT ,
		burnin , SV , save.theta , testletgroups , param ){
	#***********************************
	a <- a.chain
	b <- b.chain
	theta <- theta.chain
	N <- ncol(theta.chain )
	colnames(a) <- paste0("a[", 1:I , "]")
	colnames(b) <- paste0("b[", 1:I , "]")
	mcmcobj <- cbind( deviance.chain , b )
	colnames(mcmcobj)[1] <- "deviance"
	if ( est.slope ){
		mcmcobj <- cbind( mcmcobj , a )
					}
	if ( est.guess ){
		colnames(c.chain) <- paste0("c[", 1:I , "]")		
		mcmcobj <- cbind( mcmcobj , c.chain )
					}
	if ( ! est.slope ){
		mcmcobj <- cbind( mcmcobj , sigma.chain )
		colnames(mcmcobj)[ncol(mcmcobj)] <- "sigma"
					  }
	if ( param==3 ){
		colnames(a.testlet.chain) <- paste0("a.testlet[", 1:I , "]")		
		ind <- which( testletgroups <= TT )
		mcmcobj <- cbind( mcmcobj , a.testlet.chain[,ind]  )
					}

					  
	colnames(theta) <- paste0("theta[", 1:N , "]")
	if ( ( TT>0 ) & (param!=3 ) ){	# save only testlet effects if they are present
		colnames(sigma.testlet.chain) <- paste0("sigma.testlet[", 1:TT , "]")	
		mcmcobj <- cbind( mcmcobj , sigma.testlet.chain )
			}
	#****
	# compute marginal effects
	if (TT>0){
		sigma.testlet.chain <- cbind( sigma.testlet.chain , 0 )
		ind <- which( testletgroups <= TT )
		# intercept parameters	
		st2 <- sqrt( 1 + sigma.testlet.chain[  , testletgroups[ind] ]^2 )
#		b_marg <- b.chain[,ind] * sigma.chain / st2 	
		b_marg <- b.chain[,ind] / st2 	
		colnames(b_marg) <- paste0("b_marg[", (1:I)[ind] , "]")
		mcmcobj <- cbind( mcmcobj , b_marg  )		
		# slope parameters
#		a_marg <- a.chain[,ind] * sigma.chain / st2 	
		a_marg <- a.chain[,ind] / st2 	
		colnames(a_marg) <- paste0("a_marg[", (1:I)[ind] , "]")
		mcmcobj <- cbind( mcmcobj , a_marg  )
		}									
	if (save.theta){ mcmcobj <- cbind( mcmcobj , theta ) }
	class(mcmcobj) <- "mcmc"
	attr(mcmcobj, "mcpar") <- c( burnin+1 , burnin+SV , 1 )
	mcmcobj <- coda::as.mcmc.list( mcmcobj )
	res <- list( "mcmcobj"=mcmcobj , "theta" = theta )
	return(res)
	}
#################################################################
# description
.mcmc.description.3pno.testlet <- function( TT , param , est.guess , est.slope){
		if ( est.guess ){ mod <- "3PNO" } else { mod <- "2PNO" }
		if ( (!est.guess)& (!est.slope) ){ mod <- "1PNO" } 
		if (TT>0){  
			description <- paste0( mod , " Testlet Model (param=" , param , ")" ) 
					}
		if (TT==0){ 
			description <- paste0( mod ," Model"  )
					}	
		return(description)
		}
##################################################################

####################################################################
# calculate EAP reliability, person parameter EAPs
# and corresponding posterior SDs
.mcmc.person.3pno.testlet <- function( theta.chain , weights , 
			gamma.testlet.chain){	
	###################
	# EAP reliability
	    N <- ncol(theta.chain )
		v1 <- stats::var( colMeans( theta.chain ) )
		SV <- nrow(theta.chain)
		if ( is.null(weights) ){ 				
				h1 <- ( rowSums( theta.chain^2 ) - ( N * rowMeans( theta.chain ) )^2 ) / N 
				h1 <- mean(h1)
				EAP.rel <- v1 / h1 
								}
		if ( ! is.null(weights) ){ 
				w1 <- weights / sum(weights )
				m1 <- colMeans( theta.chain )
				v1 <- sum( m1^2 * w1 ) - ( sum( m1*w1 ) )^2
				wM <- matrix( w1 , nrow=nrow(theta.chain) , ncol=ncol(theta.chain) , byrow=TRUE )			
				h1 <- rowSums( wM * theta.chain^2 ) - ( rowSums( wM * theta.chain ) )^2
				h1 <- mean(h1)
				EAP.rel <- v1 / h1			
						}	
	###################
	# person parameter estimates
	person <- data.frame( "EAP" = colMeans( theta.chain ) ,
						  "SD" = colSds(theta.chain) )
	# compute testlet effects
    TT <- ncol(gamma.testlet.chain) / N - 1
    if (TT>0){	
	for (tt in 1:TT){
		# tt <- 1
		gtc <- gamma.testlet.chain[ , 1:N + (tt-1)*N ]
		person[ , paste0("EAP.Testlet",tt)] <- colMeans(gtc)
		person[ , paste0("SD.Testlet",tt)] <- colSds( gtc )
					}						  
			}
	# output
	res <- list( "EAP.rel" = EAP.rel , "person" = person )
	return(res)
			}
###############################################################