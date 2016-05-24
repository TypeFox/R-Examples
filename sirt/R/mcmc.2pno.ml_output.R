
######################################################				
# subfunction for calculating the DIC
.mcmc.ic.2pno.ml <- function( aM.chainsum , bM.chainsum , theta.chain , N , I , 
		dat , dat.resp ,eps , deviance.chain , groupsize ,
		sigma.res.chain , link ){
	#************
	weights <- NULL
#	aM <- matrix( colMeans( a.chain ) , nrow=N , ncol=I , byrow=TRUE )
#	bM <- matrix( colMeans( b.chain ) , nrow=N , ncol=I , byrow=TRUE )
	SV <- nrow(theta.chain )
	aM <- aM.chainsum / SV
	bM <- bM.chainsum / SV
	theta <- colMeans( theta.chain )
	if (link=="logit"){
		Dhat <- .mcmc.deviance.2pl( aM , bM , theta , dat , dat.resp , 
					weights , eps )
						}
	if (link=="normal"){
		sigma.res <- as.vector(colMeans( sigma.res.chain ) )
		Dhat <- .mcmc.deviance.normallink.2pno.ml( aM , bM , 
						theta , N , I , dat , dat.resp , sigma.res )					
						}
	Dbar <- mean( deviance.chain )
	pD <- Dbar - Dhat
	ic <- list( "Dhat"=Dhat , "Dbar"=Dbar , "pD"=pD , "DIC" = Dhat + 2*pD )
	# group statistics
	ic$G <- length(groupsize)
	ic$M.n <- mean(groupsize)
	ic$SD.n <- stats::sd(groupsize)	
	return(ic)
		}
#######################################################		
# mcmc.list output object
.mcmc.mcmclist.2pno.ml <- function( a.chain , b.chain , deviance.chain ,
	sigma1.chain , sigma2.chain , I , SV , burnin , iter , est.b.M ,
	mu.b.chain , omega.b.chain , sigma.b.chain , est.b.Var ,
	est.a.M , est.a.Var , omega.a.chain , sigma.a.chain ,
	sigma.res.chain , link	){
	#****
	a <- a.chain 
	b <- b.chain
	colnames(a) <- paste0("a[", 1:I , "]")
	colnames(b) <- paste0("b[", 1:I , "]")
	mcmcobj <- cbind( deviance.chain , b )	
	# include a parameters
	if( max( apply( a , 2 , stats::sd) ) > 0 ){
		mcmcobj <- cbind( mcmcobj , a )
				}		
	colnames(mcmcobj)[1] <- "deviance"
	mcmcobj <- cbind(  mcmcobj , "sigma1"=sigma1.chain ,
			"sigma2"=sigma2.chain )
	ICC <- sigma2.chain^2 / ( sigma2.chain^2 + sigma1.chain^2 )
	mcmcobj <- cbind( mcmcobj , "ICC"=ICC )
	if ( est.b.M == "h"){ mcmcobj <- cbind( mcmcobj , 		
				"mu.b"=mu.b.chain , "omega.b"=omega.b.chain ) 
						}
	#*** est.b.Var
    if ( est.b.Var %in% c("i") ){
		colnames(sigma.b.chain) <- paste0("sigma.b[" , 1:I,"]")
		mcmcobj <- cbind( mcmcobj , sigma.b.chain )
							}
    if ( est.b.Var == "j" ){
		mcmcobj <- cbind( mcmcobj , "sigma.b"= sigma.b.chain[,1] )
							}							
	#**** est.a.Var	
	if ( est.a.M != "f"){	
		if ( est.a.Var == "i" ){
			colnames(sigma.a.chain) <- paste0("sigma.a[" , 1:I,"]")
			mcmcobj <- cbind( mcmcobj , sigma.a.chain )
								}
		if ( est.a.Var == "j" ){
			mcmcobj <- cbind( mcmcobj , "sigma.a"= sigma.a.chain[,1] )
								}
						}		
	if ( est.a.M == "h"){ mcmcobj <- cbind( mcmcobj , 		
				 "omega.a"=omega.a.chain ) 
						}
    if ( link == "normal" ){
		colnames(sigma.res.chain) <- paste0("sigma.res[" , 1:I,"]")
		mcmcobj <- cbind( mcmcobj , sigma.res.chain )
							}						
	class(mcmcobj) <- "mcmc"
	attr(mcmcobj, "mcpar") <- c( burnin+1 , burnin+SV , 1 )
	mcmcobj <- coda::as.mcmc.list( mcmcobj )
    return(mcmcobj)
		}
######################################################
# calculating the deviance under the normal link function
.mcmc.deviance.normallink.2pno.ml <- function( aM , bM , 
		theta , N , I , dat , dat.resp , sigma.res ){							
	mexp <- aM*theta - bM							
	mres <- dat - mexp
	llres <- - sum( dat.resp * 
		log( stats::dnorm( mres , sd=matrix( sigma.res , N, I , byrow=TRUE ) )	) )
	return(llres)
		}		
		
		
####################################################################
# calculate EAP reliability, person parameter EAPs
# and corresponding posterior SDs
.mcmc.person.2pno.ml <- function( theta.chain , weights ){	
	###################
	# EAP reliability
		v1 <- stats::var( colMeans( theta.chain ) )
		if ( is.null(weights) ){ 
			weights <- rep(1,ncol(theta.chain)) }
		w1 <- weights / sum(weights )
		m1 <- colMeans( theta.chain )
		v1 <- sum( m1^2 * w1 ) - ( sum( m1*w1 ) )^2
		wM <- matrix( w1 , nrow=nrow(theta.chain) , ncol=ncol(theta.chain) , byrow=TRUE )			
		h1 <- rowSums( wM * theta.chain^2 ) - ( rowSums( wM * theta.chain ) )^2
		h1 <- mean(h1)
		EAP.rel <- v1 / h1			
	###################
	# person parameter estimates
	person <- data.frame( "EAP" = colMeans( theta.chain ) ,
						"SD" = colSds(theta.chain) )
	# output
	res <- list( "EAP.rel" = EAP.rel , "person" = person )
	return(res)
			}
###############################################################