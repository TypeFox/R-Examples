
###############################
# attach objects in the local environment
.sirt.attach.environment <- function( res , envir ){
#	e1 <- environment()
	CC <- length(res)
	for (cc in 1:CC){
		assign( names(res)[cc] , res[[cc]] , envir=envir )		
					}
			}
#####################################
# Data processing GOM JML
.gom.proc <- function( dat , seed=NULL,K ){
    dat.resp <- 1-is.na(dat)
    dat[ is.na(dat) ] <- 0
    N <- nrow(dat)
    I <- ncol(dat)
    # create starting values for person parameters
    score <- rowSums( dat ) / rowSums( dat.resp )
    theta0 <- seq( -1.5 ,1.5 , len=K )
    g <- exp( - ( outer( stats::qlogis( ( score + .1 ) / 1.2 ) , theta0 , "-" ) )^2 )
    g <- g / rowSums(g)
    # create starting values item parameters
    p.item <- colSums( dat ) / colSums( dat.resp )
    theta0 <- seq( -2 ,2 , len=K )
    lambda <- t( stats::plogis( outer( theta0 , - stats::qlogis( p.item ) , "-" ) ) )
	# random start
	if ( ! is.null(seed) ){
		g <- matrix( stats::runif( N*K ) , N , K )
		g <- g / rowSums( g)
		lambda <- matrix( stats::runif(I*K) , nrow=I , ncol=K)		
			}
	
    res <- list( "dat"=dat , "dat.resp"=dat.resp , "N"=N , "K"=K , 
        "lambda"=lambda , "g"=g , "I"=I , "score" = score , "p.item"=p.item )
    return(res)
        }
######################################
# GOM JML update g
.gom.jml.est.g <- function( lambda , g , N , I , K ,  dat , dat.resp , min.g){
    #***
    # update g
    # sum lambda
    lambda0.sum <- lambda1.sum <- 0
    for (kk in 1:K){
        # kk <- 1
        lambda1.sum <- lambda1.sum + matrix( lambda[,kk] , N , I , byrow=T ) * g[,kk ]
        lambda0.sum <- lambda0.sum + matrix( 1-lambda[,kk] , N , I , byrow=T ) * g[,kk ]    
                    }
    for (kk in 1:K){
        # kk <- 1
        lambda.kk <- matrix( lambda[,kk] , N , I , byrow=T ) 
        w1.kk <- lambda.kk * g[,kk] / lambda1.sum
        w0.kk <- ( 1 - lambda.kk ) * g[,kk]  / lambda0.sum
        g[,kk] <- rowSums( dat * dat.resp * w1.kk + ( 1- dat ) * dat.resp * w0.kk ) / rowSums( dat.resp )
                }						
	g[ g < min.g ] <- min.g
	g[ g > 1 - min.g ] <- 1 - min.g
	g <- g / rowSums(g)
    return( g )
        }
#################################################
# update lambda
.gom.jml.est.lambda <- function( lambda , g , N , I , K , dat , dat.resp , min.lambda ){
    lambda0.sum <- lambda1.sum <- 0
    for (kk in 1:K){
        # kk <- 1
        lambda1.sum <- lambda1.sum + matrix( lambda[,kk] , N , I , byrow=T ) * g[,kk ]
        lambda0.sum <- lambda0.sum + matrix( 1-lambda[,kk] , N , I , byrow=T ) * g[,kk ]    
                    }                                        
    for (kk in 1:K){
        # kk <- 1
        lambda.kk <- matrix( lambda[,kk] , N , I , byrow=T ) 
        w1.kk <- lambda.kk * g[,kk] / lambda1.sum
        w0.kk <- ( 1 - lambda.kk ) * g[,kk]  / lambda0.sum
        sum1 <- colSums( dat * dat.resp * w1.kk )
        sum0 <- colSums( (1-dat) * dat.resp * w0.kk )        
        lambda[,kk] <- sum1 / ( sum1 + sum0 )
                }				
    lambda[ lambda < min.lambda ] <- min.lambda			
	lambda[ lambda > 1 - min.lambda ] <- 1 - min.lambda	
    return(lambda)
        }
###############################################
# calculation of deviance
.gom.deviance <- function( lambda, g , K , dat , dat.resp ){
    prob0 <- prob1 <- 0
    for (kk in 1:K){
        # kk <- 1
        prob1 <- prob1 + outer( g[,kk] , lambda[,kk] )
        prob0 <- prob0 + outer( g[,kk] , 1-lambda[,kk] )
                }
    dev <- -2*sum( log( dat*dat.resp*prob1 + (1-dat)*dat.resp*prob0 ) )
    return(dev)
            }
##################################################
# Discrteizing Membership Function
.gom.discret.membership <- function(g,K,score){
	g1 <- g
	b1 <- seq(0,1,len=5)
	for (kk in 1:K){	g1[,kk] <- cut( g[,kk] , breaks=b1) }
	pattern <- apply( g1 , 1 , FUN = function(kk){ paste( kk , collapse="-") } )
	gcut.distr <- matrix(0,K , 2*(length(b1)-1) )
	for (kk in 1:K){
		for (ii in 1:(length(b1)-1) ){
			gcut.distr[kk,ii] <- sum( g1[,kk]== ii )
			gcut.distr[kk,ii+4] <- mean( score[ g1[,kk]== ii ] )
					}
			}
    gcut.distr[,1:4] <- 100*gcut.distr[,1:4] / rowSums( gcut.distr[,1:4] )
    gcut.distr <- data.frame(gcut.distr)	
	rownames(gcut.distr) <- paste0("Class",1:K)
	colnames(gcut.distr) <- c( paste0("PrGr" , seq(0,75,len=4) )	 ,
					paste0("MGr" , seq(0,75,len=4) )	 )
	gcut.distr$cor.score <- NA
	for (kk in 1:K){
		gcut.distr$cor.score[kk] <- cor( g[,kk] , score )
			}		
	g1 <- data.frame(g1)
	g1$pattern <- pattern
    g.mean <- colMeans(g)
	res <- list("g1"=g1 , "gcut.distr"=gcut.distr ,  "g.mean" = g.mean )
	return(res)
				}			

