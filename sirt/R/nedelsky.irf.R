


#####################################################################
# Item response function for the Nedelsky model
nedelsky.irf <- function( Theta , K , b , a , tau , combis , 
        thdim=1  ){
    # probabilities for category response function
    C1 <- nrow(combis)
    TP <- nrow(Theta)
    prob.cats <- matrix(NA,nrow=TP , ncol=K)
    for (cc in 1:K){
        # cc <- 1
        prob.cats[,cc] <- stats::plogis( a*Theta[,thdim ] - b[cc] )
                }
    prob.latclasses <- matrix(1 , nrow=TP , ncol=C1)
    # probabilities of latent classes
    for (kk in 1:K){
        # kk <- 1
        p1 <- prob.cats[,kk]
        p1 <- cbind( 1 - p1 , p1 )
        p1M <- p1[ , as.vector(combis[,kk])+1 ]
        prob.latclasses <- prob.latclasses * p1M
                    }
    combs <- 1 - cbind( 0 , combis )
    combs <- combs*matrix( tau , nrow=C1 , ncol=K+1 , byrow=TRUE )
    combs <- combs / rowSums(combs)
    # calculate category probabilities
    probs <- prob.latclasses %*% combs
	eps <- 1E-4
	probs[ probs < 0 ] <- eps
	probs[ probs > 1 ] <- 1 - eps 
	probs <- probs + eps
	probs <- probs / rowSums(probs )
    res <- list("probs"=probs , "prob.latclasses"=prob.latclasses )
    return(res)
        }
################################################################
