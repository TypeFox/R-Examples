multivspost2blogandtausqtot <-
function(smat, alpha, beta, D, y, By, k )

{
# function to evaluate log marginal posterior of s and generate tausqtot

smat <- as.matrix(smat)

n <- length(y)

F <- ncol(D) + 1

logpostvect <- numeric()
tausqtot <- numeric()

sumalpha <- sum(alpha)


# change to use apply 09/29/09 MKC

genst <- function( s )
{
        s0 <- 1-sum(s)

        neweigennumer <- s[1] * D[,1]
        if( F > 2)
        for(j in 2:(F-1))
                neweigennumer <- neweigennumer + s[j] * D[,j]

        neweigendenom <-  neweigennumer + s0

        neweigen <- s0 * neweigennumer / neweigendenom
# corrected to (alpha - 1)  09/18/09
        #logpostdensnumer <- sum( log(c(s0,s)) * alpha  ) + 
        logpostdensnumer <- sum( log(c(s0,s)) * (alpha-1)  ) + 
                        sum( log(neweigen[ neweigen > 0 ]) ) / 2
        whole <- sum( neweigen * By^2)

        newbeta <- whole/2 +  sum( c(s0,s) * beta )
        newalpha <- (sumalpha + (n-k)/2 )
        logpostdensdenom <-  log( newbeta)* newalpha
        #logpostvect[i] <- logpostdensnumer - logpostdensdenom
        #tausqtot[i] <- rgamma(1, newalpha, newbeta)
        c( logpostdensnumer - logpostdensdenom, 
        rgamma(1, newalpha, newbeta))
}


#if( multicoreflag )
#    output <- t( parApply(cl, smat, 1, genst ) )
#else
    output <- t( apply( smat, 1, genst ) )
output
 
}

