multivspost2blog <-
function(smat, alpha, beta, D, y, By, k )
{
# computes posterior marginal of S in He/Hodges model
# called by optimization function 
# alpha = c(alphae, alphaz_1 to alphaz_{F-1})
# beta = c(betae, betaz_1 to betaz_{F-1})

# D:  matrix of diagonals of the diagonal matrices from diag of Q matrices
# k is rank deficiency of sum(Q's); i.e., rank of Q is n-k
# smat is matrix with F-1 cols; each row is s_1 to s_{F-1}

smat <- as.matrix(smat)
n <- length(y)

F <- ncol(D) + 1

sums <- apply(D,1,sum)
#k <- length( sums[sums==0] )

logpostvect <- numeric()
tausqtot <- numeric()

sumalpha <- sum(alpha)

        s <- smat
        s0 <- 1-sum(s)

        neweigennumer <- s[1] * D[,1]
        if (F > 2)
        for(j in 2:(F-1))
                neweigennumer <- neweigennumer + s[j] * D[,j]
        neweigendenom <-  neweigennumer + s0

        neweigen <- s0 * neweigennumer / neweigendenom
# corrected to (alpha-1) 09/18/09
        #logpostdensnumer <- sum( log(c(s0,s)) * alpha ) + 
        logpostdensnumer <- sum( log(c(s0,s)) * (alpha-1) ) + 
                        sum( log(neweigen[ neweigen > 0 ]) ) / 2
#       whole <- sum( neweigen * Bysq )
       whole <- sum( neweigen * By^2 )

        newbeta <- whole/2 +  sum( c(s0,s) * beta )
        newalpha <- (sumalpha + (n-k)/2 )
        logpostdensdenom <-  log( newbeta)* newalpha
        logpostvect <- logpostdensnumer - logpostdensdenom

logpostvect
 
}

