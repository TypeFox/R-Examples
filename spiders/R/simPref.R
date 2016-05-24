##' @title simulate data
##'
##' @description  simulate data for predator preferences model
##'
##' @details Both lambda and gamma must be specified as a matrix with rows indexing
##' time and columns indexing the number of species.
##'
##' @return A list consisting of two dataframes, eaten and caught, made
##' specifically for the function \code{predPref}.
##'
##' @param S number of prey species
##' @param T number of time periods
##' @param J scalar or vector (of length T) number of predators caught at each time
##' @param I scalar or vector (of length T) effective number of traps at each time
##' @param lambda matrix of rates at which predator eats prey species; TxS
##' @param gamma matrix of rates at which prey species is seen in habitat; TxS
##' @param EM boolean specifying test of EM algorithm
##' 
##' @seealso \code{\link{predPref}}
##' @export
simPref <- function(S, T, J, I, lambda, gamma, EM=F) {

    ## checks
    if ( !isTRUE(all.equal(J, round(J))) ) {
        J <- round(J)
    }
    if ( !isTRUE(all.equal(I, round(I))) ) {
        I <- round(I)
    }
    
    ## some numbers
    ns <- seq_len(S)                    # index prey species
    nt <- seq_len(T)                    # index times
    lJ <- length(J)
    lI <- length(I)
    if ( lJ != 1 && lJ != T ) {
        stop("J not a scalar or vector of appropriate size.")
    }
    if ( lI != 1 && lI != T ) {
        stop("J not a scalar or vector of appropriate size.")
    }
    if (lJ == 1) {
        J <- rep(J[1], T)
    }
    if (lI == 1) {
        I <- rep(I[1], T)
    }

    ## initialize data frame
    eaten <- as.data.frame(matrix(NA, nrow=sum(J), ncol=S+1))
    caught <- as.data.frame(matrix(NA, nrow=sum(I), ncol=S+1))

    ## fill times
    colnames(eaten) <- colnames(caught) <- c("time", paste("preySpecies", ns, sep=""))
    eaten[,1] <- unlist(lapply(nt, function(x) rep(paste("time", x, sep=""), J[x])))
    caught[,1] <- unlist(lapply(nt, function(x) rep(paste("time", x, sep=""), I[x])))

    ## fill data frame
    times <- unique(eaten$time)
    for ( i in ns ) {
        for ( j in nt ) {
            jdx <- which(eaten$time == times[j])
            idx <- which(caught$time == times[j])
            eaten[jdx,i+1] <- rpois(J[j], lambda[j,i])
            caught[idx,i+1] <- rpois(I[j], gamma[j,i])
        }
    }

    ## if EM set, produce only binary observations
    if (EM) {
        jdx <- 2:(S+1)
        eaten[,jdx][which(eaten[,jdx]>0, arr.ind=T)] <- 1
    }
    
    out <- list(eaten=eaten,
                caught=caught)
    out
}
