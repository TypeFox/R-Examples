# simulateMSAR.R
#
# simulate.msar() returns list of simulated data from an
# MS(h)AR(p) d.g.p.
# first element in list is simulated data
# second element is simulated regimes/states
#
# 20120113 : Initial version by Ryan Davis
# 20120120 : Renamed function so that is does not conflict with S3
#            generic simulate.*

simulateMSAR <- function(bigt, Q, theta, st1, y1)
{

# Q = h x h transition matrix

    h <- ncol(Q)
    p <- ncol(theta) - 2 # subtract off intercept and variance

# sanity checks on transitions matrix
    if ((nrow(Q) != h) || (ncol(Q) != h)) stop('Number of rows or columns does not equal h')
    if (sum(((Q >= 0) & (Q <= 1))) != h*h) stop('Probabilities are not between zero and one')
    if (sum(rowSums(Q)) != h) stop('Rows of transition matrix Q do not sum to unity).')
# sanity checks on theta
    if (nrow(theta) != h) stop('Number of parameters in theta does not match number of states')
    if (ncol(theta) != (p+2)) stop('Number of parameters in theta does not match AR(p)')


# setup vector to store states/regimes
    st.sim <- rep(NA, bigt)
    st.sim[1] <- st1  # starting regime

    for (i in 2:bigt) {

  # given the previous state was j, the only relevant
  # probabilities are those for the j'th row of transition matrix Q
  # use cumulative probs so very easy to use one random number
  # to determine state

        p.cumsum <- cumsum(Q[st.sim[i-1],])  # get cumulative probs
        u <- runif(1)                        # draw random number

  #sum(.) sums up an h x 1 vector of T/F based on cum probs
        st.sim[i] <- h - sum(u <= p.cumsum) + 1
}

# draw disturbances all at oncee
# bigt x h matrix
    emat <- NULL
    for (i in 1:h) {
  # rnorm takes st.dev., so sqrt(.)
        emat <- cbind(emat, rnorm(bigt, 0, sqrt(theta[i,p+2])))
    }

# setup vector for data
    y.sim <- rep(NA, bigt)
    y.sim[1:p] <- y1         # starting value

#
    for (i in (p+1):bigt) {
  # get current state, and then choose parameter
  # values based on regime
        s <- st.sim[i]
        ylag <- matrix(c(1, y.sim[(i-1):(i-p)]))  # bind intercept
        y.sim[i] <- crossprod(theta[s,1:(p+1)], ylag) + emat[i,s]
    }

    return(list(Y=y.sim, st=st.sim))

}
