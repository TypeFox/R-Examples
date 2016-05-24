# simulateMSVAR.R
#
# simulate.msvar() returns list of simulated data from an
# MS(h)VAR(p) d.g.p.
# first element in list is simulated data
# second element is simulated regimes/states

# 20120113 : Initial version in package, by Ryan Davis.
# 20120120 : Renamed function so that is does not conflict with S3
#            generic simulate.*

# e.vcv: m x m x h matrix containing variance-covariance
#        matrices for each regime

simulateMSVAR <- function(bigt, m, p, var.beta0, var.betas, e.vcv, Q,
                           seed=214)
{

# Q = h x h transition matrix

    h <- ncol(Q)
    st1 <- 1  # start in this state

# sanity checks on transitions matrix
    if ((nrow(Q) != h) || (ncol(Q) != h)) stop('Number of rows or columns does not equal h')
    if (sum(((Q >= 0) & (Q <= 1))) != h*h) stop('Probabilities are not between zero and one')
    if (sum(rowSums(Q)) != h) stop('Rows of transition matrix Q do not sum to unity).')

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
# bigt x m x h array
    e.rmvn <- array(NA, c(bigt, m, h))
    for (i in 1:h) {
  # rmvnorm takes var-covar matrix (not st.dev.)
        e.rmvn[,,i] <- rmvnorm(bigt, rep(0,m), e.vcv[,,i])
    }

# setup vector for data
    Y.sim <- matrix(NA, bigt, m)

# set first p observations to intercept
    Y.sim[1:p,] <- matrix(rep(var.beta0[,,st1], p), p, byrow=TRUE)

    for (i in p:(bigt-1)) {
  # get current state, and then choose parameter
  # values based on regime
        s <- st.sim[i]
        Ylag <- matrix(c(t(Y.sim[i:(i-p+1),])), 1)
        Y.sim[i+1,] <- var.beta0[,,s] + tcrossprod(var.betas[,,s], Ylag) + c(e.rmvn[i,,s])
}

    return(list(Y=Y.sim, st=st.sim))

}
