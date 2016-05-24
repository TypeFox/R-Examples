# Simulation Functions
Gibbs_Sampler <- function(GERGM_Object,
                          theta,
                          MCMC.burnin,
                          num.draws,
                          thin = 1,
                          start = NULL,
                          num.nodes = NULL,
                          directed,
                          possible.stats) {
  # formula specifies which variables to use MCMC.burnin is the number of
  # discarded draws n is the total number of draws to be reported dh is a
  # function that takes as arguments, i, j, theta, net and returns the partial
  # of the hamiltonian theta is the vector-valued parameter thin reduces
  # autocorrelation in the simulations, every thinth iteration is returned start
  # is the initial network, if not supplied, a random uniform nodesXnodes
  # network is used num.nodes is the number of nodes in the network dir is a
  # logical indicator of whether the network is directed

  net <- GERGM_Object@network
  if (is.null(num.nodes) == TRUE) {
    num.nodes <- GERGM_Object@num_nodes
  }

  statistics <- GERGM_Object@stats_to_use

  sims <- num.draws * thin
  netarray <- array(NA, dim = c(num.nodes, num.nodes, num.draws + 1))
  if (is.null(start))
    start <- matrix(rdisp(num.nodes * num.nodes), num.nodes, num.nodes)
  net <- start
  diag(net) = 0
  netarray[, , 1] <- net
  for (t in 1:(num.draws + MCMC.burnin)) {
    if (directed == TRUE) {
      for (i in 1:num.nodes) {
        for (j in (1:num.nodes)[-i]) {
          net[i, j] <- rtexp(1, t(theta) %*% dh(net, statistics, i, j))
        }
      }
    }
    if (directed == FALSE) {
      for (i in 2:num.nodes) {
        for (j in (1:(i - 1))) {
          net[i, j] <- rtexp(1, t(theta) %*% dh(net, statistics, i, j))
        }
      }
    }
    if (t > MCMC.burnin) {
      netarray[, , t - MCMC.burnin] <- net
    }
  }
  return(netarray[, , round(seq(1, sims, length = num.draws))])
}
