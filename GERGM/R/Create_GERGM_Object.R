# gergm object constructor (formerly gergm.object)
Create_GERGM_Object <- function(network,
                                bounded.network = network,
                                formula,
                                thetas,
                                lambda,
                                alpha,
                                together = together,
                                possible.stats,
                                thresholds = NULL) {
  if(is.null(thresholds)){
    thresholds <- rep(0, length(possible.stats))
  }
  num_stats <- length(possible.stats)
  num.nodes <- nrow(network)
  triples <- t(combn(1:num.nodes, 3))
  h.statistics1 <- h2(network,
                      triples = triples,
                      statistics = rep(1, num_stats),
                      together = together,
                      threshold = thresholds)
  h.statistics2 <- h2(bounded.network,
                      triples = triples,
                      statistics = rep(1, num_stats),
                      together = together,
                      threshold = thresholds)
  statistics <- rbind(h.statistics1, h.statistics2)
  colnames(statistics) <- possible.stats
  rownames(statistics) <- c("network", "bounded.network")
  # Check whether or not the weights are NULL
  if(is.null(alpha) == TRUE){
    alpha = 0
  }
  new("gergm", network = network,
      bounded.network = bounded.network,
      formula = formula,
      stats = statistics,
      theta.coef = thetas,
      lambda.coef = lambda,
      weights = alpha,
      num_nodes = num.nodes,
      thresholds = thresholds)
}
