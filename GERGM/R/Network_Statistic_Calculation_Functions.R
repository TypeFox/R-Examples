# Network Statistics in R

# Statistics out2stars
out2star <- function(net, triples, alpha = 1, together = together) {
  if(together == 0){
    st1 <- sum(net[triples[, c(1, 2)]]^(alpha) * net[triples[, c(1, 3)]]^(alpha))
    st2 <- sum(net[triples[, c(2, 1)]]^(alpha) * net[triples[, c(2, 3)]]^(alpha))
    st3 <- sum(net[triples[, c(3, 1)]]^(alpha) * net[triples[, c(3, 2)]]^(alpha))
    return(st1 + st2 + st3)
  }
  if(together == 1){
    st1 <- sum(net[triples[, c(1, 2)]] * net[triples[, c(1, 3)]])
    st2 <- sum(net[triples[, c(2, 1)]] * net[triples[, c(2, 3)]])
    st3 <- sum(net[triples[, c(3, 1)]] * net[triples[, c(3, 2)]])
    return((st1 + st2 + st3)^alpha)
  }
}

#-------------------------------------------------------
# in2stars
in2star <- function(net, triples, alpha = 1, together = together) {
  if(together == 0){
    st1 <- sum(net[triples[, c(3, 1)]]^(alpha) * net[triples[, c(2, 1)]]^(alpha))
    st2 <- sum(net[triples[, c(3, 2)]]^(alpha) * net[triples[, c(1, 2)]]^(alpha))
    st3 <- sum(net[triples[, c(1, 3)]]^(alpha) * net[triples[, c(2, 3)]]^(alpha))
    return(st1 + st2 + st3)
  }
  if(together == 1){
    st1 <- sum(net[triples[, c(3, 1)]] * net[triples[, c(2, 1)]])
    st2 <- sum(net[triples[, c(3, 2)]] * net[triples[, c(1, 2)]])
    st3 <- sum(net[triples[, c(1, 3)]] * net[triples[, c(2, 3)]])
    return((st1 + st2 + st3)^alpha)
  }
}

#-------------------------------------------------------
# transitive triads
ttriads <- function(net, triples, alpha = 1, together) {
  if(together == 0){
    t2 <- sum(net[triples[, c(1, 2)]]^(alpha) * net[triples[, c(2, 3)]]^(alpha) *
                net[triples[, c(1, 3)]]^(alpha))
    t3 <- sum(net[triples[, c(1, 2)]]^(alpha) * net[triples[, c(3, 2)]]^(alpha) *
                net[triples[, c(3, 1)]]^(alpha))
    t4 <- sum(net[triples[, c(1, 2)]]^(alpha) * net[triples[, c(3, 2)]]^(alpha) *
                net[triples[, c(1, 3)]]^(alpha))
    t5 <- sum(net[triples[, c(2, 1)]]^(alpha) * net[triples[, c(2, 3)]]^(alpha) *
                net[triples[, c(3, 1)]]^(alpha))
    t6 <- sum(net[triples[, c(2, 1)]]^(alpha) * net[triples[, c(2, 3)]]^(alpha) *
                net[triples[, c(1, 3)]]^(alpha))
    t7 <- sum(net[triples[, c(2, 1)]]^(alpha) * net[triples[, c(3, 2)]]^(alpha) *
                net[triples[, c(3, 1)]]^(alpha))
    return(t2 + t3 + t4 + t5 + t6 + t7)
  }
  if(together == 1){
    t2 <- sum(net[triples[, c(1, 2)]] * net[triples[, c(2, 3)]] *
                net[triples[, c(1, 3)]])
    t3 <- sum(net[triples[, c(1, 2)]] * net[triples[, c(3, 2)]] *
                net[triples[, c(3, 1)]])
    t4 <- sum(net[triples[, c(1, 2)]] * net[triples[, c(3, 2)]] *
                net[triples[, c(1, 3)]])
    t5 <- sum(net[triples[, c(2, 1)]] * net[triples[, c(2, 3)]] *
                net[triples[, c(3, 1)]])
    t6 <- sum(net[triples[, c(2, 1)]] * net[triples[, c(2, 3)]] *
                net[triples[, c(1, 3)]])
    t7 <- sum(net[triples[, c(2, 1)]] * net[triples[, c(3, 2)]] *
                net[triples[, c(3, 1)]])
    return((t2 + t3 + t4 + t5 + t6 + t7)^alpha)
  }
}

#-------------------------------------------------------
# cyclic triads
ctriads <- function(net, triples, alpha = 1, together) {
  if(together == 0){
    t1 <- sum(net[triples[, c(1, 2)]]^(alpha) * net[triples[, c(2, 3)]]^(alpha) *
                net[triples[, c(3, 1)]]^(alpha))
    t8 <- sum(net[triples[, c(2, 1)]] * net[triples[, c(3, 2)]] *
                net[triples[, c(1, 3)]]^(alpha))
    return(t1 + t8)
  }
  if(together == 1){
    t1 <- sum(net[triples[, c(1, 2)]] * net[triples[, c(2, 3)]] *
                net[triples[, c(3, 1)]])
    t8 <- sum(net[triples[, c(2, 1)]] * net[triples[, c(3, 2)]] *
                net[triples[, c(1, 3)]])
    return((t1 + t8)^alpha)
  }
}

#-------------------------------------------------------
# reciprocity
recip <- function(net, alpha = 1, together) {
  pairs <- t(combn(1:nrow(net), 2))
  if(together == 0){
    return(sum(net[pairs]^(alpha) * net[pairs[, c(2, 1)]]^(alpha)))
  }
  if(together == 1){
    return(sum(net[pairs] * net[pairs[, c(2, 1)]])^(alpha))
  }
}

#-------------------------------------------------------
# edgeweight
edgeweight <- function(net, alpha = 1, together) {
  pairs <- t(combn(1:nrow(net), 2))
  if(together == 0){
    return(sum((net[pairs]^(alpha) + net[pairs[, c(2, 1)]])^(alpha)))
  }
  if(together == 1){
    return(sum((net[pairs] + net[pairs[, c(2, 1)]]))^(alpha))
  }
}

#-------------------------------------------------------

absdiff <- function(attrname, pow, directed){
  # sum of the absolute difference of the nodal attribute "attrname" between
  # every pair of nodes
  diff <- numeric()
  for(i in 1:length(attrname)){
    for(j in 1:length(attrname)){
      diff <- diff + abs(attrname[i] - attrname[j])^pow
    }
  }
  if(directed == FALSE){
    diff <- diff/2
  }
  return(diff)
}
#-------------------------------------------------------

atleast <- function(network, threshold, directed){
  # number of edges that have weight that exceed specified threshold if
  # threshold = 0, this is just the number of edges
  stat <- length(which(network > threshold | network == threshold))
  if(directed == FALSE){
    stat <- stat/2
  }
  return(stat)
}

#-------------------------------------------------------

degree_atleast = function(network, threshold, directed){
  # the sum of edgeweights that exceed specified threshold if threshold = 0,
  # this is just the total edgeweight of the network
  directed <- directed[1]
  stat <- sum(network[which(network > threshold | network == threshold)])
  if(directed == FALSE){
    stat <- stat/2
  }
  return(stat)
}

#-------------------------------------------------------

density_atleast = function(network, threshold = 0, directed = c(TRUE, FALSE)){
  # the sum of edgeweights that exceed specified threshold if threshold = 0,
  # this is just the total edgeweight of the network
  directed <- directed[1]
  stat <- sum(network[which(network > threshold | network == threshold)])/length(which(network > threshold | network == threshold))
  return(stat)
}

#-------------------------------------------------------

mutual = function(network, threshold, form){
  # mutuality statistic calculated for each pair of mutual edges that have
  # edgeweights that exceed threshold. The form specifies which statistic is
  # calculated and can take the following values: form = c("min", "product",
  # "geometric") DEFAULT is "min". This can only be used for a directed
  # network.
  form <- form[1]
  indx <- which(network < threshold)
  network[indx] <- 0
  values <- c(0,0)
  for(i in 1:dim(network)[1]){
    for(j in 1:dim(network)[1]){
      if(i != j){
        values <- rbind(values, c(network[i,j], network[j,i]))
      }
    }
  }
  if(form == "min"){
    stat <- sum(apply(values, 1, min))
  }
  if(form == "product"){
    stat <- sum(values[, 1]*values[, 2])
  }
  if(form == "geometric"){
    stat <- sum(sqrt(values[, 1])*sqrt(values[, 2]))
  }
  return(stat)
}

#-------------------------------------------------------

nodecov = function(network, attrname, threshold){
  # statistic that adds the sum of attribute[i] and attribute[j] for all edges
  # (i,j) such that its corresponding edgeweight exceeds the specified
  # threshold NOTE: attrname must be numeric (not categorical)
  indx <- which(network > threshold | network == threshold, arr.ind = TRUE)
  stat <- 0
  for(i in 1:dim(indx)[1]){
    stat <- stat + attrname[indx[i,1]] + attrname[indx[i,2]]
  }
  return(stat)
}

#-------------------------------------------------------

nodefactor = function(attrname, base = 1){
  # adds several statistics where each statistic counts the number of instances
  # attrname takes a discrete (categorical) value. The base specifies which
  # categorical value should not be counted. For instance, base = 1, there will
  # be no statistic for the number of times 1 occurs. It is recommended to
  # avoid counting all values NOTE: attrname must be categorical (not numeric)
  attrname.factor <- as.factor(attrname)
  names <- unique(attrname.factor)
  num.vars <- length(unique(attrname.factor))
  stat <- numeric()
  for(i in 1:num.vars){
    stat[i] <- length(which(attrname.factor == names[i]))
  }
  indx <- which(1:num.vars != base)
  result <- data.frame(stat[indx])
  rownames(result) <- unique(attrname)[indx]
  colnames(result) <- "Count"
  result <- t(result)
  return(result)
}

#-------------------------------------------------------

#undirected statistic for triads
undirected_triads = function(net, triples, alpha = 1, together) {
    if(together == 0){
      stat <- sum(net[triples[, c(1, 2)]]^(alpha) * net[triples[, c(2, 3)]]^(alpha) *
                  net[triples[, c(1, 3)]]^(alpha))
      return(stat)
    }
    if(together == 1){
      stat <- sum(net[triples[, c(1, 2)]] * net[triples[, c(2, 3)]] *
                  net[triples[, c(1, 3)]])
      return(stat^alpha)
    }
  }

#-------------------------------------------------------

# Calculate the statistics of a formula object
h <- function(possible.stats,
              alpha = NULL,
              theta = NULL,
              together = together,
              GERGM_Object = GERGM_Object) {
  net <- GERGM_Object@observed_bounded_network
  alphas <- GERGM_Object@weights
  thresholds <- GERGM_Object@thresholds
  num.nodes <- GERGM_Object@num_nodes
  statistics <- GERGM_Object@stats_to_use
  triples = t(combn(1:num.nodes, 3))
  temp <- c(out2star(net, triples, alphas[1], together),
            in2star(net, triples, alphas[2], together),
            ctriads(net, triples, alphas[3], together),
            recip(net, alphas[4], together),
            ttriads(net, triples, alphas[5], together),
            edgeweight(net, alphas[6], together))
  # check to make sure we did not mess things up
  if(length(temp) != length(statistics)){
    stop("Development ERROR! Please email mdenny@psu.edu! The h() internal function in Network_Statistic_Calculation_Functions.R has been supplied an incorrect number of statistics.")
  }
  value <- temp[statistics > 0]
  result <- rbind(round(value, 3), round(alphas[statistics > 0], 3))
  colnames(result) <- possible.stats[statistics > 0]
  rownames(result) <- c("value", "alpha")
  return(result)
}

#h.corr for correlation matrices, we want to use the unbounded networks
h.corr <- function(possible.stats,
              alpha = NULL,
              theta = NULL,
              together = together,
              GERGM_Object = GERGM_Object) {
  net <- GERGM_Object@observed_network
  alphas <- GERGM_Object@weights
  thresholds <- GERGM_Object@thresholds
  num.nodes <- GERGM_Object@num_nodes
  statistics <- GERGM_Object@stats_to_use
  triples = t(combn(1:num.nodes, 3))
  temp <- c(out2star(net, triples, alphas[1], together),
            in2star(net, triples, alphas[2], together),
            ctriads(net, triples, alphas[3], together),
            recip(net, alphas[4], together),
            ttriads(net, triples, alphas[5], together),
            edgeweight(net, alphas[6], together))
  # check to make sure we did not mess things up
  if(length(temp) != length(statistics)){
    stop("Development ERROR! Please email mdenny@psu.edu! The h.corr() internal function in Network_Statistic_Calculation_Functions.R has been supplied an incorrect number of statistics.")
  }
  value <- temp[statistics > 0]
  result <- rbind(round(value, 3), round(alphas[statistics > 0], 3))
  colnames(result) <- possible.stats[statistics > 0]
  rownames(result) <- c("value", "alpha")
  return(result)
}

# A second version of the h function used for calculation in estimation
# This function calculates the network statistics associated with net
h2 <- function(net,
               triples,
               statistics,
               alphas = NULL,
               together = 1,
               directed = TRUE,
               threshold = NULL) {
  if(is.null(threshold)){
    thresholds <- rep(0, length(statistics))
  }
  if(is.null(alphas)){
    alphas <- rep(1,length(statistics))
  }
  temp = c(out2star(net, triples, alphas[1], together),
           in2star(net, triples, alphas[2], together),
           ctriads(net, triples, alphas[3], together),
           recip(net, alphas[4], together),
           ttriads(net, triples, alphas[5], together),
           edgeweight(net, alphas[6], together))
  # check to make sure we did not mess things up
  if(length(temp) != length(statistics)){
    stop("Development ERROR! Please email mdenny@psu.edu! The h2() internal function in Network_Statistic_Calculation_Functions.R has been supplied an incorrect number of statistics.")
  }
  return(temp[which(statistics > 0)])
  #return(temp)
}
