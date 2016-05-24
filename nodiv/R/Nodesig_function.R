
### INTERNAL FUNCTIONS

pval_sig <- function(pval) 
  2*(0.5-abs(pval-0.5))

logit <- function(p)
{
  if(!(sum(p >= 1, na.rm = T) + sum(p <= 0, na.rm = T) == 0))
    warning("logit is only defined between the interval 0 and 1")
  p[p == 0] <- NA
  p[p == 1] <- NA
  return(log(p/(1-p)))
}

inv_logit <- function(a)
{
  return(exp(a)/(1 + exp(a)))
}


quasiswap_nodesig <- function(simcom, Node_sp, repeats, show )
{
  ll <- 0
  if (show) pb <- txtProgressBar(min = 0, max = repeats-1, style = 3)
  replicate(repeats-1,
  {
    if (show) setTxtProgressBar(pb, ll <<- ll + 1)
    simcom <- commsimulator(simcom, method = "quasiswap")
    rowSums(simcom[, Node_sp])
    })
}

rdtable_nodesig <- function(simcom, Node_sp, repeats, show )
{
  tempcom <- cbind(rowSums(simcom[,Node_sp]), rowSums(simcom[,-Node_sp]))
  r <- rowSums(tempcom)
  c <- colSums(tempcom)
  ll <- 0
  if (show) pb <- txtProgressBar(min = 0, max = repeats-1, style = 3)
  
  # Use the rdtable algorithm to created random matrices, basing each new matrix on the previous
  nodereps <- replicate(repeats-1,
  {
    # perform the randomization
    simcom <- r2dtable(1, r, c)[[1]]
    # return the simulated species richness of sites
    if (show) setTxtProgressBar(pb, ll <<- ll + 1)
    simcom[, 1]
    
  })
  return(nodereps)
}

### EXPORTED FUNCTIONS ########

Nodesig <- function(nodiv_data, Node_sp = NULL, repeats = 100, method = c("rdtable", "quasiswap"), show = T)
{
  if(!inherits(nodiv_data, "distrib_data"))
    stop("nodiv_data must be an object of type nodiv_data or distrib_data")
  
  if(is.null(Node_sp)) 
    if(inherits(nodiv_data, "nodiv_data"))
      Node_sp <- Node_species(nodiv_data, Descendants(basal_node(nodiv_data), nodiv_data)[1], names = FALSE) else
      stop("Node_sp must be defined if nodiv_data has type distrib_data")
  # a boolean vector indicating which of species descending from the parent node that descend from the focal node
  if(is.character(Node_sp))
    Node_sp <- which(species(nodiv_data) %in% Node_sp)  # make a boolean vector
  
  if(length(Node_sp) == 1 | length(Node_sp) == Nspecies(nodiv_data)-1 | Nsites(nodiv_data) < 2) 
    return(list(SR = rep(NA,Nsites(nodiv_data)),  rval = rep(NA,Nsites(nodiv_data)), nodeemp = rep(NA,Nsites(nodiv_data)), nodemeans = rep(NA,Nsites(nodiv_data)), nodesds = rep(NA,Nsites(nodiv_data)))) #if one of the descendant clades is a single species
 
  method = match.arg(method)
  # A global variable to count the number of repeats

  simcom <- nodiv_data$comm
  nodereps <- switch(method,
         quasiswap = quasiswap_nodesig(simcom, Node_sp, repeats, show),
         rdtable = rdtable_nodesig(simcom, Node_sp, repeats, show)
  )
  nodeemp <- rowSums(nodiv_data$comm[, Node_sp])
  nodereps <- cbind(nodeemp, nodereps)
  ord <- apply(nodereps, 1,  rank)[1,]
  ord[ord == repeats] <- repeats - 1
  ord[rowSums(nodiv_data$comm) == 0] <- NA
  
  #calculates the effect size (take care of sd = 0)!
  nodemeans <- rowMeans(nodereps)
  nodesds <- apply(nodereps, 1, sd)
  
  # all sites occupied by the parent node get an SES value - the rest retain NA
  ses <- (nodeemp - nodemeans )/ nodesds
  rval <-  ord/repeats
  
  return(list(SR = ses,  rval = rval, nodeemp = nodeemp, nodemeans = nodemeans, nodesds = nodesds))
}
