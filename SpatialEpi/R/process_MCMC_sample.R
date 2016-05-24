process_MCMC_sample <-
function(sample, param, RR.area, cluster.list, cutoffs){

n <- length(RR.area)
n.sim <- length(sample)
  
# Store outputs here. These are averages weighted by the probability of 
# each configuration
high.area <- low.area <- RR.est.area <- rep(0, n)

# Zone based measures
high.z <- pgamma(cutoffs$high, param$shape, param$rate, lower.tail=FALSE)
low.z <- pgamma(cutoffs$low, param$shape, param$rate, lower.tail=TRUE)
RR.z <- param$shape / param$rate


#-------------------------------------------------------------------------------
# Find probabilities of cluster and anti-cluster membership & estimated 
# relative risk
#-------------------------------------------------------------------------------
# Get probabilities of each configuration
unique <- unique(sample)
indices <- match(sample, unique)
prob <- table(indices)/n.sim
prob.names <- as.numeric(names(prob))

# Over all sampled configurations
for(i in 1:length(prob.names)){
  config <- unique[[prob.names[i]]]
  k <- length(config)
  
  if(k==0){
    RR.est.area <- RR.est.area + RR.area * prob[i]
  }else{
    # Used to identify all areas included in configuration
    config.areas <- NULL    
    for(j in 1:k){        
      # All areas in this single zone in this configuration
      zone <- config[j]
      areas <- cluster.list[[zone]]
      config.areas <- c(config.areas, areas)
      
      # Update all areas inside configuration
      high.area[areas] <- high.area[areas] + high.z[zone] * prob[i]
      low.area[areas] <- low.area[areas] + low.z[zone] * prob[i]
      RR.est.area[areas] <- RR.est.area[areas] + RR.z[zone] * prob[i] 
    }
    
    # Update all areas outside configuration
    outside.areas <- setdiff(1:n, config.areas)
    RR.est.area[outside.areas] <- 
      RR.est.area[outside.areas] + RR.area[outside.areas]*prob[i]
  }
}

return(list(
  high.area=high.area,
  low.area=low.area,
  RR.est.area=RR.est.area)
  )
}
