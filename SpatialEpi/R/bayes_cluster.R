bayes_cluster <-
function(y, E, population, sp.obj, centroids, max.prop, shape, rate, J, pi0,
         n.sim.lambda, n.sim.prior, n.sim.post, burnin.prop = 0.1,
         theta.init = vector(mode="numeric", length=0)
){
print(paste("Algorithm started on:", date()))
  
# MCMC parameters
pattern <- c(0, 1)
p.moves <- normalize(c(1, 1, 1, 1, 1))
  
# Geographic info and create geographical objects to use
geo.objects <- create_geo_objects(max.prop, population, centroids, sp.obj)
overlap <- geo.objects$overlap
cluster.coords <- geo.objects$cluster.coords    
n <- length(E)
n.zones <- nrow(cluster.coords)

# Set cutoffs of the relative risk to determine low and high risk
RR.range <- seq(0.7, 1/0.7, length=100000)
pN <- dgamma(RR.range, shape[1], rate[1])
pW <- dgamma(RR.range, shape[2], rate[2])
cutoffs <- list(
  low=RR.range[min(which(pN/pW>1))],
  high=RR.range[max(which(pN/pW>1))]
)


#-------------------------------------------------------------------------------
# Obtain values of lambda using importance sampling
#-------------------------------------------------------------------------------
# Set uniform prior on single zones.
prior.z <- rep(1/n.zones, n.zones)
log_prior.z <- log(prior.z) - log(sum(prior.z))

# Estimate q and lambda
results <- estimate_lambda(n.sim.lambda, J, prior.z, overlap, pi0)
lambda <- results$lambda
prior.j <- results$prior.j
print(paste("Importance sampling of lambda complete on:", date()))


#-------------------------------------------------------------------------------
# Obtain prior map
#-------------------------------------------------------------------------------
# Generate MCMC samples from prior
prior.chain <- MCMC_simulation(n.sim.prior, pattern, theta.init, overlap, 
                               cluster.coords, p.moves, J, prior.z, lambda)

# Trim burn-in
burnin <- n.sim.prior * burnin.prop
prior.sample <- prior.chain$sample[-c(1:burnin)]

# Prior Probs of Cluster Membership for each Area
param.prior.zone <- list(shape=rep(shape[2],n.zones), rate=rep(rate[2],n.zones))
RR.prior.area <- rep(shape[1]/rate[1], n)
prior.map <- process_MCMC_sample(prior.sample, param.prior.zone, RR.prior.area, 
                                 overlap$cluster.list, cutoffs)
print(paste("Prior map MCMC complete on:", date()))


#-------------------------------------------------------------------------------
# Obtain posterior map
#-------------------------------------------------------------------------------
# Single Zone Values
yz <- sapply(overlap$cluster.list, function(x){sum(y[x])})
Ez <- sapply(overlap$cluster.list, function(x){sum(E[x])})
log_BF.z <- coeff(y, E, shape, rate, overlap$cluster.list)
BF.z <- exp(log_BF.z)
log_post.z <- log_prior.z + log_BF.z
# Do NOT normalize this quantity
post.z <- exp(log_post.z)

# Generate MCMC samples from prior
post.chain <- MCMC_simulation(n.sim.post, pattern, theta.init, overlap, 
                              cluster.coords, p.moves, J, post.z, lambda)

# Trim burn-in
burnin <- n.sim.post * burnin.prop
post.sample <- post.chain$sample[-c(1:burnin)]

# Prior Probs of Cluster Membership for each area
param.post.zone <- list(shape=c(yz+shape[2]), rate=c(Ez+rate[2]))
RR.post.area <- (y+shape[1])/(E+rate[1])
post.map <- process_MCMC_sample(post.sample, param.post.zone, RR.post.area, 
                                overlap$cluster.list, cutoffs)

# Posterior probs of k cluster/anti-clusters for k=0,...,J
k.vector <- table(sapply(post.chain$sample,length))
k.names <- as.numeric(names(k.vector))
k.vector <- normalize(k.vector)
pk.y <- rep(0, J+1)
pk.y[k.names+1] <- k.vector
print(paste("Posterior estimation complete on:", date()))



#-------------------------------------------------------------------------------
# Output
#-------------------------------------------------------------------------------
return(list(
  prior.map=prior.map, 
  post.map=post.map, 
  pk.y=pk.y))
}
