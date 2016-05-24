# # C++ code to update weights and pi's
# update_weights <- function(sums, alpha, nc) {
#   .Call('sourceR_update_weights', PACKAGE = 'sourceR', sums, alpha, nc)
# }
# 
# # C++ code to update clusters
# update_probs <- function(no_J, no_I, gmax, no_T, no_L, r, a, prev, theta, human_data, pi) {
#   .Call('sourceR_update_probs', PACKAGE = 'sourceR', no_J, no_I, gmax, no_T, no_L, r, a, prev, theta, human_data, pi)
# }
# 
# ## C++ code to calculate lambda i #################################################################### 
# calc_li <- function(no_J, no_I, no_T, no_L, r, a, prev, q) {
#   .Call('sourceR_calc_li', PACKAGE = 'sourceR', no_J, no_I, no_T, no_L, r, a, prev, q)
# }
# 
# set.seed(15653)
# no_T <- 1
# no_L <- 1
# no_J <- 5
# no_I <- 50
# gmax <- 20
# 
# ##########################################################################
# ############################# Simulate Data ##############################
# ##########################################################################
# a_fixed <- list()
# for(t in 1:times){
#   a_fixed[[t]] <- list()
#   for(l in 1:locations){
#     a_fixed[[t]][[l]] <- rexp(sources, 0.01)*4
#   }
# }
# 
# alpha_fixed <- 1
# 
# # This gives a Gamma base distribution for the Dirichlet Process governing 
# # the type effects
# q_prior <- c(1,1)
# 
# # Generate Poisson(lambda_i) distributed human cases for each type i
# # Model parameters are fixed for the source effects (a), the strength parameter for the DP 
# # (alpha) and the dispersion parameter (d)
# # The type effects are not fixed
# priors <- list(q=q_prior)
# fixed <- list(a=a_fixed, alpha=alpha_fixed)
# 
# # Generate Poisson(lambda_i) distributed human cases for each type i
# # Model parameters are fixed for the source effects (a), the strength parameter for the DP 
# # (alpha) and the dispersion parameter (d)
# # The type effects are not fixed
# priors <- list(q=q_prior, a=a_prior)
# fixed <- list(alpha=alpha_fixed)
# sim_data <- sim_data(times=no_T, locations=no_L, types=no_I, sources=no_J, 
#                      likelihood_dist="pois", params_prior=priors, 
#                      params_fixed=fixed)
# 
# real_data_sums <- apply(sim_data$data[,-1], 1, sum)
# real_data <- sim_data$data[which(real_data_sums != 0),]
# 
# r <- list()
# prev <- list()
# human_data <- list()
# for(t in 1:no_T){
#   human_data[[t]] <- list()
#   prev[[t]] <- rep(1,no_J)
#   r[[t]] <- apply(real_data[,-1], 2, function(x) x/sum(x))
#   for(l in 1:no_L){
#     human_data[[t]][[l]] <- real_data[,1]
#   }
# }
# q <- sim_data$params$q
# a <- sim_data$params$a
# alpha <- 1
# theta <- c(unique(q), rep(1, gmax-length(unique(q))))
# 
# plot(log(calc_li(no_J, no_I, no_T, no_L, r, a, prev, q)[[t]][[l]]), log(human_data[[t]][[l]]))
#   
# weight <- numeric(gmax)
# pis <- rep(1 / gmax, gmax)
# 
# weight[1] <- rbeta(1, 1, alpha)
# sticklen <- 1 - weight[1]
# for (k in 2 : (gmax - 1)) {
#   weight[k] <- rbeta(1, 1, alpha) # weight
#   pis[k] <- weight[k] * sticklen
#   sticklen <- sticklen - pis[k]
# }
# pis[gmax] <- sticklen
# 
# # prob test
# prob <- update_probs(no_J = no_J, no_I = no_I, gmax = gmax, no_T = no_T, no_L = no_L, r = r, a = a, prev = prev, theta = theta, human_data = human_data, pi = pis)
# 
# cluster <- sapply(1 : no_I, function(i) {probs <- prob[,i]; sample(as.factor(1 : gmax), 1, prob = probs)})
# nc <- table(cluster)
# sums <- sum(nc[2 : gmax])
# 
# # Weight test
# weight_pi <- update_weights(sums = sums, alpha = alpha, nc = as.vector(nc))
# 
# weight <<- weight_pi[, 1]
# pis <<- weight_pi[, 2]
# 
# # Update each cluster parameter
# params <- function(a1, b1, k) {
#   yk <- 0
#   sums <- 0
#   for (t in 1 : no_T) {
#     for (l in 1 : no_L) {
#       yk <- yk + sum(human_data[[t]][[l]][cluster == k])
#       sums <- sums + sum((r[[t]] %*% (a[[t]][[l]] * prev[[t]]))[which(cluster == k)])
#     }
#   }
#   a_star <- a1 + yk
#   b_star <- b1 + sums
#   return(c(a_star, b_star))
# }
# theta <- sapply(1 : gmax, function(k) {params <- params(a1 = q_prior[1], b1 = q_prior[2], k); rgamma(1, params[1], params[2])})
# q <- theta[cluster]
#          
# 
# 
# theta
# q
# weight
# pis
# cluster
