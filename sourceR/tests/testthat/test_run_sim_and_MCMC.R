# library(sourceR)
# 
# # set.seed(19127)
# set.seed(15653)
# times <- 2
# locations <- 2
# sources <- 5
# types <- 50
# 
# ##########################################################################
# ############################# Simulate Data ##############################
# ##########################################################################
# 
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
# sim_data <- sim_data(times=times, locations=locations, types=types, sources=sources, 
#                      likelihood_dist="pois", params_prior=priors, 
#                      params_fixed=fixed)
# 
# sim_data$params$a
# unique(sim_data$params$q)
# sim_data$data[,1]
# sim_data$lambda_j
# 
# for (t in 1 : times) {
#   for (l in 1 : locations) {
#     print(sim_data$lambda_j[[t]][[l]] / sum(sim_data$lambda_j[[t]][[l]]))
#   }
# }
# 
# 
# 
# # times <- 1
# # locations <- 1
# # sources <- 6#5
# # types <- 91 #183
# # # real_data <- read.csv(file="/media/poppy/SP UFD U2/Honours Project/MyCode/NewCodeWeek10/Real Data/FullData_other.csv", row.names=1) # full data
# # real_data <- read.csv(file="/media/poppy/SP UFD U2/sourceR_jss_paper/Data_case_study/real_data.csv", row.names=1) # data from 2010 petra paper
# # real_data_sums <- apply(real_data[,-1], 1, sum)
# # 
# # real_data <- real_data[which(real_data_sums != 0),]
# 
# ##########################################################################
# #### Initial values ######################################################
# ##########################################################################
# 
# a1 <- list()
# for(t in 1:times){
#   a1[[t]]<- list()
#   for(l in 1:locations){
# #     a1[[t]][[l]] <- sim_data$params$a[[t]][[l]] #rexp(sources, 0.01)*3
#     a1[[t]][[l]] <- rexp(sources, 0.01)*3
#   }
# }
# 
# q1 <- c(rep(1, types))
# 
# initials <- list(a=a1, q=q1, alpha=1)
# # initials <- list(q=q1, alpha=1)
# 
# ##########################################################################
# #### Tuning ##############################################################
# ##########################################################################
# a_cent <- list()
# a_noncent <- list()
# for(t in 1:times){
#   a_cent[[t]] <- list()
#   a_noncent[[t]] <- list()
#   for(l in 1:locations){
#     a_cent[[t]][[l]] <- 0.4
#     a_noncent[[t]][[l]] <- 0.4 # usually larger than the noncentered
#   }
# }
# #OR
# tunes <- list(a_cent=a_cent, a_noncent=a_noncent, r=0.01)
# 
# ##########################################################################
# #### Priors ##############################################################
# ##########################################################################
# a_prior <- list()
# r_prior <- list()
# for(t in 1:times){
#   a_prior[[t]] <- list()
#   r_prior[[t]] <- 0.01
#   for(l in 1:locations){
#     a_prior[[t]][[l]] <- rep(0.01, sources) # Exp(0.01) prior
#   }
# }
# priors <- list(a=a_prior, r=r_prior, q=c(1,1), alpha=c(1,1))
# 
# ##########################################################################
# #### Run model ###########################################################
# ##########################################################################
# n.iter=1000
# res <- list()
# ptime <- system.time({
# res[[1]] <- doMCMC(raw_data=sim_data$data,
# #               raw_data=real_data, 
#               params_init=initials, 
#               params_tune=tunes, 
#               params_prior=priors, 
#               n_iter=n.iter, 
#               likelihood_dist="pois", 
#               sources=sources,
#               types=types, 
#               times=times, 
#               locations=locations,
#               gmax=25,
#               fix_a1=FALSE,
#               fix_rij=FALSE,
#               fix_alpha=FALSE,
#               fix_type_effects=FALSE,
#               save_lambda=TRUE,
#               n_rij=200)
# })
# print(ptime)
# # output(BI=1, thin=1, res_data=res, s_data=sim_data, prior=priors, prev=c(rep(1,sources)))# output(BI=100, thin=1, res_data=res, prior=priors, prev=prev)
# output(BI=1000, thin=100, res_data=res, prior=priors, prev=c(rep(1,sources)))
# 
# lambda_j_proportion <- list()
# for(t in 1:times){
#   lambda_j_proportion[[t]] <- list()
#   for(l in 1:locations){
#     lambda_j_proportion[[t]][[l]] <- t(apply(X=res$posterior$lj[[t]][[l]], MARGIN=1:1, FUN=function(x) x/sum(x)))
# #     lambda_j_proportion[[t]][[l]] <- as.ff(t(apply(X=res$posterior$lj[[l]][[l]][], MARGIN=1:1, FUN=function(x) x/sum(x))))
#   }
# }
# 
# X11()
# par(mfrow=c(3,3))
# for(t in 1:times){
#   for(l in 1:locations){
#     for(j in 1:sources){
# #       plot(lambda_j_proportion[[t]][[l]][,j], type="l")
#       plot(density(lambda_j_proportion[[t]][[l]][,j]), xlim=c(0,1), main=paste("t",t, "l",l, "j", j))
#       abline(v=sim_data$lambda_j[[t]][[l]][j]/sum(sim_data$lambda_j[[t]][[l]]),col="red")
# #       plot(density(res$posterior$lj[[t]][[l]][,j]))
# #       abline(v=sim_data$lambda_j[[t]][[l]][j],col="red")
#     }
#   }
# }
# 
# plot(density(res$posterior$alpha), main="alpha")
# # plot(res$posterior$alpha[], type="l", main="alpha")
# 
# X11()
# par(mfrow=c(5,5))
# for(i in 1:types){
#   plot(density(res$posterior$q[,i]))
# }
# 
# sums <- summary(res, alpha=0.05, burnin=1000, thin=1)
# 
# X11()
# plot(log(sums$li[[1]][[1]][,'median']), log(sim_data$data[((t - 1) * no_I + 1) : (t * no_I),1]))
# abline(0,1)
# for (t in 1 : times) {
#   for (l in 1 : locations) {
#     points(log(sums$li[[t]][[l]][,'median']), log(sim_data$data[((t - 1) * no_I + 1) : (t * no_I),l]))
#   }
# }
# 
# 
# 
