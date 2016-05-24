# library(ECctmc)
# library(GillespieSSA)
# library(reshape2)
# library(ggplot2)
# set.seed(52787)
#
# # Generate sample paths for SIS dynamics with GillespieSSA ---------------------------------
# obstimes <- seq(0,10,by=0.05); niter <- 500000
# params <- c(mu = rnorm(1, 0.5, 1e-3), beta = rnorm(1, 0.5, 1e-3))
# a = c("beta * S", "mu * I")
# nu <- matrix(c(-1, 1, 1, -1),nrow=2)
# init_state = c(S = 1, I = 0)
#
# SSA_paths <- matrix(0, nrow = length(obstimes), ncol = niter)
#
# k <- 1
#
# while(k <= niter) {
#         if(k %% 1000 == 0) print(k)
#
#         # simulate paths
#         path <- ssa(init_state, a, nu, params, tf = 10)$data[,c(1,3)]
#         path <- path[-nrow(path),]
#         path <- path[path[,1] < 10,, drop = FALSE]
#         path <- rbind(path, c(10, path[nrow(path), 2]))
#
#         # accept the paths that end in infection
#         if(path[nrow(path),2] == 1) {
#                 path_long <- path[findInterval(obstimes, path[,1]), 2]
#                 SSA_paths[, k] <- path_long
#                 k <- k+1
#         }
# }
#
# # Generate sample paths using ECctmc ------------------
#
# Q <- matrix(c(-params[1], params[1], params[2], -params[2]), nrow = 2, byrow = T)
# MR_paths <- matrix(0, nrow = length(obstimes), ncol = niter)
# Unif_paths <- matrix(0, nrow = length(obstimes), ncol = niter)
#
# t0_mr_start <- Sys.time()
# for(k in 1:niter) {
#
#         # simulate paths
#         path <- sample_path_mr(1, 2, 0, 10, Q)
#         path_long <- path[findInterval(obstimes, path[,1]), 2] - 1
#         MR_paths[,k] <- path_long
#
# }
# t0_mr_end <- Sys.time()
#
# t0_unif_start <- Sys.time()
# for(k in 1:niter) {
#
#         path <- sample_path_unif(1, 2, 0, 10, Q)
#         path_long <- path[findInterval(obstimes, path[,1]), 2] - 1
#         Unif_paths[,k] <- path_long
# }
# t0_unif_end <- Sys.time()
#
# # Collect results ---------------------------------------------------------
#
# SSA_res <- cbind(time = obstimes, I_mean = rowMeans(SSA_paths),
#                  lower = rowMeans(SSA_paths) - 1.96 * apply(SSA_paths, 1, sd)/sqrt(niter),
#                  upper = rowMeans(SSA_paths) + 1.96 * apply(SSA_paths, 1, sd)/sqrt(niter))
# MR_res <- cbind(time = obstimes, I_mean = rowMeans(MR_paths),
#                  lower = rowMeans(MR_paths) - 1.96 * apply(MR_paths, 1, sd)/sqrt(niter),
#                  upper = rowMeans(MR_paths) + 1.96 * apply(MR_paths, 1, sd)/sqrt(niter))
# Unif_res <- cbind(time = obstimes, I_mean = rowMeans(Unif_paths),
#                  lower = rowMeans(Unif_paths) - 1.96 * apply(Unif_paths, 1, sd)/sqrt(niter),
#                  upper = rowMeans(Unif_paths) + 1.96 * apply(Unif_paths, 1, sd)/sqrt(niter))
#
# res <- data.frame(rbind(SSA_res, MR_res, Unif_res)); res <- cbind(method = rep(c("SSA", "MR", "Unif"), each = length(obstimes)), res)
# res_melt <- melt(res, id.vars = c("time", "method"))
#
# ggplot(res, aes(x = as.numeric(time)))+ geom_line(aes(y = I_mean, colour = method)) + geom_ribbon(aes(ymin = lower, ymax = upper, fill = method), alpha = 0.5) + labs(x = "time", y = "Pr(Xt=2") + theme_bw()
