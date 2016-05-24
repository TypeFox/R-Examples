# Here, examine a scenario where you start with a fixed quantity of `b` to
# invest, different starting values of `b`, and `b` declines in absolute
# quantities across the simulations.
#
# Strategies to test:
# (1) Invest evenly in all streams
# (2) Triage: pick those that start above a threshold and distribute
# `b` among them. Put a tiny bit of `b` (5) in the rest.
# (3) Panicker: pick those that start below a threshold and concentrate on
# them, put some `b` (a bit above extinction by the end) in the bigger ones.
#
# Have `b` decline continuously by some constant magnitude across stocks. E.g.
# have `b` decline by 200 units over 100 years.
#
# is is respon. diversity, initial capacity a better preditor...
#
# simpler idea: start with X quantity of b and you can distribute it however...
# you can either distribute your b as follows:
# in 2, 4, 8, or 16 populations with the b divided and 5 in the rest
# so in the fewer Ns, they will all go to extinction or near extinction
# at higher Ns, they will stay far away from extinction

set.seed(1)
USE_CACHE <- FALSE

n_trials <- 500 # number of trials at each n conservation plan
num_pops <- c(2, 4, 8, 12, 16) # n pops to conserve
b_conserve <- 2000 / num_pops
n_plans <- length(num_pops)
w <- list()
for(i in 1:n_plans) { # loop over number conserved
 w[[i]] <- list()
 for(j in 1:n_trials) { # loop over trials
   w[[i]][[j]] <- matrix(rep(b_conserve[i], 16), nrow = 1)
   # conserve num_pops[i] populations; wipe out rest:
   w[[i]][[j]][-sample(1:16, num_pops[i])] <- 5
 }
}
plans_name_n <- paste(num_pops, "populations")
cols <- RColorBrewer::brewer.pal(6, "Greys")[c(2:6)]

##########
w <- list()
for(i in 1:n_plans) { # loop over number conserved
 w[[i]] <- list()
 for(j in 1:n_trials) { # loop over trials
   w[[i]][[j]] <- matrix(rep(b_conserve[i], 16), nrow = 1)
   # conserve num_pops[i] populations; wipe out rest:
   w[[i]][[j]][-sample(1:16, num_pops[i])] <- 5
 }
}
plans_name_n <- paste(num_pops, "populations")

# linear-arma version:

linear_arma_env_params <- list(min_value = 15, max_value = 19,
  start_t = 30, mean_value = 16, ar = 0.1, sigma_env = 2, ma = 0)

set.seed(123)
pdf_eps("n-linear-arma-sim-16-squeeze", width = 5, height = 7, type = TYPE)
# try a version with ARMA and linear change:
eg_linear_arma <- meta_sim(b = w[[5]][[2]], n_pop = 16, env_params =
 linear_arma_env_params, env_type = "linear_arma",
  assess_freq = 5, max_a = thermal_integration(16), decrease_b = 0.85)
plot_sim_ts(eg_linear_arma, years_to_show = 100, burn = 30)
dev.off()

pdf_eps("n-linear-arma-sim-2-squeeze", width = 5, height = 7, type = TYPE)
# try a version with ARMA and linear change:
eg_linear_arma <- meta_sim(b = w[[1]][[2]], n_pop = 16, env_params =
 linear_arma_env_params, env_type = "linear_arma",
  assess_freq = 5, max_a = thermal_integration(16), decrease_b = 0.85)
plot_sim_ts(eg_linear_arma, years_to_show = 100, burn = 30)
dev.off()

if(!USE_CACHE) {
  x_linear_arma_n <- run_cons_plans(w, env_type = "linear_arma", env_params =
    linear_arma_env_params, max_a = thermal_integration(16), decrease_b = 0.85)
  x_linear_arma_n$plans_port <- NULL
  save(x_linear_arma_n, file = "x_linear_arma_n.rda")
} else {
  load("x_linear_arma_n.rda")
# only has "plans_mv"
}

# some metapop time series panels to plot:
for(i in 1:n_plans) { # loop over number conserved
 w[[i]] <- list()
   w[[i]] <- matrix(rep(b_conserve[i], 16), nrow = 1)
}

w[[1]][-c(8, 9)] <- 5 # conserve 2, wipe out the rest
w[[2]][-c(7:10)] <- 5 # conserve 4, wipe out the rest
w[[3]][-c(5:12)] <- 5 # conserve 8, wipe out the rest
w[[4]][-c(3:14)] <- 5 # conserve 12, wipe out the rest
w[[5]][-c(1:16)] <- 5 # conserve 16, wipe out the rest

set.seed(1279)
cons_linear_arma_ts <- list()
for(i in 1:length(w)) {
  use_cache <- ifelse(i == 1, FALSE, TRUE)
  cons_linear_arma_ts[[i]] <- meta_sim(b = w[[i]], n_pop =
    ncol(w[[i]]), env_params = linear_arma_env_params, env_type =
    "linear_arma", assess_freq = 5, use_cache = use_cache, cache_env =
    use_cache, decrease_b = 0.85)
}

## figure:
pdf_eps("cons-plans-squeeze", width = 4.0, height = 7, type = TYPE)
layout(rbind(
  c(1),
  c(1),
  c(1),
  c(4),
  c(2),
  c(2),
  c(3),
  c(3)))

xlim <- c(0.08, 0.9)
ylim <- c(-0.038, 0.028)
#par(family = "Times")
par(las = 1, cex = 0.8, mar = c(0, 0, 0, 0), oma = c(4, 5.2, 1.8, .5),
  tck = -0.02, mgp = c(2, .5, 0))
plot_cons_plans(x_linear_arma_n$plans_mv, plans_name = plans_name_n, cols = cols,
  add_all_efs = FALSE, xlim = xlim, ylim = ylim, add_legend = TRUE, add_poly = TRUE,
  legend_pos = "bottomright")
mtext("(a) Reduction in stream flow", side = 3, line = .4,
  cex = 0.8, adj = 0.05)
mtext("Variance of metapopulation growth rate", side = 1, line = 2.25,
  outer = FALSE, cex = 0.8)
par(las = 0)
mtext("Mean of metapopulation growth rate", side = 2, line = 3,
  outer = FALSE, cex = 0.8)
par(las = 1)

# ticks need to be a bit bigger here to match:
par(tck = -0.03)
plot_sp_A_ts(list(cons_linear_arma_ts[[1]], cons_linear_arma_ts[[5]]), ylim = c(-1.25, 1.25), rate = TRUE, x_axis = FALSE, labels = "(b) \n", cols = cols[c(2, 5)], add_lm = FALSE)
par(las = 0)
mtext("Metapopulation\ngrowth rate", side = 2, line = 3, outer = FALSE, cex = 0.8)
par(las =1)

plot_sp_A_ts(list(cons_linear_arma_ts[[1]], cons_linear_arma_ts[[5]]), ylim = c(0, 5000), rate = FALSE, x_axis = TRUE, labels = "(c)\n", cols = cols[c(2, 5)], add_lm = FALSE)
par(las = 0)
mtext("Metapopulation\nabundance", side = 2, line = 3, outer = FALSE, cex = 0.8)
par(las =1)
par(xpd = NA)
mtext("Generation", side = 1, line = 2, outer = FALSE, cex = 0.8)

text(35, 400, "16 populations", pos = 4, col = cols[5])
text(11, 4400, "2 populations", pos = 4, col = cols[2])


par(xpd = FALSE)
dev.off()

# results to quantify:
# - reduction in mean variance with increasing N
# - but also reduction in mean growth rate with increasing N
# - note that this forms the efficient fontier
#
mean.v <- plyr::ldply(x_linear_arma_n$plans_mv, function(x) mean(x$v))
message("mean variance of 12 compared to 4")
message(round(mean.v$V1[2] / mean.v$V1[4], 1))

mean.m <- plyr::ldply(x_linear_arma_n$plans_mv, function(x) mean(x$m))
message("mean mean of 16 compared to 8")
message(round(mean.v$V1[3] / mean.v$V1[5], 1))
