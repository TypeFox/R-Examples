# This file looks at different spatial conservation strategies

set.seed(2)
USE_CACHE <- FALSE

w_plans <- list()
w_plans[[1]] <- c(5, 1000, 5, 1000, 5, 5, 1000, 5, 1000, 5)
w_plans[[2]] <- c(5, 5, 5, 1000, 1000, 1000, 1000, 5, 5, 5)
w_plans[[3]] <- c(rep(1000, 4), rep(5, 6))
w_plans[[4]] <- rev(w_plans[[3]])
plans_name_sp <- c("Full response range", "Most stable only",
  "Lower half", "Upper half")

n_trials <- 500 # number of trials at each n conservation plan
num_pops <- c(10, 10, 10, 10) # n pops to conserve
n_plans <- length(num_pops) # number of plans
w <- list()
for(i in 1:n_plans) { # loop over plans
 w[[i]] <- list()
 for(j in 1:n_trials) { # loop over trials
   w[[i]][[j]] <- matrix(w_plans[[i]], nrow = 1)
 }
}

## ARMA:
arma_env_params <- list(mean_value = 16, ar = 0.1, sigma_env = 3, ma = 0)

pdf_eps("spatial-arma-sim-full", width = 5, height = 7, type = TYPE)
#par(family = "Times")
set.seed(1)
eg_arma <- meta_sim(b = w[[1]][[1]], n_pop = 10, env_params = arma_env_params,
  env_type = "arma", assess_freq = 5)
col_pal <- grey(12:0/12)[3:12]
plot_sim_ts(eg_arma, years_to_show = 100, burn = 30, yticks = list(NA, NA, NA,
    NA, NA, c(0, 100), c(0, 20), c(-2, 0, 2), c(1, 2), c(0, 1000)), pal = col_pal, oma = c(4, 4.5, 1.7, 1))
mtext("Shaded lines represent individual populations", side = 3, line = 0.55, outer = TRUE, adj = 0.05, cex = 0.7)
dev.off()

pdf_eps("spatial-arma-sim-full-colour", width = 5, height = 7, type = TYPE)
#par(family = "Times")
set.seed(1)
eg_arma <- meta_sim(b = w[[1]][[1]], n_pop = 10, env_params = arma_env_params,
  env_type = "arma", assess_freq = 5)
plot_sim_ts(eg_arma, years_to_show = 100, burn = 30, yticks = list(NA, NA, NA,
    NA, NA, NA, NA, c(2, 0, 2), NA, NA))
dev.off()

eg_arma <- meta_sim(b = rep(1000, 10), n_pop = 10, env_params = arma_env_params,
  env_type = "arma", assess_freq = 5)

pdf_eps("example-return-correlations", width = 5, height = 5, type = TYPE)
plot_correlation_between_returns(eg_arma)
dev.off()

pdf_eps("spatial-arma-sim-onehalf", width = 5, height = 7, type = TYPE)
eg_arma <- meta_sim(b = w[[4]][[1]], n_pop = 10, env_params = arma_env_params,
  env_type = "arma", assess_freq = 5)
plot_sim_ts(eg_arma, years_to_show = 100, burn = 30)
dev.off()

if(!USE_CACHE) {
x_arma_sp <- run_cons_plans(w, env_type = "arma", env_params =
  arma_env_params)
x_arma_sp$plans_port <- NULL
  save(x_arma_sp, file = "x_arma_sp.rda")
} else {
  load("x_arma_sp.rda")
}

## Linear:
linear_env_params <- list(min_value = 15, max_value = 19, sigma_env = 0.001,
  start_t = 30)

pdf_eps("spatial-linear-sim-full", width = 5, height = 7, type = TYPE)
eg_linear <- meta_sim(b = w[[1]][[1]], n_pop = 10, env_params =
    linear_env_params, env_type = "linear", assess_freq = 5)
plot_sim_ts(eg_linear, years_to_show = 100, burn = 30)
dev.off()

pdf_eps("spatial-linear-sim-onehalf", width = 5, height = 7, type = TYPE)
eg_linear <- meta_sim(b = w[[3]][[1]], n_pop = 10, env_params =
    linear_env_params, env_type = "linear", assess_freq = 5)
plot_sim_ts(eg_linear, years_to_show = 100, burn = 30)
dev.off()

if(!USE_CACHE) {
  x_linear_sp <- run_cons_plans(w, env_type = "linear", env_params =
    linear_env_params, max_a = thermal_integration(10))
  x_linear_sp$plans_port <- NULL
  save(x_linear_sp, file = "x_linear_sp.rda")
} else {
  load("x_linear_sp.rda")
}

cols <- RColorBrewer::brewer.pal(5, "Dark2")
#cols <- RColorBrewer::brewer.pal(5, "Greys")[c(2:5)]

pdf_eps("spatial-mv", width = 6.5, height = 6.8, type = TYPE)
layout(rbind(
  c(1, 2),
  c(1, 2),
  c(1, 2),
  c(1, 2),
  c(7, 8), # padding
  c(3, 5),
  c(3, 5),
  c(4, 6),
  c(4, 6)))

xlim <- c(0.18, 0.80)
ylim <- c(-0.027, 0.027)
par(las = 1, cex = 0.8, mar = c(0, 0, 0, 0), oma = c(4, 5.2, 1.5, .5),
  tck = -0.02, mgp = c(2, .6, 0))
#par(family = "Times")
plot_cons_plans(x_arma_sp$plans_mv, plans_name = plans_name_sp, cols = cols,
  add_all_efs = FALSE, xlim = xlim, ylim = ylim, add_legend = FALSE)
#add_inset_env(eg_arma$env_ts[-c(1:30)], x = 0.12, y = -0.013, size = c(1, .5))

mtext("(a) Short-term environmental fluctuations", side = 3, line = 0.2, cex =
  0.8, adj = 0.05)
par(las = 0)
mtext("Mean of metapopulation growth rate", side = 2, line = 3, outer = FALSE,
  cex = 0.8)
par(las = 1)

plot_cons_plans(x_linear_sp$plans_mv, plans_name = plans_name_sp, cols = cols,
  add_all_efs = FALSE, xlim = xlim, ylim = ylim, y_axis = FALSE, add_legend = TRUE, legend_pos = "bottomright")
#add_inset_env(eg_linear$env_ts[-c(1:30)], x = 0.12, y = -0.013, size = c(1, .5))

mtext("(b) Long-term environmental change", side = 3, line = 0.2, cex = 0.8, adj = 0.05)
mtext("Variance of metapopulation growth rate", side = 1, line = 2.25, outer =
  FALSE, cex = 0.8, adj = -3)

# time series plots:

set.seed(4)
par(tck = -0.035)
cons_arma_ts <- list()
for(i in 1:4) {
  use_cache <- ifelse(i == 1, FALSE, TRUE)
  cons_arma_ts[[i]] <- meta_sim(b = w[[i]][[1]], n_pop = 10, env_params =
    arma_env_params, env_type = "arma", assess_freq = 5,
    use_cache = use_cache, sigma_v = 0.2, skip_saving_cache = FALSE)
}
cons_linear_ts <- list()
for(i in 1:4) {
  use_cache <- ifelse(i == 1, FALSE, TRUE)
  cons_linear_ts[[i]] <- meta_sim(b = w[[i]][[1]], n_pop = 10, env_params =
    linear_env_params, env_type = "linear", assess_freq = 5,
    use_cache = use_cache, skip_saving_cache = FALSE)
}
burn <- 1:30

plot_sp_A_ts(cons_arma_ts, ylim = c(0000, 12400),
  start_new_plots = c(1, 3),
  labels = c("(c) Response diversity dampens\n     short-term risk",
    "ignore", "(e)\n", "ignore"), cols = cols)

####### temp
# plot_sp_A_ts(cons_arma_ts, ylim = c(-2, 2),
#   start_new_plots = c(1, 3),
#   labels = c("(c) Response diversity dampens\n     short-term risk",
#     "ignore", "(e)\n", "ignore"), cols = cols, rate = TRUE)
# plot_sp_A_ts(cons_linear_ts, ylim = c(-2, 2), y_axis = FALSE,
#   start_new_plots = c(1, 3), labels =
#   c("(d) Response diversity ensures\n      long-term persistence",
#     "ignore", "(f)\n", "ignore"), cols = cols, rate = TRUE)
######

par(las = 0)
mtext("Metapopulation abundance", side = 2, line = 3, outer = FALSE, cex = 0.8, adj = -2)
par(las =1)

plot_sp_A_ts(cons_linear_ts, ylim = c(0000, 12400), y_axis = FALSE,
  start_new_plots = c(1, 3), labels =
  c("(d) Response diversity ensures\n      long-term persistence",
    "ignore", "(f)\n", "ignore"), cols = cols)

par(xpd = NA)
mtext("Generation", side = 1, line = 2, outer = FALSE, cex = 0.8, adj = -.2)
par(xpd = FALSE)
dev.off()

## report summary statistics:
mean.v <- plyr::ldply(x_arma_sp$plans_mv, function(x) mean(x$v))
message(round(mean(mean.v$V1[3:4]) / mean(mean.v$V1[1:2]), 1))

mean.m <- plyr::ldply(x_linear_sp$plans_mv, function(x) mean(x$m))
print(mean.m)
#message(round(mean(mean.m$V1[3:4]) / mean(mean.m$V1[1:2]), 1))

#############################
#
# set.seed(2)
# cons_linear_ts <- list()
# for(i in 1:4) {
#   use_cache <- ifelse(i == 1, FALSE, TRUE)
#   cons_linear_ts[[i]] <- meta_sim(b = w[[i]][[1]], n_pop = 10, env_params =
#     linear_env_params, env_type = "linear", assess_freq = 5,
#     use_cache = use_cache, skip_saving_cache = FALSE)
# }
# burn <- 1:30
# par(mfrow = c(2, 1), mar = c(2, 4, 1, 1))
# plot_sp_A_ts(cons_linear_ts, ylim = c(0000, 12400),
#   start_new_plots = c(1, 3), cols = cols)
