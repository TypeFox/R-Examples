# Show difference of increasing the number of streams you maintain
# with unknown response diversity
#
# In this version, the response diversity is randomly drawn

set.seed(1)
USE_CACHE <- FALSE

# moved this out of the package to avoid Teaching Demos import
add_inset_env <- function(env, x = 0.12, y = -0.013, size = c(1, .5), ...) {
  TeachingDemos::subplot(plot(seq(0.1, 0.15, length.out =
        length(env)), env, axes = FALSE, ann = FALSE, type =
      "l", col = "grey10"), x = x, y = y, size = size, ...)
  par(xpd = NA)
  text(x, y, "Environment", adj = c(.5, 3.0), col = "black", cex = 1)
  par(xpd = FALSE)
}

# in this version we start with a set amount of b and can split it up among many
# or invest it in a few
n_trials <- 500 # number of trials at each n conservation plan
num_pops <- c(2, 4, 8, 16) # n pops to conserve
b_conserve <- 2000 / num_pops
n_plans <- length(num_pops) # number of plans
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

## ARMA:
arma_env_params <- list(mean_value = 16, ar = 0.1, sigma_env = 2, ma = 0)

pdf_eps("n-arma-sim-2", width = 5, height = 7, type = TYPE)
eg_arma <- meta_sim(b = w[[1]][[2]], n_pop = 16, env_params =
  arma_env_params, env_type = "arma", assess_freq = 5, max_a =
  thermal_integration(16))
plot_sim_ts(eg_arma, years_to_show = 70, burn = 30)
dev.off()

pdf_eps("n-arma-sim-16", width = 5, height = 7, type = TYPE)
eg_arma <- meta_sim(b = w[[4]][[2]], n_pop = 16, env_params =
  arma_env_params, env_type = "arma", assess_freq = 5, max_a =
  thermal_integration(16))
plot_sim_ts(eg_arma, years_to_show = 70, burn = 30)
dev.off()

if(!USE_CACHE) {
  x_arma_n <- run_cons_plans(w, env_type = "arma", env_params =
    arma_env_params, max_a = thermal_integration(16))
  x_arma_n$plans_port <- NULL # so save space
  save(x_arma_n, file = "x_arma_n.rda")
} else {
  load("x_arma_n.rda")
# has "plans_mv" only
}

## Linear:
linear_env_params <- list(min_value = 15, max_value = 19, sigma_env = 0.001,
  start_t = 30)

pdf_eps("n-linear-sim-2", width = 5, height = 7, type = TYPE)
eg_linear <- meta_sim(b = w[[1]][[1]], n_pop = 16, env_params =
  linear_env_params, env_type = "linear", assess_freq = 5, max_a =
  thermal_integration(16))
plot_sim_ts(eg_linear, years_to_show = 70, burn = 30)
dev.off()

pdf_eps("n-linear-sim-16", width = 5, height = 7, type = TYPE)
eg_linear <- meta_sim(b = w[[4]][[1]], n_pop = 16, env_params =
  linear_env_params, env_type = "linear", assess_freq = 5, max_a =
  thermal_integration(16))
plot_sim_ts(eg_linear, years_to_show = 100, burn = 30)
dev.off()

if(!USE_CACHE) {
  x_linear_n <- run_cons_plans(w, env_type = "linear",
    env_params = linear_env_params, max_a = thermal_integration(16))
  x_linear_n$plans_port <- NULL
  save(x_linear_n, file = "x_linear_n.rda")
} else {
  load("x_linear_n.rda")
# has "plans_mv" only
}

#cols <- RColorBrewer::brewer.pal(5, "Spectral")
cols <- RColorBrewer::brewer.pal(5, "Greys")[c(2:5)]
set.seed(2)
pdf_eps("cons-plans-n", width = 6.5, height = 6.8, type = TYPE)

layout(rbind(
  c(1, 1, 1, 2, 2, 2),
  c(1, 1, 1, 2, 2, 2),
  c(1, 1, 1, 2, 2, 2),
  c(1, 1, 1, 2, 2, 2),
  c(5, 5, 5, 5, 5, 5), # padding
  c(6, 3, 3, 3, 3, 7),
  c(6, 3, 3, 3, 3, 7),
  c(6, 4, 4, 4, 4, 7),
  c(6, 4, 4, 4, 4, 7)))

#1/9 padding
#2/9 botton
#4/9 top

xlim <- c(0.008, 0.90)
ylim <- c(-0.034, 0.027)
#par(family = "Times")
par(las = 1, cex = 0.8, mar = c(0, 0, 0, 0), oma = c(4, 5.2, 1.5, .5),
  tck = -0.02, mgp = c(2, .6, 0))
plot_cons_plans(x_arma_n$plans_mv, plans_name = plans_name_n, cols = cols,
  add_all_efs = FALSE, xlim = xlim, ylim = ylim, add_legend = FALSE)

add_inset_env(eg_arma$env_ts[-c(1:30)], x = 0.7, y = -0.026, size = c(1, .5))

mtext("(a) Short-term environmental fluctuations", side = 3, line = 0.2,
  cex = 0.8, adj = 0.05)

par(las = 0)
mtext("Mean of metapopulation growth rate", side = 2, line = 3,
  outer = FALSE, cex = 0.8)
par(las = 1)
plot_cons_plans(x_linear_n$plans_mv, plans_name = plans_name_n, cols = cols,
  add_all_efs = FALSE, xlim = xlim, ylim = ylim, y_axis = FALSE, add_legend = TRUE)

add_inset_env(eg_linear$env_ts[-c(1:30)], x = 0.7, y = -0.026, size = c(1, .5))

mtext("(b) Long-term environmental change", side = 3, line = 0.2, cex = 0.8,
  adj = 0.05)
mtext("Variance of metapopulation growth rate", side = 1, line = 2.25,
  outer = FALSE, cex = 0.8, adj = -3)

## time series plots:
par(tck = -0.035)
w <- list()
for(i in 1:n_plans) { # loop over number conserved
 w[[i]] <- list()
   w[[i]] <- matrix(rep(b_conserve[i], 16), nrow = 1)
}

   w[[1]][-c(8, 9)] <- 5 # conserve n, wipe out the rest
   w[[2]][-c(7:10)] <- 5 # conserve n, wipe out the rest
   w[[3]][-c(5:12)] <- 5 # conserve n, wipe out the rest
   w[[4]][-c(1:16)] <- 5 # conserve n, wipe out the rest

cons_arma_ts <- list()
set.seed(18182)
for(i in 1:4) {
  use_cache <- ifelse(i == 1, FALSE, TRUE)
  cons_arma_ts[[i]] <- meta_sim(b = w[[i]], n_pop =
    ncol(w[[i]]), env_params = arma_env_params, env_type =
    "arma", assess_freq = 5, use_cache = use_cache, cache_env =
    use_cache, skip_saving_cache = FALSE)
}
cons_linear_ts <- list()
for(i in 1:4) {
  use_cache <- ifelse(i == 1, FALSE, TRUE)
  cons_linear_ts[[i]] <- meta_sim(b = w[[i]], n_pop =
    ncol(w[[i]]), env_params = linear_env_params, env_type =
    "linear", assess_freq = 5, use_cache = use_cache, cache_env =
    use_cache, skip_saving_cache = FALSE)
}

plot_sp_A_ts(list(cons_arma_ts[[1]], cons_arma_ts[[4]]), ylim = c(-1.2, 1.2), rate = TRUE, x_axis = FALSE,
  labels = "(c)\n", cols = cols[c(1, 4)])

par(las = 0)
mtext("Metapopulation\ngrowth rate", side = 2, line = 3, outer = FALSE, cex = 0.8)
par(las =1)

par(xpd = NA)
text(72, 0.5, "2 populations", pos = 4, col = cols[2])
text(72, 0.1, "16 populations", pos = 4, col = cols[4])
par(xpd = FALSE)

plot_sp_A_ts(list(cons_arma_ts[[1]], cons_arma_ts[[4]]), ylim = c(500, 6000), rate = FALSE, log = "y",
  y_axis = TRUE, y_axis_ticks = c(1000, 2000, 5000),
  labels = "(d)", cols = cols[c(1, 4)])
mtext("(d)", side = 3, line = -1.2, cex = 0.8, adj = 0.025)

par(xpd = NA)
mtext("Generation", side = 1, line = 2, outer = FALSE, cex = 0.8)
par(xpd = FALSE)

par(las = 0)
mtext("Metapopulation\nabundance", side = 2, line = 3, outer = FALSE, cex = 0.8)
par(las =1)

dev.off()

# report summary statistics:
# short term summaries:

mean.v <- plyr::ldply(x_arma_n$plans_mv, function(x) mean(x$v))

message("mean variance in growth rate for short term: 4 vs. 16")
message(round(mean.v$V1[2]/mean.v$V1[4], 1))

message("mean variance in growth rate for short term: 8 vs. 16")
message(round(mean.v$V1[3]/mean.v$V1[4], 1))

quant.v <- plyr::laply(x_arma_n$plans_mv, function(x) {
  q.l <- as.numeric(quantile(x$v, probs = 0.125))
  q.u <- as.numeric(quantile(x$v, probs = 0.875))
  q.u - q.l
    }
  )
message("width of 75% quantile of variance for short term: 16 vs. 4")
message(round(quant.v[2]/quant.v[4], 1))

message("width of 75% quantile of variance for short term: 16 vs. 8")
message(round(quant.v[3]/quant.v[4], 1))

# long term summaries:

quant.m <- plyr::laply(x_linear_n$plans_mv, function(x) {
  q.l <- as.numeric(quantile(x$m, probs = 0.125))
  q.u <- as.numeric(quantile(x$m, probs = 0.875))
  q.u - q.l
    }
  )
message("width of 75% quantile of mean for long term: 16 vs. 4")
message(round(quant.m[2]/quant.m[4], 1))

message("width of 75% quantile of mean for long term: 16 vs. 8")
message(round(quant.m[3]/quant.m[4], 1))

#sink("cons_plans_n_random_resp_div.tex")
#n_arma_75q_var_width_16_v_8 <- round(quant.m[3]/quant.m[4], 1)
#cat(paste0("\\newcommand{\\narma75qvarwidth16vs8}{", sprintf("%.1f", round(quant.m[3]/quant.m[4], 1))), "}")
#sink()

message("mean variance in growth rate for long term")
mean.v <- plyr::ldply(x_linear_n$plans_mv, function(x) mean(x$v))
message("16 vs. 4")
message(round(mean.v$V1[2]/mean.v$V1[4], 1))
message("16 vs. 2")
message(round(mean.v$V1[1]/mean.v$V1[4], 1))
message("16 vs. 8")
message(round(mean.v$V1[3]/mean.v$V1[4], 1))
message("8 vs. 4")
message(round(mean.v$V1[2]/mean.v$V1[3], 1))
