## ----setup, include=FALSE, cache=FALSE--------------------------------------------------
library(knitr)
# set global chunk options
opts_chunk$set(fig.path='figure/', fig.align='center', fig.show='hold')
options(replace.assign=TRUE, width=90)
opts_knit$set(out.format = "latex")
opts_chunk$set(warning=FALSE, message=FALSE, comment=NA, tidy=FALSE,
refresh=TRUE, cache=FALSE, autodep=TRUE)

## ----my-load-package, include=FALSE-----------------------------------------------------
library(metafolio)
set.seed(1)

## ----load-package, eval=FALSE-----------------------------------------------------------
#  library(metafolio)

## ----help, eval=FALSE-------------------------------------------------------------------
#  vignette("metafolio")
#  ?metafolio
#  help(package = "metafolio")

## ----setup-arma-------------------------------------------------------------------------
arma_env_params <- list(mean_value = 16, ar = 0.1, sigma_env = 2,
  ma = 0)

## ----base-case-run, cache=FALSE---------------------------------------------------------
base1 <- meta_sim(n_pop = 10, env_params = arma_env_params,
  env_type = "arma", assess_freq = 5)

## ----plot-base-case-ts, fig.width=5, fig.height=6.5, out.width="4in", fig.cap="An example simulation with stationary environmental stochasticity and the base-case parameter values."----
plot_sim_ts(base1, years_to_show = 70, burn = 30)

## ----spatial-plans----------------------------------------------------------------------
w_plans <- list()
w_plans[["balanced"]] <- c(5, 1000, 5, 1000, 5, 5, 1000, 5,
    1000, 5)
w_plans[["one_half"]] <- c(rep(1000, 4), rep(5, 6))
w <- list()
for(i in 1:2) { # loop over plans
 w[[i]] <- list()
 for(j in 1:80) { # loop over iterations
   w[[i]][[j]] <- matrix(w_plans[[i]], nrow = 1)
 }
}

## ----spatial-runs, fig.width=5, fig.height=5, out.width="4in", cache=FALSE, fig.cap="Two spatial conservation strategies shown in risk-return space with stationary environmental stochasticity. The dots show simulated metapopulations and the contours show 25\\% and 75\\% quantiles across 80 simulations per scenario. The grey line indicates the efficient frontier across all simulated metapopulations. The efficient frontier represents the minimum expected mean growth rate for a given expected variance in growth rate."----
set.seed(1)
arma_sp <- run_cons_plans(w, env_type = "arma", env_params =
  arma_env_params)
plot_cons_plans(arma_sp$plans_mv,
  plans_name = c("Balanced", "One half"),
  cols = c("#E41A1C", "#377EB8"), xlab = "Variance of growth rate",
  ylab = "Mean growth rate")

## ----env-eg, fig.width=5, fig.height=6, out.width="3in", fig.cap="Example environmental time series.", fig.pos="hpb"----
types <- c("sine", "arma", "regime", "linear", "constant")
x <- list()
for(i in 1:5) x[[i]] <- generate_env_ts(n_t = 100, type = types[i])
par(mfrow = c(5, 1), mar = c(3,3,1,0), cex = 0.7)
for(i in 1:5) plot(x[[i]], type = "o", main = types[i])

## ----plot-rickers, fig.width=7, fig.height=4.5, out.width="5in", fig.pos="hpb", fig.cap="Ricker curves from a sample of 40 years in the example simulation. Each panel represents a different stream population. Population 1 is more productive in cool conditions and population 10 is more producitive in warm conditions. The colour of the Ricker curves represents the relative temperatue in that year (warm: red; cool: blue). The grey shaded area represents the variation in spawners observed within the simulation."----
plot_rickers(base1, pal = rep("black", 10))

## ----return-correlations, fig.width=5, fig.height=5, out.width="5in", fig.cap="A plot comparing the log(returns) between each population. The population IDs are coloured from warm tolerant (warm colours) to cool tolerant (cool colours). Note how populations 1 and 10 have asynchronous returns whereas populations with more similar thermal-tolerance curves (say populations 9 and 10) have more synchronous dynamics. Populations with thermal tolerance curves in the middle (e.g.~population 6) are less correlated with other populations. Their population dynamics end up primarily driven by demographic stochasticity and less so by temperature-induced systematic changes in productivity."----
plot_correlation_between_returns(base1)


## ----monte-carlo-eg, eval=FALSE---------------------------------------------------------
#  set.seed(1)
#  weights_matrix <- create_asset_weights(n_pop = 6, n_sims = 3000,
#    weight_lower_limit = 0.001)
#  mc_ports <- monte_carlo_portfolios(weights_matrix = weights_matrix,
#    n_sims = 3000, mean_b = 1000)

## ----monte-carlo-eg-load, eval=TRUE, echo=FALSE-----------------------------------------
# To make the vignette compile more quickly:
# port_vals <- mc_ports$port_vals
# save(port_vals, file = "port_vals.rda")
weights_matrix <- create_asset_weights(n_pop = 6, n_sims = 3000,
  weight_lower_limit = 0.001)
load("port_vals.rda")
mc_ports <- list()
mc_ports$port_vals <- port_vals

## ----monte-carlo-plot, fig.width=7, fig.height=4, out.width="5in", fig.cap="Efficient frontier of metapopulation portfolios (red dots). Each dot represents a different set of weights of the Ricker $b$ parameters. The colours in the right panel correspond to the five populations with warm tolerant populations in warmer colours and cool tolerant populations in cooler colours."----
col_pal <- rev(gg_color_hue(6))
ef_dat <- plot_efficient_portfolios(port_vals = mc_ports$port_vals,
  pal = col_pal, weights_matrix = weights_matrix)

