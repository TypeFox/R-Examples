arma_env_params <- list(mean_value = 16, ar = 0.1, sigma_env = 2, ma = 0)
n_pop <- 4
col_pal <- rev(gg_color_hue(n_pop - 1)) # throw one away to simplify
seed_val <- 123

eg <- list()
ex <- list()

set.seed(seed_val)
eg[["Base case"]] <- meta_sim(
  n_pop       = n_pop,
  env_params  = arma_env_params,
  env_type    = "arma",
  assess_freq = 5
  )
ex[["base"]] <- expression("Base case, see main text or Table S1 for parameter values")

set.seed(seed_val)
eg[["No temperature fluctuations"]] <- meta_sim(
  n_pop       = n_pop,
  env_params  = list(mean_value = 16, ar = 0.1, sigma_env = 0.00000001, ma = 0),
  env_type    = "arma",
  assess_freq = 5
  )
ex[["No temperature fluctuations"]] <- expression(No~temperature~fluctuations~sigma[d]==0)

set.seed(seed_val)
eg[["Larger temperature fluctuations"]] <- meta_sim(
  n_pop       = n_pop,
  env_params  = list(mean_value = 16, ar = 0.1, sigma_env = 5, ma = 0),
  env_type    = "arma",
  assess_freq = 5
  )
ex[["Larger temperature fluctuations"]] <- expression(Larger~temperature~fluctuations~sigma[d]==5)

set.seed(seed_val)
eg[["No straying"]] <- meta_sim(
  n_pop       = n_pop,
  env_params  = arma_env_params,
  env_type    = "arma",
  assess_freq = 5,
  stray_fraction = 0.00000001
  )
ex[["No straying"]] <- expression(No~straying~f[stray]==0)

set.seed(seed_val)
eg[["More straying"]] <- meta_sim(
  n_pop       = n_pop,
  env_params  = arma_env_params,
  env_type    = "arma",
  assess_freq = 5,
  stray_fraction = 0.2
  )
ex[["More straying"]] <- expression(More~straying~f[stray]==0.2)

set.seed(seed_val)
eg[["No spawner-return variability"]] <- meta_sim(
  n_pop       = n_pop,
  env_params  = arma_env_params,
  env_type    = "arma",
  assess_freq = 5,
  sigma_v = 0.000000001
  )
ex[["No spawner-return variability"]] <- expression(No~recruitment~variability~sigma[r]==0)

set.seed(seed_val)
eg[["Larger spawner-return variability"]] <- meta_sim(
  n_pop       = n_pop,
  env_params  = arma_env_params,
  env_type    = "arma",
  assess_freq = 5,
  sigma_v = 1.0
  )
ex[["Larger spawner-return variability"]] <- expression(Larger~recruitment~variability~sigma[r]==1.0)

set.seed(seed_val)
eg[["No implementation error"]] <- meta_sim(
  n_pop       = n_pop,
  env_params  = arma_env_params,
  env_type    = "arma",
  assess_freq = 5,
  sigma_impl = 0.0000001
  )
ex[["No implementation error"]] <- expression(No~implementation~error~sigma[h]==0)

set.seed(seed_val)
eg[["Larger implementation error"]] <- meta_sim(
  n_pop       = n_pop,
  env_params  = arma_env_params,
  env_type    = "arma",
  assess_freq = 5,
  sigma_impl = 0.25
  )
ex[["Larger implementation error"]] <- expression(Larger~implementation~error~sigma[h]==0.25)

bg.plot <- function(colour = "#00000019") rect(par("usr")[1],
  par("usr")[3], par("usr")[2], par("usr")[4], col = colour, border =
  FALSE)

source("plot_sim_ts_simple.R")
pdf_eps("plot-various-options-ts-3pops", width = 4, height = 6.4, type = TYPE)

par(mfrow = c(length(eg)-1, 1), mar = c(0,2,0,0), oma = c(4, 2.5, 1, 1), cex = 0.7, las = 1, xpd = FALSE)
for(i in (1:length(eg))[-7]) {
  plot_sim_ts_simple(eg[[i]], text = ex[[i]], pal = col_pal)
  #if(i == 1) bg.plot()
}
axis(1, col = "grey50", tck = -0.05, padj = -1)
mtext("Generation", side = 1, outer = TRUE, line = 1.8, cex = 0.75)
mtext("Metapopulation return abundance", side = 2, outer = TRUE, line = 1.2, cex = 0.75, las = 0)

dev.off()

