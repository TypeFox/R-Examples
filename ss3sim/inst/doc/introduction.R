## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ---- install-cran, eval=FALSE-------------------------------------------
#  install.packages("ss3sim")

## ---- install-and-load, eval=FALSE---------------------------------------
#  # install.packages("devtools") # if needed
#  devtools::install_github("ss3sim/ss3sim", dependencies = TRUE)
#  
#  # If you would like the vignettes available with the GitHub development version:
#  devtools::install_github("ss3sim/ss3sim", dependencies = TRUE, build_vignettes = TRUE)
#  
#  # If you would like to run simulations in parallel, also install:
#  install.packages(c("doParallel", "foreach"))

## ---- load-package, eval=FALSE-------------------------------------------
#  library("ss3sim")

## ---- help, eval=FALSE---------------------------------------------------
#  ?ss3sim
#  help(package = "ss3sim")
#  browseVignettes("ss3sim")

## ---- citation, cache=FALSE----------------------------------------------
citation("ss3sim")

## ---- locate-folders-----------------------------------------------------
library(ss3sim)
d <- system.file("extdata", package = "ss3sim")
case_folder <- paste0(d, "/eg-cases")
om <- paste0(d, "/models/cod-om")
em <- paste0(d, "/models/cod-em")

## ---- case-file-checks, eval=FALSE---------------------------------------
#  run_ss3sim(iterations = 1, scenarios =
#    c("D0-E0-F0-M0-cod",
#      "D1-E0-F0-M0-cod",
#      "D0-E1-F0-M0-cod",
#      "D1-E1-F0-M0-cod"),
#    case_files = list(F = "F", D = c("index", "lcomp", "agecomp"),
#      E = "E", M = "M"),
#    case_folder = case_folder, om_dir = om,
#    em_dir = em)

## ------------------------------------------------------------------------
recdevs_det <- matrix(0, nrow = 100, ncol = 20)

## ---- deterministic-runs, eval=FALSE-------------------------------------
#  run_ss3sim(iterations = 1:20,
#    scenarios = c("D100-E100-F0-M0-cod", "D100-E101-F0-M0-cod"),
#    case_files = list(F = "F", D = c("index", "lcomp", "agecomp"),
#      E = "E", M = "M"),
#    case_folder = case_folder, om_dir = om, em_dir = em,
#    bias_adjust = TRUE, user_recdevs = recdevs_det)

## ---- deterministic-runs-expand, eval=FALSE------------------------------
#  x <- expand_scenarios(list(D = 100, E = 100:101, F = 0, M = 0),
#    species = "cod")
#  run_ss3sim(iterations = 1:20, scenarios = x,
#    case_folder = case_folder, om_dir = om, em_dir = em,
#    bias_adjust = TRUE, user_recdevs = recdevs_det,
#    case_files = list(F = "F", D = c("index", "lcomp", "agecomp"),
#      E = "E", M = "M"))

## ---- stochastic-runs, eval=FALSE----------------------------------------
#  run_ss3sim(iterations = 1:100, scenarios =
#    c("D0-E0-F0-M0-cod",
#      "D1-E0-F0-M0-cod",
#      "D0-E1-F0-M0-cod",
#      "D1-E1-F0-M0-cod"),
#    case_files = list(F = "F", D = c("index", "lcomp", "agecomp"),
#      E = "E", M = "M"),
#    case_folder = case_folder, om_dir = om,
#    em_dir = em, bias_adjust = TRUE)

## ---- get-results, eval=FALSE--------------------------------------------
#  get_results_all(user_scenarios =
#    c("D100-E100-F0-M0-cod",
#      "D100-E101-F0-M0-cod",
#      "D0-E0-F0-M0-cod",
#      "D1-E0-F0-M0-cod",
#      "D0-E1-F0-M0-cod",
#      "D1-E1-F0-M0-cod"))

## ---- read-output, eval=FALSE--------------------------------------------
#  scalar_dat <- read.csv("ss3sim_scalar.csv")
#  ts_dat <- read.csv("ss3sim_ts.csv")

## ---- load-output--------------------------------------------------------
data("ts_dat", package = "ss3sim")
data("scalar_dat", package = "ss3sim")

## ---- transform-output---------------------------------------------------
scalar_dat <- transform(scalar_dat,
  steep = (SR_BH_steep_om - SR_BH_steep_em)/SR_BH_steep_om,
  logR0 = (SR_LN_R0_om - SR_LN_R0_em)/SR_LN_R0_om,
  depletion = (depletion_om - depletion_em)/depletion_om,
  SSB_MSY = (SSB_MSY_em - SSB_MSY_om)/SSB_MSY_om,
  SR_sigmaR = (SR_sigmaR_em - SR_sigmaR_om)/SR_sigmaR_om,
  NatM = (NatM_p_1_Fem_GP_1_em - NatM_p_1_Fem_GP_1_om)/
     NatM_p_1_Fem_GP_1_om)

ts_dat <- transform(ts_dat,
  SpawnBio = (SpawnBio_em - SpawnBio_om)/SpawnBio_om,
  Recruit_0 = (Recruit_0_em - Recruit_0_om)/Recruit_0_om)
ts_dat <- merge(ts_dat, scalar_dat[,c("scenario", "replicate",
    "max_grad")])

scalar_dat_det <- subset(scalar_dat, E %in% c("E100", "E101"))
scalar_dat_sto <- subset(scalar_dat, E %in% c("E0", "E1"))
ts_dat_det <- subset(ts_dat, E %in% c("E100", "E101"))
ts_dat_sto <- subset(ts_dat, E %in% c("E0", "E1"))

## ---- reshape-scalars----------------------------------------------------
scalar_dat_long <- reshape2::melt(scalar_dat[,c("scenario", "D", "E",
  "replicate", "max_grad", "steep", "logR0", "depletion", "SSB_MSY",
  "SR_sigmaR", "NatM")], id.vars = c("scenario", "D", "E",
  "replicate", "max_grad"))
scalar_dat_long <- plyr::rename(scalar_dat_long,
  c("value" = "relative_error"))

## ---- relative-error-boxplots-det, fig.height=7, fig.width=5, fig.cap="Relative error box plots for deterministic runs. In case E100, *M* is fixed at the true value; in E101 we estimate *M*. In case D100, the standard deviation on the survey index observation error is 0.001."----
library("ggplot2")
p <- ggplot(subset(scalar_dat_long, E %in% c("E100", "E101") &
       variable != "SR_sigmaR"), aes(D, relative_error)) +
     geom_boxplot() +
     geom_hline(aes(yintercept = 0), lty = 2) +
     facet_grid(variable~E) +
     theme_bw() + ylim(-0.4, 0.4)
print(p)

## ---- plot-sto-ts, fig.height=5, fig.width=7, fig.cap="Time series of relative error in spawning stock biomass."----
p <- ggplot(ts_dat_sto, aes(x = year)) + xlab("Year") +
    theme_bw() + geom_line(aes(y = SpawnBio, group = replicate,
    colour = max_grad), alpha = 0.3, size = 0.15) + facet_grid(D~E) +
    scale_color_gradient(low = "gray", high = "red")
print(p)

## ---- ssb-ts-plots, fig.height=5, fig.width=7, cache=TRUE, fig.cap="Spawning stock biomass time series."----
p <- ggplot(ts_dat_sto, aes(year, SpawnBio_em, group = replicate)) +
  geom_line(alpha = 0.3, aes(colour = max_grad)) + facet_grid(D~E) +
  scale_color_gradient(low = "darkgrey", high = "red") + theme_bw()
print(p)

## ---- relative-error-boxplots-sto, fig.height=7, fig.width=5, cache=TRUE, fig.cap="Relative error box plots for stochastic runs. In case E0, *M* is fixed at the true value; in E1 we estimate *M*. In case D1, the standard deviation on the survey index observation error is 0.4. In case D0, the standard deviation is quartered representing an increase in survey sampling effort."----
p <- ggplot(subset(scalar_dat_long, E %in% c("E0", "E1")),
       aes(D, relative_error)) +
     geom_boxplot() + geom_hline(aes(yintercept = 0), lty = 2) +
     facet_grid(variable~E) +
     geom_jitter(aes(colour = max_grad),
       position = position_jitter(height = 0, width = 0.1),
       alpha = 0.4, size = 1.5) +
     scale_color_gradient(low = "darkgrey", high = "red") +
     theme_bw()
print(p)

## ---- ss3sim-base-eg, eval=FALSE-----------------------------------------
#  d <- system.file("extdata", package = "ss3sim")
#  om <- paste0(d, "/models/cod-om")
#  em <- paste0(d, "/models/cod-em")
#  
#  F0 <- list(years = 1913:2012, years_alter = 1913:2012, fvals =
#    c(rep(0, 25), rep(0.114, 75)))
#  
#  index1 <- list(fleets = 2, years = list(seq(1974, 2012, by = 2)),
#    sds_obs = list(0.1))
#  
#  lcomp1 <- list(fleets = c(1, 2), Nsamp = list(100, 100), years =
#    list(1938:2012, seq(1974, 2012, by = 2)), lengthbin_vector = NULL,
#    cpar = c(1, 1))
#  
#  agecomp1 <- list(fleets = c(1, 2), Nsamp = list(100, 100), years =
#    list(1938:2012, seq(1974, 2012, by = 2)), agebin_vector = NULL,
#    cpar = c(1, 1))
#  
#  E0 <- list(natM_type = "1Parm", natM_n_breakpoints = NULL,
#    natM_lorenzen = NULL, natM_val = c(NA,-1), par_name =
#    "LnQ_base_3_CPUE", par_int = NA, par_phase = -1, forecast_num = 0)
#  
#  M0 <- list(NatM_p_1_Fem_GP_1 = rep(0, 100))
#  
#  ss3sim_base(iterations = 1:20, scenarios = "D1-E0-F0-M0-cod",
#    f_params = F0, index_params = index1, lcomp_params = lcomp1,
#    agecomp_params = agecomp1, estim_params = E0, tv_params = M0,
#    om_dir = om, em_dir = em)

## ---- alternative-case-lists, eval=FALSE---------------------------------
#  case_folder <- system.file("extdata", "eg-cases", package = "ss3sim")
#  om <- system.file("extdata", "models/cod-om", package = "ss3sim")
#  em <- system.file("extdata", "models/cod-em", package = "ss3sim")
#  files <- list.files(case_folder, pattern = "0-cod")
#  
#  temp_case_folder <- "ss3sim-example-cases"
#  dir.create(temp_case_folder)
#  file.copy(paste0(case_folder, "/", files), temp_case_folder)
#  
#  # now make X, Y, and Z case files:
#  setwd(temp_case_folder)
#  file.copy("index0-cod.txt", "X0-cod.txt")
#  file.copy("lcomp0-cod.txt", "Y0-cod.txt")
#  file.copy("agecomp0-cod.txt", "Z0-cod.txt")
#  setwd("..")
#  
#  # our custom specification:
#  case_files <- list(F = "F", X = "index", Y = "lcomp", Z = "agecomp")
#  
#  # and use run_ss3sim() with our new case_file list:
#  run_ss3sim(iterations = 1,
#    scenarios = "X0-Y0-Z0-F0-cod",
#    case_folder = temp_case_folder,
#    om_dir = om, em_dir = em,
#    case_files = case_files)

## ---- custom-case-eg, eval=FALSE-----------------------------------------
#  case_files = list(D = c("index", "lcomp", "agecomp"), F = "F", S = "S")

## ---- parallel-one, eval=FALSE-------------------------------------------
#  require(doParallel)
#  registerDoParallel(cores = 4)

## ---- parallel-two, eval=FALSE-------------------------------------------
#  require(foreach)
#  getDoParWorkers()
#  
#  #> [1] 4

## ---- parallel-three, eval=FALSE-----------------------------------------
#  run_ss3sim(iterations = 1, scenarios =
#    c("D1-E0-F0-cod",
#      "D2-E0-F0-cod",
#      "D1-E1-F0-cod",
#      "D2-E1-F0-cod"),
#    case_files =
#      list(F = "F", D = c("index", "lcomp", "agecomp"), E = "E"),
#    case_folder = case_folder, om_dir = om, em_dir = em,
#    bias_adjust = TRUE, parallel = TRUE)

## ---- parallel-iterations1, eval=FALSE-----------------------------------
#  run_ss3sim(iterations = 1:2, scenarios = "D0-F0-cod",
#    case_folder = case_folder, om_dir = om, em_dir = em,
#    parallel = TRUE, parallel_iterations = TRUE)

## ---- parallel-iterations2, eval=FALSE-----------------------------------
#  run_ss3sim(iterations = 1:2, scenarios = "D0-F0-cod",
#    case_folder = case_folder, om_dir = om, em_dir = em,
#    parallel = TRUE, parallel_iterations = TRUE, bias_nsim = 2,
#    bias_adjust = TRUE)

