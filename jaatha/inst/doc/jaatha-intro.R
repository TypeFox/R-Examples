## ----eval=FALSE----------------------------------------------------------
#  vignette("jaatha-evolution", package = jaatha)

## ----load_jaatha---------------------------------------------------------
library(jaatha)
set.seed(112233)

## ----data_obs------------------------------------------------------------
data_obs <- c(2, 8, 0, 6, 1, 3, 2, 2, 0, 7)
data_obs

## ----sim_func------------------------------------------------------------
sim_func <- function(x) rpois(10, x)
sim_func(c(p1 = 1, p2 = 10))

## ----sum_stats-----------------------------------------------------------
sum_stats <- list(create_jaatha_stat("id", function(x, opts) x))

## ----par_ranges----------------------------------------------------------
par_ranges <- matrix(c(0.1, 0.1, 10, 10), 2, 2)
rownames(par_ranges) <- c("x", "y")
colnames(par_ranges) <- c("min", "max")
par_ranges

## ----create_model--------------------------------------------------------
jaatha_model <- create_jaatha_model(sim_func, par_ranges, sum_stats)

## ----create_data---------------------------------------------------------
jaatha_data <- create_jaatha_data(data_obs, jaatha_model)

## ----execute_jaatha------------------------------------------------------
estimates <- jaatha(jaatha_model, jaatha_data, 
                    sim = 50, repetitions = 2, verbose = FALSE)

## ----print_estimates-----------------------------------------------------
estimates

