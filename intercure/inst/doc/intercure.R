## ---- echo = FALSE, message = FALSE, warning = FALSE---------------------
require(intercure)

## ---- eval = FALSE-------------------------------------------------------
#  sim_bch(10)

## ---- echo = FALSE-------------------------------------------------------
set.seed(2)
knitr::kable(sim_bch(10), align = 'c')

## ---- message = FALSE, warning = FALSE, results = 'hide'-----------------
# hypothetical interval censored dataset with a cure fraction
set.seed(2)
cureset <- sim_bch(100)

# allocates the estimated parameters and covariance matrix
output <- inter_bch(cureset, 
                        cureset$L, cureset$R, 
                        c("xi1", "xi2"))

## ------------------------------------------------------------------------
output$par

## ------------------------------------------------------------------------
indiv <- c(1, 1, 0)
est_effect <- output$par
cf <- exp(-exp(est_effect[1:3]%*%indiv))
cf

## ---- eval = FALSE-------------------------------------------------------
#  sim_frailty(10)

## ---- echo = FALSE-------------------------------------------------------
set.seed(3)
knitr::kable(sim_frailty(10), align = 'c')

## ---- message = FALSE, warning = FALSE, results = 'hide'-----------------
# hypothetical interval censored dataset with a cure fraction
set.seed(2)
cureset <- sim_frailty(100)

# allocates the estimated parameters and covariance matrix
output <- inter_frailty(cureset, 
                        cureset$L, cureset$R, cureset$delta, 
                        c("xi1", "xi2"), c("xi1", "xi2"),
                        M = 10, max_n = 30, burn_in = 10)



## ------------------------------------------------------------------------
output$par

## ------------------------------------------------------------------------
indiv <- c(1, 1, 0)
est_effect <- output$par
cf <- exp(-exp(est_effect[1:3]%*%indiv)/2)
cf

## ---- eval = FALSE-------------------------------------------------------
#  require(parallel)
#  require(doParallel)
#  cl <- makeCluster(4)
#  registerDoParallel(cl)
#  output <- inter_frailty(cureset,
#                          cureset$L, cureset$R, cureset$delta,
#                          c("xi1", "xi2"), c("xi1", "xi2"),
#                          M = 10, par_cl = cl)
#  
#  
#  
#  stopCluster(cl)

## ---- eval = FALSE-------------------------------------------------------
#  sim_frailty_cl(15, nclus = 3)

## ---- echo = FALSE-------------------------------------------------------
set.seed(2)
knitr::kable(sim_frailty_cl(15, nclus = 3), align = 'c')

## ---- message = FALSE, warning = FALSE, results = 'hide'-----------------
# hypothetical interval censored dataset with a cure fraction
set.seed(2)
cureset <- sim_frailty_cl(120, nclus = 3)

# allocates the estimated parameters and covariance matrix
output <- inter_frailty_cl(cureset, 
                        cureset$L, cureset$R, cureset$delta, 
                        c("xi1", "xi2"), c("xi1", "xi2"),
                        grp = cureset$clus, M = 30, max_n = 30,
                        burn_in = 10)

## ------------------------------------------------------------------------
output$par

