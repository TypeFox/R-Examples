## ------------------------------------------------------------------------
knitr::opts_chunk$set(comment=NA, warning=FALSE, message=FALSE)

## ----eval=FALSE----------------------------------------------------------
#  library("devtools")
#  install_github("cboettig/pmc2")

## ------------------------------------------------------------------------
library("pmc")
library("geiger")
library("ouch")
library("ggplot2")
library("tidyr")
library("dplyr")
library("gridExtra")

## ------------------------------------------------------------------------
cores <- 1 # parallel cores available
nboot <- 50 # too small, but vignette gotta build faster
phy <- sim.bdtree(n=10)
dat <- sim.char(rescale(phy, "lambda", .5), 1)[,1,]
out <- pmc(phy, dat, "BM", "lambda", nboot = nboot, mc.cores = cores)

## ----fig.width=7---------------------------------------------------------
dists <- data.frame(null = out$null, test = out$test)
dists %>% 
  gather(dist, value) %>%
  ggplot(aes(value, fill = dist)) + 
  geom_density(alpha = 0.5) + 
  geom_vline(xintercept = out$lr)


## ------------------------------------------------------------------------
data(geospiza)
bm_v_lambda <- pmc(geospiza$phy, geospiza$dat[, "wingL"],
  "BM", "lambda", nboot = nboot, mc.cores = cores) 


## ------------------------------------------------------------------------
lambdas <- bm_v_lambda$par_dists %>% filter(comparison=="BB", parameter=="lambda") 

## ------------------------------------------------------------------------
est <- coef(bm_v_lambda[["B"]])[["lambda"]]

## ------------------------------------------------------------------------
ggplot(lambdas) + geom_histogram(aes(value)) +
      geom_vline(xintercept=est)

## ------------------------------------------------------------------------
bm_v_lambda$par_dists %>% filter(comparison=="BB", parameter=="sigsq") %>% 
ggplot() + geom_histogram(aes(sqrt(value))) 

## ----eval=FALSE----------------------------------------------------------
#  reshape2::cast(lambdas, comparison ~ parameter, function(x)
#       quantile(x, c(.05, .95)), value=c("lower", "upper"))

## ------------------------------------------------------------------------
data(anoles)
tree <- with(anoles, ouchtree(node, ancestor, time / max(time), species))

ou3v4 <- pmc(tree, log(anoles["size"]), modelA = "hansen", modelB = "hansen", 
             optionsA = list(regimes = anoles["OU.LP"]), 
             optionsB = list(regimes = anoles["OU.4"]),
             nboot = nboot, sqrt.alpha = 1, sigma = 1, mc.cores = cores)

ou3v15 <- pmc(tree, log(anoles["size"]), "hansen", "hansen", 
             list(regimes = anoles["OU.LP"]), 
             list(regimes = anoles["OU.15"]),
             nboot = nboot, sqrt.alpha = 1, sigma = 1, mc.cores = cores)
                   
ou1v3 <- pmc(tree, log(anoles["size"]), "hansen", "hansen", 
             list(regimes = anoles["OU.1"]), 
             list(regimes = anoles["OU.LP"]),
             nboot = nboot, sqrt.alpha = 1, sigma = 1, mc.cores = cores)
 
ou0v1 <- pmc(tree, log(anoles["size"]), "brown", "hansen", 
             list(), 
             list(regimes = anoles["OU.1"], sqrt.alpha = 1, sigma = 1),
             nboot = nboot, mc.cores = cores)

## ----fig.width=7, fig.height=7-------------------------------------------
results <- bind_rows(
  data.frame(comparison = "ou3v4", null = ou3v4$null, test = ou3v4$test, lr = ou3v4$lr),
  data.frame(comparison = "ou3v15", null = ou3v15$null, test = ou3v15$test, lr = ou3v15$lr),
  data.frame(comparison = "ou1v3", null = ou1v3$null, test = ou1v3$test, lr = ou1v3$lr),
  data.frame(comparison = "ou0v1", null = ou0v1$null, test = ou0v1$test, lr = ou0v1$lr)) %>%
gather(variable, value, - comparison, - lr) 
ggplot(results) + 
  geom_density(aes(value, fill = variable), alpha=0.7) + 
  geom_vline(aes(xintercept=lr)) +
  facet_wrap(~ comparison, scales="free")

## ----fig.width=7, fig.height=7, eval=FALSE-------------------------------
#  a <- plot(ou3v4, A = "OU.3", B ="OU.4") + labs(title = "OU.3 vs OU.4")
#  b <- plot(ou3v15, A = "OU.3", B ="OU.15") + labs(title = "OU.3 vs OU.15")
#  c <- plot(ou1v3, A = "OU.1", B = "OU.3")  + labs(title = "OU.1 vs OU.3")
#  d <- plot(ou0v1, A = "BM", B = "OU.1")  + labs(title = "BM vs OU.1")
#  grid.arrange(a, b, c, d, ncol = 2)

