### R code from vignette source 'rstan_vignette.Rnw'

###################################################
### code chunk number 1: rstan_vignette.Rnw:77-79
###################################################
options(width=80)
library(rstan)


###################################################
### code chunk number 2: lookup
###################################################
lookup("dnorm")
tail(lookup("~")) # looks up all Stan sampling statements
lookup(dwilcox)   # no corresponding Stan function


###################################################
### code chunk number 3: expose_stan_functions
###################################################
model_code <-
'
functions {
  real standard_normal_rng() {
    return normal_rng(0,1);
  }
}
model {}
'
expose_stan_functions(stanc(model_code = model_code))
standard_normal_rng(seed = 1)


###################################################
### code chunk number 4: rstan_vignette.Rnw:333-337
###################################################
schools_data <- 
  list(J=8, 
  y=c(28,  8, -3,  7, -1,  1, 18, 12),
  sigma=c(15, 10, 16, 11,  9, 11, 10, 18))


###################################################
### code chunk number 5: callstan
###################################################
J <- 8
y <- c(28,  8, -3,  7, -1,  1, 18, 12)
sigma <- c(15, 10, 16, 11,  9, 11, 10, 18)
library(rstan)
fit1 <- stan(file="schools.stan",
             # better to add explicitly include: data=schools_data, 
             iter=2000, chains=4, cores=1)


###################################################
### code chunk number 6: rstan_vignette.Rnw:392-394
###################################################
print(fit1, pars=c("theta", "mu", "tau", "lp__"), 
      probs=c(.1,.5,.9))


###################################################
### code chunk number 7: rstan_vignette.Rnw:463-464
###################################################
y <- as.array(y)


###################################################
### code chunk number 8: stanfit_plot
###################################################
plot(fit1)


###################################################
### code chunk number 9: stanfit_tplot
###################################################
traceplot(fit1, pars = "tau")


###################################################
### code chunk number 10: rstan_vignette.Rnw:557-564
###################################################
s <- extract(fit1, pars = c("theta", "mu"), permuted = TRUE)
names(s)
dim(s$theta)
dim(s$mu)
s2 <- extract(fit1, pars = "theta", permuted = FALSE)
dim(s2)
dimnames(s2)


###################################################
### code chunk number 11: rstan_vignette.Rnw:614-619
###################################################
# all chains combined
summary(do.call(rbind, args = get_sampler_params(fit1, inc_warmup = TRUE)),
        digits = 2)
# each chain separately
lapply(get_sampler_params(fit1, inc_warmup = TRUE), summary, digits = 2)


###################################################
### code chunk number 12: pairs_plot
###################################################
pairs(fit1, pars = c("eta", "theta"), include = FALSE, las = 1)


###################################################
### code chunk number 13: optimizer
###################################################
ocode <- "
  data {
    int<lower=1> N;
    real y[N];
  } 
  parameters {
    real mu;
  } 
  model {
    y ~ normal(mu, 1);
  } 
"

sm <- stan_model(model_code = ocode)
y2 <- rnorm(20)
mean(y2)
message("This will crash on SPARC")
optimizing(sm, data = list(y = y2, N = length(y2)), hessian = TRUE)


