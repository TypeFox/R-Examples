## =============================================================================
## Test of the Baranyi Growth model, cf. Baranyi, L (1995),
## Int. J. Food Microbiology, doi:10.1016/0168-1605(94)00121-L
##
## Note: the original model formulation works in log space, so that
##       y(t), y_0 and y_max are all given as natural log.
##  This is advantageous for model fitting, but here we use y_0 and K in
##  untransformed space to be compatible with the growth parameters of
##  most other models. The downside is, that we need box constraints in most
##  cases and possibly more iterations.
##
## Author: Thomas Petzoldt, TU Dresden
## License: GPL >= 2, https://www.gnu.org/licenses/
## Please cite our work when using this package.
## =============================================================================



library("growthrates")
library("lattice")

data(bactgrowth)
splitted.data <- multisplit(value ~ time|strain + conc + replicate, data = bactgrowth)

## use a single data set
dat <- splitted.data[[23]]

## initial parameters and bocx constraints
p   <- c(y0=0.03, mumax=.1, K=0.1, h0=1)

lower   <- c(y0=0.001, mumax=1e-2, K=0.005, h0=0)
upper   <- c(y0=0.1,   mumax=1,    K=0.5,   h0=10)

## fit model
fit <- fit_growthmodel(FUN=grow_baranyi, p=p, time=dat$time, y=dat$value,
                       lower = lower, upper = upper,
                       control=list(trace=TRUE)
)

## coefficients and plot
coef(fit)
plot(fit)


## fit growth models to all data using (log transformed residuals)
system.time(
  L <- all_growthmodels(grow_baranyi, p=p, data=bactgrowth,
                        grouping = c("strain", "conc", "replicate"),
                        lower = lower, upper=upper,
                        log="y")
)

## same with formula interface
system.time(
  L <- all_growthmodels(value ~ grow_baranyi(time, parms) | strain + conc + replicate,
                        data = bactgrowth,
                        p=p, lower = lower, upper=upper,
                        log = "y"
                        )
)

par(mfrow=c(4,3))
par(mar=c(2.5,4,2,1))
plot(L, log="y")

par(mfrow=c(4,3))
plot(L)

res <- results(L)
xyplot(mumax ~ log(conc + 1)| strain, data=res, layout=c(3,1))
