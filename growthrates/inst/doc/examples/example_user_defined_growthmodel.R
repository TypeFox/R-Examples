## =============================================================================
## The following example shows how to create user-defined growth models
##
## Author: Thomas Petzoldt, TU Dresden
## License: GPL >= 2, https://www.gnu.org/licenses/
## Please cite our work when using this package.
## =============================================================================


library("growthrates")


## =============================================================================
## create a "growthmodel" with interfaces compatible to package growthrates
## see ?growthmodel for details
## Note:
##   It is essential to use consistent names for parameters and initial values!
## =============================================================================

grow_logistic_yshift <- function(time, parms) {
  with(as.list(parms), {
    y <- (K * y0) / (y0 + (K - y0) * exp(-mumax * time)) + y_shift
    return(as.matrix(data.frame(time=time, y=y, log_y=log(y))))
  })
}

grow_logistic_yshift <- growthmodel(grow_logistic_yshift,
                                    c("y0", "mumax", "K", "y_shift"))


## =============================================================================
## Fit the model
## =============================================================================

x <- seq(5, 100, 5)

y <- c(2.1, 2.3, 5, 4.7, 4.3, 6.9, 8.2, 11.5, 8.8, 10.2, 14.5, 12.5,
       13.6, 12.7, 14.2, 12.5, 13.8, 15.1, 12.7, 14.9) + 5

fit <- fit_growthmodel(grow_logistic_yshift,
                       p=c(y0=1, mumax=0.1, K=10, K = 10, y_shift=1),
                       time=x, y=y)
plot(fit)
summary(fit)
