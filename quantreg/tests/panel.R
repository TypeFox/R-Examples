# Example rqss vs rq.fit.panel vs rq.fit.lasso  
# Slightly modified version of a test problem of Stefan Bache and co.
require(quantreg)
source("rq.fit.panel.R")
set.seed(1917)

tt <- 3   # time periods
nn <- 3  # individuals

#  Generate some data:
some.data <- data.frame(
    id=rep(1:nn, each=tt)             #ids
   ,ic=1                              #intercept
   ,x1=rnorm(nn*tt)                   #regressors
   ,x2=runif(nn*tt)
   ,alpha=rep(runif(nn)*2 ,each=tt)   #fixed effects
)

# response:
some.data$y <- some.data$x1 + some.data$x2 + some.data$alpha + rnorm(nn*tt)

lambda <- .2
tau <- 0.25


# Fit with rq.fit.panel
fit1 <- rq.fit.panel(cbind(1, some.data$x1, some.data$x2), some.data$y
            ,rep(1:nn, each=tt)
            ,tau=tau
            ,w=1
            ,lambda=lambda)

# fit with rqss (using the global debug variable from the panel function
# so remember to run that too! :-) :
fit2 <- rqss(y ~ ic + x1 + x2 + as.factor(id) - 1 
            ,data=some.data
            ,method="lasso"
            ,tau=tau
            ,lambda=c(0, 0, 0, rep(lambda, nn)))

# fit with rq.lasso
fit3 <- rq(y ~ ic + x1 + x2 + as.factor(id) - 1 
            ,data=some.data
            ,tau=tau
            ,method="lasso"
            ,lambda=c(0, 0, 0, rep(lambda, nn)))

# Print coefficients for comaparison.
comparefit <- function()
{
  compmat <- cbind(fit1$coef,fit2$coef, fit3$coef)
  colnames(compmat) <- c("rq.fit.panel", "rqss", "rq.fit.lasso")
  noquote(formatC(compmat, format="f", digits=8, width=12))
}

cat("all(rq.fit.panel$coef==rq$coef): ", all.equal(fit1$coef,fit3$coef), "\n")
cat("all(rq.fit.$coef==rqss$coef): ", all.equal(fit2$coef,fit3$coef), "\n")
cat("all(rqss$coef==rq.fit.panel$coef): ", all.equal(fit2$coef,fit1$coef), "\n")

