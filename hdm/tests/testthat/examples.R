#' @example
#' ## DGP
#' set.seed(2)
#' n <- 250
#' p <- 100
#' px <- 10
#' X <- matrix(rnorm(n*p), ncol=p)
#' beta <- c(rep(2,px), rep(0,p-px))
#' intercept <- 1
#' P <- exp(intercept + X %*% beta)/(1+exp(intercept + X %*% beta))
#' y <- numeric(length=250)
#' for(i in 1:n){
#'   y[i] <- sample(x=c(1,0), size=1, prob=c(P[i],1-P[i]))
#' }
#' ## fit rlogisticlasso object
#'  rlogisticlasso.reg <- rlogisticlasso(x=X, y=y)
#' 
#'  ## methods
#' summary(rlogisticlasso.reg, all=F)
#' print(rlogisticlasso.reg)
#' predict(rlogisticlasso.reg, type="response")
#' X3 <- matrix(rnorm(n*p), ncol=p)
#' predict(rlogisticlasso.reg, newdata=X3)



#' @examples
#' ## DGP
#' n <- 250
#' p <- 100
#' px <- 10
#' X <- matrix(rnorm(n*p), ncol=p)
#' beta <- c(rep(2,px), rep(0,p-px))
#' intercept <- 1
#' y <- intercept + X %*% beta + rnorm(n)
#' ## fit rlassoLM object with inference on three variables
#' rlassoLM.reg <- rlassoLM(x=X, y=y, index=c(1,7,20))
#' ## methods
#' summary(rlassoLM.reg)
#' print(rlassoLM.reg)
#' confint(rlassoLM.reg, level=0.9)


set.seed(2)
n <- 250
p <- 100
px <- 10
X <- matrix(rnorm(n*p), ncol=p)
beta <- c(rep(2,px), rep(0,p-px))
intercept <- 1
P <- exp(intercept + X %*% beta)/(1+exp(intercept + X %*% beta))
y <- numeric(length=250)
for(i in 1:n){
  y[i] <- sample(x=c(1,0), size=1, prob=c(P[i],1-P[i]))
}
## fit rlogisticlasso object
rlogisticlasso.reg <- rlogisticlasso(x=X, y=y)

## methods
summary(rlogisticlasso.reg, all=F)
print(rlogisticlasso.reg)
head(predict(rlogisticlasso.reg, type="response"))
X3 <- matrix(rnorm(n*p), ncol=p)
head(predict(rlogisticlasso.reg, newdata=X3))



library(hdm)
## DGP
set.seed(2)
n <- 250
p <- 100
px <- 10
X <- matrix(rnorm(n*p), ncol=p)
beta <- c(rep(2,px), rep(0,p-px))
intercept <- 1
P <- exp(intercept + X %*% beta)/(1+exp(intercept + X %*% beta))
y <- numeric(length=250)
for(i in 1:n){
   y[i] <- sample(x=c(1,0), size=1, prob=c(P[i],1-P[i]))
}
 ## fit rlogisticlasso object
  rlogisticlasso.reg <- rlogisticlasso(x=X, y=y)
  ## methods
summary(rlogisticlasso.reg, all=F)
print(rlogisticlasso.reg)
predict(rlogisticlasso.reg, type="response")
X3 <- matrix(rnorm(n*p), ncol=p)
predict(rlogisticlasso.reg, newdata=X3)
