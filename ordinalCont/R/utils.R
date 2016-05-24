set.beta_start <- function(x,v){
  vv=ifelse(v<median(v),0,1)
  #as.numeric(-coef(glm(vv~0+x,family=gaussian(link="logit"))))
  as.numeric(-coef(suppressWarnings(glm(vv~0+x,family=binomial(link="logit")))))
  #print(str(x))
  #fit=glm(vv~x,family=binomial(link="logit"))
  #as.numeric(-tail(coef(fit),-1))
  #rep(0, ncol(x))
}


set.glf_start <- function(x,v){
  c(0,1,1)
}



inv.logit <- function(x){
  ifelse(is.finite(x),exp(x)/(1+exp(x)),sign(x)*Inf)
}


##Functions used in plot.ocm to bootstrapping data (random-x or fixed-x resampling) and find CIs.
rnd.x.bootstrap <- function(data, indices, fit){
  data <- data[indices,]
  mod <- update(fit, .~., data = data)
  coefficients(mod)
}

fix.x.bootstrap <- function(data, indices, fit){
  W = as.numeric(-fit$x %*% fit$coefficients[1:fit$len_beta])
  data$new_v <- g_glf_inv(W + residuals(fit)[indices], tail(fit$coefficients,3))
  mod <- update(fit, new_v ~., data = data)
  coefficients(mod)
}

param.bootstrap <- function(data, indices, fit){
  W = as.numeric(-fit$x %*% fit$coefficients[1:fit$len_beta])
  data$new_v <- g_glf_inv(W + rlogis(fit$sample.size), tail(fit$coefficients,3))
  mod <- update(fit, new_v ~., data = data)
  coefficients(mod)
}



## returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)



## rnd generation - multivariate normal
mvrnormR <- function(n, mu, sigma) {
    ncols <- ncol(sigma)
    mu <- rep(mu, each = n) ## not obliged to use a matrix (recycling)
    mu + matrix(rnorm(n * ncols), ncol = ncols) %*% chol(sigma)
}


## my version of (a)ghQaud. It is the same as in the fastGHQuad package, but deals nicely with arrays (i.e. more than one record per cluster).
##FIXME not used
my.ghQuad <- function (f, rule, ...) 
{
  apply(rule$w * f(rule$x, ...), 2, sum)
}


my.aghQuad <- function (g, muHat, sigmaHat, rule, ...) 
{
  z <- muHat + sqrt(2) * sigmaHat * rule$x
  wStar <- exp(rule$x * rule$x + log(rule$w))
  val <- sqrt(2) * sigmaHat * apply(wStar * g(z, ...), 2, sum)
  return(val)
}


format.perc <- function(probs, digits)
    ## Not yet exported, maybe useful in other contexts:
    ## quantile.default() sometimes uses a version of it
    paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), "%")