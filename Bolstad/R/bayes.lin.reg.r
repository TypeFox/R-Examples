#' Bayesian inference for simple linear regression
#' 
#' This function is used to find the posterior distribution of the simple
#' linear regression slope variable \eqn{\beta}{beta} when we have a random
#' sample of ordered pairs \eqn{(x_{i}, y_{i})} from the simple linear
#' regression model: \deqn{ }{y_i = alpha_xbar + beta*x_i+epsilon_i}\deqn{
#' y_{i} = \alpha_{\bar{x}} + \beta x_{i}+\epsilon_{i} }{y_i = alpha_xbar +
#' beta*x_i+epsilon_i}\deqn{ }{y_i = alpha_xbar + beta*x_i+epsilon_i} where the
#' observation errors are, \eqn{\epsilon_i}{epsilon_i}, independent
#' \eqn{normal(0,\sigma^{2})}{normal(0,sigma^2)} with known variance.
#' 
#' 
#' @param y the vector of responses.
#' @param x the value of the explantory variable associated with each response.
#' @param slope.prior use a ``flat'' prior or a ``normal'' prior. for
#' \eqn{\beta}{beta}
#' @param intcpt.prior use a ``flat'' prior or a ``normal'' prior. for
#' \eqn{\alpha_[\bar{x}]}{alpha_xbar}
#' @param mb0 the prior mean of the simple linear regression slope variable
#' \eqn{\beta}{beta}. This argument is ignored for a flat prior.
#' @param sb0 the prior std. deviation of the simple linear regression slope
#' variable \eqn{\beta}{beta} - must be greater than zero. This argument is
#' ignored for a flat prior.
#' @param ma0 the prior mean of the simple linear regression intercept variable
#' \eqn{\alpha_{\bar{x}}}{alpha_xbar}. This argument is ignored for a flat
#' prior.
#' @param sa0 the prior std. deviation of the simple linear regression variable
#' \eqn{\alpha_{\bar{x}}}{alpha_xbar} - must be greater than zero. This
#' argument is ignored for a flat prior.
#' @param sigma the value of the std. deviation of the residuals. By default,
#' this is assumed to be unknown and the sample value is used instead. This
#' affects the prediction intervals.
#' @param alpha controls the width of the credible interval.
#' @param plot.data if true the data are plotted, and the posterior regression
#' line superimposed on the data.
#' @param pred.x a vector of x values for which the predicted y values are
#' obtained and the std. errors of prediction
#' @return A list will be returned with the following components:
#' \item{post.coef}{the posterior mean of the intecept and the slope}
#' \item{post.coef}{the posterior standard deviation of the intercept the
#' slope} \item{pred.x}{the vector of values for which predictions have been
#' requested. If pred.x is NULL then this is not returned} \item{pred.y}{the
#' vector predicted values corresponding to pred.x. If pred.x is NULL then this
#' is not returned} \item{pred.se}{The standard errors of the predicted values
#' in pred.y. If pred.x is NULL then this is not returned}
#' @keywords misc
#' @examples
#' 
#' ## generate some data from a known model, where the true value of the
#' ## intercept alpha is 2, the true value of the slope beta is 3, and the
#' ## errors come from a normal(0,1) distribution
#' x = rnorm(50)
#' y = 22+3*x+rnorm(50)
#' 
#' ## use the function with a flat prior for the slope beta and a
#' ## flat prior for the intercept, alpha_xbar.
#' 
#' bayes.lin.reg(y,x)
#' 
#' ## use the function with a normal(0,3) prior for the slope beta and a
#' ## normal(30,10) prior for the intercept, alpha_xbar.
#' 
#' bayes.lin.reg(y,x,"n","n",0,3,30,10)
#' 
#' ## use the same data but plot it and the credible interval
#' 
#' bayes.lin.reg(y,x,"n","n",0,3,30,10, plot.data = TRUE)
#' 
#' ## The heart rate vs. O2 uptake example 14.1
#' O2 = c(0.47,0.75,0.83,0.98,1.18,1.29,1.40,1.60,1.75,1.90,2.23)
#' HR = c(94,96,94,95,104,106,108,113,115,121,131)
#' plot(HR,O2,xlab="Heart Rate",ylab="Oxygen uptake (Percent)")
#' 
#' bayes.lin.reg(O2,HR,"n","f",0,1,sigma=0.13)
#' 
#' ## Repeat the example but obtain predictions for HR = 100 and 110
#' 
#' bayes.lin.reg(O2,HR,"n","f",0,1,sigma=0.13,pred.x=c(100,110))
#' 
#' @export bayes.lin.reg

bayes.lin.reg = function(y, x, slope.prior = "flat", 
                        intcpt.prior = "flat", 
                        mb0 = 0, sb0 = 0, ma0 = 0, sa0 = 0, 
                        sigma = NULL, alpha = 0.05, plot.data = FALSE, 
                        pred.x = NULL) {

  if(sum(is.na(y)) > 0 || sum(is.na(x)) > 0)
    stop("Error: x and y may not contain missing values")

  if(length(y) != length(x))
    stop("Error: x and y are unequal lengths")

  if(!is.null(sigma) && sigma <= 0){
    stop("Error: the std. deviation of the resisuals, sigma, must be greater than or equal to zero")
  }

  patn = "n((orm)*al)*"
  patf = "f(lat)*"
  pat = paste0("(", patn, "|", patf, ")")
  
  if(!grepl(pat, slope.prior))
    stop("The slope prior must be one of \"normal\" or \"flat\"")
  
  if(grepl(patn, slope.prior))
    slope.prior = "normal"
  else if(grepl(patf, slope.prior))
    slope.prior = "flat"

  if(!grepl(pat, intcpt.prior))
    stop("The intercept prior must be one of \"normal\" or \"flat\"")
  
  if(grepl(patn, intcpt.prior))
    intcpt.prior = "normal"
  else if(grepl(patf, intcpt.prior))
    intcpt.prior = "flat"

  if(slope.prior == "normal" && sb0 <= 0)
    stop("Error: the prior std. devation sb0 must be greater than zero")

  if(intcpt.prior == "normal" && sa0 <= 0)
    stop("Error: the prior std. devation sa0 must be greater than zero")

  if(alpha <= 0 || alpha > 0.5)
    stop("Error: alpha must be in the range (0, 0.5]")

  if(length(y) <= 2)
    stop("Error: you really should have more than 2 points for a regression!")

  n = length(y)
  x.bar = mean(x)
  y.bar = mean(y)
  x2.bar = mean(x^2)
  xy.bar = mean(x * y)
  y2.bar = mean(y^2)

  b.ls = (xy.bar - x.bar * y.bar) / (x2.bar - x.bar^2)
  fitted = y.bar + b.ls * (x - x.bar)
  residuals = y - fitted

  A0 = y.bar - b.ls * x.bar
  Ax.bar = y.bar

  sigma.known = TRUE
  if(is.null(sigma)){
    sigma.known = FALSE
    sigma = sqrt(sum((y - (Ax.bar + b.ls * (x - x.bar)))^2)/ (n - 2))
    cat(paste("Standard deviation of residuals: ", signif(sigma, 3), "\n"))
  } else {
    cat(paste("Known standard deviation: ", signif(sigma, 3), "\n"))
  }

  SSx = n * (x2.bar - x.bar^2)
  lb = 0
  ub = 0
  prior.b = rep(0, 1001)
  beta = prior.b
  likelihood.b = prior.b
  posterior.b = prior.b
  
  d = as.data.frame(cbind(y, x))
  
  if (slope.prior == "flat") {
    prior.prec.b = 0
    mb0 = 0
    sb0 = 1 # ???
    bnd.mult.b = 4
  } else { 
    prior.prec.b = 1 / sb0^2
    bnd.mult.b = 3
  }
  
  if (intcpt.prior == "flat") {
    prior.prec.a = 0
    ma0 = 0
    sa0 = 1 # ???
    bnd.mult.a = 4
  } else {
    prior.prec.a = 1 / sa0^2
    bnd.mult.a = 3
  }
  
  mdl = bayes.lm(y ~ x, data = d, model = FALSE, prior = list(b0 = c(ma0, mb0), V0 = diag(c(sa0^2, sb0^2))))
  
  ################
  # SLOPE 
  ################
  
  prec.ls = SSx / sigma^2
  sd.ls = sqrt(1 / prec.ls)
  post.prec.b = prior.prec.b + prec.ls

  post.var.b = mdl$post.var[2, 2]
  post.sd.b = sqrt(post.var.b)
  post.mean.b = mdl$post.mean[2]
  
  lb = post.mean.b - bnd.mult.b * post.sd.b
  ub = post.mean.b + bnd.mult.b * post.sd.b
  
  beta = seq(lb, ub, length = 1001)
  
  if (slope.prior == "flat") {
    prior.b = rep(1, 1001)
    norm.const = 0.5 * (2 * sum(prior.b) - prior.b[1] - prior.b[1001] * ((ub - lb) * 0.001))
    prior.b = prior.b / norm.const
  } else {
    prior.b = dnorm(beta, mb0, sb0)
  }
  
  likelihood.b = dnorm(beta, b.ls, sd.ls)
  posterior.b = dnorm(beta, post.mean.b, post.sd.b)
  
  old.par = par(mfrow = c(2, 2))

  y.max = max(c(prior.b, likelihood.b, posterior.b))
  plot(beta, prior.b, type = "l", col = "black", lty = 1, 
       ylim = c(0, 1.1 * y.max), xlab = expression(beta), 
       ylab = "", 
       main = expression(paste("Prior, likelihood and posterior for ", beta, 
           sep = "")), 
       sub = "(slope)")
  lines(beta, likelihood.b, lty = 2, col = "red")
  lines(beta, posterior.b, lty = 3, col = "blue")
  legend("topleft", bty = "n", cex = 0.7, 
         lty = 1:3, col = c("black", "red", "blue"), 
         legend = c("Prior", "Likelihood", "Posterior"))

  ####################################################################################
  
  alpha.xbar = rep(0, 1001)
  prior.a = alpha.xbar
  likelihood.a = alpha.xbar
  posterior.a = alpha.xbar
  
  prec.ls = n / (sigma^2)
  sd.ls = sqrt(1 / prec.ls)
  post.prec.a = prior.prec.a + prec.ls
  
  post.var.a = mdl$post.var[1, 1]
  post.sd.a = sqrt(post.var.a)
  post.mean.a = mdl$post.mean[1]
  
  lb = post.mean.a - bnd.mult.a * post.sd.a
  ub = post.mean.a + bnd.mult.a * post.sd.a

  alpha.xbar = seq(lb, ub, length = 1001)
  
  if(intcpt.prior == "flat") {
    prior.a = rep(1, 1001)
    norm.const = (2 * sum(prior.a) - prior.a[1] - prior.a[1001] * ((ub - lb) / 1000)) / 2
    prior.a = prior.a / norm.const
  } else {
    prior.a = dnorm(alpha.xbar, ma0, sa0)
  }
  
  likelihood.a = dnorm(alpha.xbar, y.bar, sd.ls)
  posterior.a = dnorm(alpha.xbar, post.mean.a, post.sd.a)

  cat(sprintf("%-11s %-14s %-24s\n", " ", "Posterior Mean", "Posterior Std. Deviation"))
  cat(sprintf("%-11s %-14s %-24s\n", " ", "--------------", "------------------------"))
  cat(sprintf("Intercept:  %-14.6g %-24.6g\n", signif(post.mean.a, 4), signif(post.sd.a, 5)))
  cat(sprintf("Slope:      %-14.6g %-24.6g\n", signif(post.mean.b, 4), signif(post.sd.b, 5)))

  y.max = max(c(prior.a, likelihood.a, posterior.a))
  plot(alpha.xbar, prior.a, type = "l", col = "black", lty = 1, 
       ylim = c(0, 1.1 * y.max), xlab = expression(alpha), 
       ylab = "", 
       main = expression(paste("Prior, likelihood and posterior for ", alpha[bar(x)], 
           sep = "")), 
       sub = "(intercept)")
  lines(alpha.xbar, likelihood.a, lty = 2, col = "red")
  lines(alpha.xbar, posterior.a, lty = 3, col = "blue")
  legend("topleft", cex = 0.7, lty = 1:3, col = c("black", "red", "blue"), 
         legend = c("Prior", "Likelihood", "Posterior"), bty = "n")

  if(sigma.known){
    s.e = sqrt(x2.bar - x.bar^2)
    x.lwr = x.bar - 3 * s.e
    x.upr = x.bar + 3 * s.e
    x.values = seq(x.lwr, x.upr, length = 1001)
    pred.y = post.mean.b * (x.values - x.bar) + post.mean.a

    se.pred = sqrt(post.var.a + (x.values - x.bar)^2 * post.var.b + sigma^2)
    t.crit = qt(1 - alpha * .5, n - 2)
    pred.lb = pred.y - t.crit * se.pred
    pred.ub = pred.y + t.crit * se.pred
  } else{
    s.e = sqrt(x2.bar - x.bar^2)
    x.lwr = x.bar - 3 * s.e
    x.upr = x.bar + 3 * s.e
    x.values = seq(x.lwr, x.upr, length = 1001)
    pred.y = post.mean.b * (x.values - x.bar) + post.mean.a

    se.pred = sqrt(post.var.a + (x.values - x.bar)^2 * post.var.b + sigma^2)
    z.crit = qnorm(1 - alpha * 0.5)
    pred.lb = pred.y - z.crit * se.pred
    pred.ub = pred.y + z.crit * se.pred
  }

  y.max = max(pred.ub)
  y.min = min(pred.lb)


  if(plot.data){
    plot(y~x, main = paste("Predicitions with ", round(100 * (1 - alpha))
               ,"% bounds", sep = ""), xlab = "x", ylab = "y", ylim = 1.1 * c(y.min, y.max))
    lines(x.values, pred.y, lty = 1, col = "black")
  } else{
    plot(x.values, pred.y, type = "l", lty = 1, col = "black", 
         main = paste("Predicitions with ", round(100 * (1 - alpha)), 
           "% bounds", sep = ""), xlab = "x", ylab = "y", 
         ylim = 1.1 * c(y.min, y.max))
  }

  lines(x.values, pred.lb, lty = 2, col = "red")
  lines(x.values, pred.ub, lty = 3, col = "blue")

  legend("topleft", lty = 1:3, col = c("black", "red", "blue"), 
         legend = c("Predicted value", 
           paste(round(100 * (1 - alpha)), "% lower bound", sep = ""), 
           paste(round(100 * (1 - alpha)), "% upper bound", sep = "")), 
         cex = 0.7, bty = "n")

  pred.y = NULL
  pred.se = NULL
  if(!is.null(pred.x)){
    pred.y = post.mean.a + post.mean.b * (pred.x - x.bar)
    pred.se = sqrt(post.var.a + (pred.x - x.bar)^2 * post.var.b + sigma^2)
    predicted.values = cbind(pred.x, pred.y, pred.se)
    fmt = "%-8.4g  %-12.4g %-11.5g\n"
    fmtS = "%-6s  %-12s %-11s\n"
    
    cat(sprintf(fmtS, "x", "Predicted y", "SE"))
    cat(sprintf(fmtS, "------", "-----------", "-----------"))
    n.pred.x = length(pred.x)
    for(i in 1:n.pred.x){
      cat(sprintf(fmt, signif(predicted.values[i, 1], 4), 
                       signif(predicted.values[i, 2], 4), 
                       signif(predicted.values[i, 3], 5)))
    }
  }

  par(old.par)
  interceptResults = list(name = 'alpha[0]',
                          param.x = alpha.xbar,
                          prior = prior.a, likelihood = likelihood.a, posterior = posterior.a,
                          mean = post.mean.a,
                          var = post.var.a)
  
  slopeResults = list(name = 'beta',
                      param.x = beta,
                      prior = prior.b, likelihood = likelihood.b, posterior = posterior.b,
                      mean = post.mean.b,
                      var = post.var.b)
  
  class(interceptResults) = "Bolstad"
  class(slopeResults) = "Bolstad"
  
  if(!is.null(pred.x)){
      invisible(list(intercept = interceptResults, 
                     slope = slopeResults,
                     post.coef = c(post.mean.a, post.mean.b), 
                     post.coef.sd = c(post.sd.a, post.sd.b), 
                  pred.x = pred.x, pred.y = pred.y, pred.se = pred.se))
  } else{
      invisible(list(intercept = interceptResults, 
                     slope = slopeResults,
                     post.coef = c(post.mean.a, post.mean.b), 
                     post.coef.sd = c(post.sd.a, post.sd.b)))
  }
}

