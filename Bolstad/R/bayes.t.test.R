#' Bayesian t-test
#' 
#' @description Performs one and two sample t-tests (in the Bayesian hypothesis testing framework) on vectors of data
#' @param x a (non-empty) numeric vector of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @param alternative a character string specifying the alternative hypothesis, must be one of 
#' \code{"two.sided"} (default), \code{"greater"} or \code{"less"}. You can specify just the initial 
#' letter.
#' @param mu a number indicating the true value of the mean (or difference in means if you are performing a two sample test).
#' @param paired a logical indicating whether you want a paired t-test.
#' @param var.equal a logical variable indicating whether to treat the two variances as being equal. 
#' If \code{TRUE} (default) then the pooled variance is used to estimate the variance otherwise the 
#' Welch (or Satterthwaite) approximation to the degrees of freedom is used. The unequal variance case is
#' implented using Gibbs sampling.
#' @param conf.level confidence level of interval.
#' @param prior a character string indicating which prior should be used for the means, must be one of
#' \code{"jeffreys"} (default) for independent Jeffreys' priors on the unknown mean(s) and variance(s), 
#' or \code{"joint.conj"} for a joint conjugate prior.
#' @param m if the joint conjugate prior is used then the user must specify a prior mean in the one-sample
#' or paired case, or two prior means in the two-sample case. Note that if the hypothesis is that there is no difference
#' between the means in the two-sample case, then the values of the prior means should usually be equal, and if so, 
#' then their actual values are irrelvant.This parameter is not used if the user chooses a Jeffreys' prior.
#' @param n0 if the joint conjugate prior is used then the user must specify the prior precision 
#' or precisions in the two sample case that represent our level of uncertainty
#' about the true mean(s). This parameter is not used if the user chooses a Jeffreys' prior.
#' @param sig.med if the joint conjugate prior is used then the user must specify the prior median
#' for the unknown standard deviation. This parameter is not used if the user chooses a Jeffreys' prior.
#' @param kappa if the joint conjugate prior is used then the user must specify the degrees of freedom
#' for the inverse chi-squared distribution used for the unknown standard deviation. Usually the default
#' of 1 will be sufficient. This parameter is not used if the user chooses a Jeffreys' prior.
#' @param sigmaPrior If a two-sample t-test with unequal variances is desired then the user must choose between
#' using an chi-squared prior ("chisq") or a gamma prior ("gamma") for the unknown population standard deviations.
#' This parameter is only used if \code{var.equal} is set to \code{FALSE}.
#' @param nIter Gibbs sampling is used when a two-sample t-test with unequal variances is desired.
#' This parameter controls the sample size from the posterior distribution.
#' @param nBurn Gibbs sampling is used when a two-sample t-test with unequal variances is desired.
#' This parameter controls the number of iterations used to burn in the chains before the procedure 
#' starts sampling in order to reduce correlation with the starting values. 
#' @param formula a formula of the form \code{lhs ~ rhs} where lhs is a numeric variable giving the data values and rhs a factor with two 
#' levels giving the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see \code{\link{model.frame}}) containing 
#' the variables in the formula formula. By default the variables are taken from \code{environment(formula)}.
#' @param subset currently ingored.
#' @param na.action currently ignored.
#' @param ... any additional arguments
#' @return A list with class "htest" containing the following components:
#'  \item{statistic}{the value of the t-statistic.}
#'  \item{parameter}{the degrees of freedom for the t-statistic.}                                                                                        
#'   \item{p.value}{the p-value for the test.}"                                                                                                              
#'   \item{conf.int}{a confidence interval for the mean appropriate to the specified alternative hypothesis.}
#'   \item{estimate}{the estimated mean or difference in means depending on whether it was a one-sample test or a two-sample test.}
#'   \item{null.value}{the specified hypothesized value of the mean or mean difference depending on whether it was a one-sample test or a two-sample test.}
#'   \item{alternative}{a character string describing the alternative hypothesis.}
#'   \item{method}{a character string indicating what type of t-test was performed.}
#'   \item{data.name}{a character string giving the name(s) of the data.}
#'   \item{result}{an object of class \code{Bolstad}}
#' @examples
#' bayes.t.test(1:10, y = c(7:20))      # P = .3.691e-01
#' 
#' ## Same example but with using the joint conjugate prior
#' ## We set the prior means equal (and it doesn't matter what the value is)
#' ## the prior precision is 0.01, which is a prior standard deviation of 10
#' ## we're saying the true difference of the means is between [-25.7, 25.7]
#' ## with probability equal to 0.99. The median value for the prior on sigma is 2
#' ## and we're using a scaled inverse chi-squared prior with 1 degree of freedom
#' bayes.t.test(1:10, y = c(7:20), var.equal = TRUE, prior = "joint.conj", 
#'              m = c(0,0), n0 =  rep(0.01, 2), sig.med = 2)
#' 
#' ##' Same example but with a large outlier. Note the assumption of equal variances isn't sensible
#' bayes.t.test(1:10, y = c(7:20, 200)) # P = .1979    -- NOT significant anymore
#' 
#' ## Classical example: Student's sleep data
#' plot(extra ~ group, data = sleep)
#' 
#' ## Traditional interface
#' with(sleep, bayes.t.test(extra[group == 1], extra[group == 2]))
#' 
#' ## Formula interface
#' bayes.t.test(extra ~ group, data = sleep)
#' @author R Core with Bayesian internals added by James Curran
#' @export 
bayes.t.test = function(x,  ...){
  UseMethod("bayes.t.test")
}

#' @describeIn bayes.t.test Bayesian t-test
#' @export
bayes.t.test.default = function(x, y = NULL, alternative = c("two.sided", "less", "greater"),
       mu = 0, paired = FALSE, var.equal = TRUE,
       conf.level = 0.95, prior = c("jeffreys", "joint.conj"), 
       m = NULL, n0 = NULL, sig.med = NULL, kappa = 1, 
       sigmaPrior = "chisq", nIter = 10000, nBurn = 1000, ...){
  
  prior = match.arg(prior)
  
  if(prior == "joint.conj" & (is.null(m) | is.null(n0) | is.null(sig.med) | kappa < 1)){
    m1 = "If you are using the joint conjugate prior, you need so specify:"
    m2 = "the prior mean(s), the prior precision(s), the prior median standard deviation," 
    m3 = "and the degrees of freedom associated with the prior for the standard deviation"
    stop(paste(m1, m2, m3, sep = "\n"))
  }
    
  ## Shamelessly copied from t.test.default
  
  alternative = match.arg(alternative)
  
  if (!missing(mu) && (length(mu) != 1 || is.na(mu))) 
    stop("'mu' must be a single number")
  
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
                                 conf.level < 0 || conf.level > 1)) 
    stop("'conf.level' must be a single number between 0 and 1")
  
  param.x = NULL
  
  tstat = 0
  df = 0
  pval = 0 
  cint = 0
  estimate = 0
  method = NULL
    
  if (!is.null(y)) {
    dname = paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    if (paired) 
      xok = yok = complete.cases(x, y)
    else {
      yok = !is.na(y)
      xok = !is.na(x)
    }
    y = y[yok]
  } else {
    dname = deparse(substitute(x))
    if (paired) 
      stop("'y' is missing for paired test")
    xok = !is.na(x)
    yok = NULL
  }
  
  x = x[xok]
  
  if (paired) {
    x = x - y
    y = NULL
  }
  
  nx = length(x)
  mx = mean(x)
  SSx = sum((x-mx)^2)
  vx = var(x)
 
  bolstadResult = NULL
   
  if (is.null(y)) { ## one sample or paired
    if (nx < 2) 
      stop("not enough 'x' observations")
    
    stderr = sqrt(vx/nx)
    if (stderr < 10 * .Machine$double.eps * abs(mx)) 
      stop("data are essentially constant")
    
    name = 'mu'
    name = if(!paired) 'mu' else 'mu[d]'
    
    if(prior == "jeffreys"){
      S1 = SSx
      kappa1 = nx - 1
      npost = nx
      
      mpost = mx
      se.post = sqrt(S1 / kappa1 / nx)
      df = kappa1
      
      param.x = seq(mx - 4 * sqrt(vx), mx + 4 * sqrt(vx), length = 200)
      prior = 1 / diff(range(x))
      likelihood = dnorm(mx, param.x, se.post)
      std.x = (param.x - mpost) / se.post
      posterior = dt(std.x, df = df)
      
      bolstadResult = list(name = name, param.x = param.x, 
                           prior = prior, likelihood = likelihood, posterior = posterior,
                           mean = mpost,
                           var = se.post^2,
                           cdf = function(x)pt((x - mpost) / se.post, df = df),
                           quantileFun = function(probs, ...){se.post * qt(probs, df = df, ...) + mpost})
      class(bolstadResult) = 'Bolstad'
      
      tstat = (mpost - mu) / se.post
      estimate = mpost
    }else{
      S0 = qchisq(0.5, kappa) * sig.med^2
      S1 = SSx + S0
      kappa1 = nx + kappa
      npost = n0 + nx
      sigma.sq.B = (S1 + (n0 * nx / kappa1) * (mx - m)^2)/npost
      
      mpost = (nx * mx + n0 * m) / npost
      se.post = sqrt(sigma.sq.B / kappa1)
      df = kappa1
      
      estimate = mpost
      tstat = (mpost - mu) / se.post
      
      lb = min(mpost - 4 * se.post, m - 4 * sqrt(1 / n0))
      ub = max(mpost + 4 * se.post, m + 4 * sqrt(1 / n0))
      param.x = seq(lb, ub, length = 200)
      prior = dnorm(param.x, m, sqrt(1 / n0))
      likelihood = dnorm(mx, param.x, se.post)
      std.x = (param.x - mpost) / se.post
      posterior = dt(std.x, df = df)
    }
    
    method = if (paired) 
      "Paired t-test"
    else 
      "One Sample t-test"
    
  } else { ## two sample
    ny = length(y)
    if (nx < 1 || (!var.equal && nx < 2)) 
      stop("not enough 'x' observations")
    
    if (ny < 1 || (!var.equal && ny < 2)) 
      stop("not enough 'y' observations")
    
    if (var.equal && nx + ny < 3) 
      stop("not enough observations")
    
    name = 'mu[d]'
    
    my = mean(y)
    vy = var(y)
    stderr = sqrt((sum((x - mx)^2) + sum((y - my)^2) / (nx + ny - 2)) * (1 / nx + 1 / ny))
    
    if (stderr < 10 * .Machine$double.eps * max(abs(mx), abs(my))) 
      stop("data are essentially constant")
    
    SSp = sum((x - mx)^2) + sum((y - my)^2)
    
    method = paste(if (!var.equal) 
      "Gibbs", "Two Sample t-test")
    
    estimate = c(mx, my) ## this may get changed elsewhere
    names(estimate) = c("posterior mean of x", "posterior mean of y")
    
    lb = mx - my - 4 * sqrt(vx/nx + vy/ny)
    ub = mx - my + 4 * sqrt(vx/nx + vy/ny)
    param.x = seq(lb, ub, length = 1000)
    name = 'mu[1]-mu[2]'

    if (var.equal) {
      if(prior == "jeffreys"){
        kappa1 = nx + ny -2
        kappa1 = nx + ny - 2
        sigma.sq.Pooled = SSp / kappa1 
        
        mpost = mx - my
        se.post = sqrt(sigma.sq.Pooled * (1/nx + 1/ny))
        df = kappa1
        
        prior = 1 / diff(range(param.x))
        likelihood = dnorm(mx - my, param.x, se.post)
        posterior = dt((param.x - mpost)/se.post, df)
        
        bolstadResult = list(name = name, param.x = param.x, 
                             prior = prior, likelihood = likelihood, posterior = posterior,
                             mean = mpost,
                             var = se.post^2,
                             cdf = function(x)pt((x - mpost) / se.post, df = df),
                             quantileFun = function(probs, ...){se.post * qt(probs, df = df, ...) + mpost})
        class(bolstadResult) = 'Bolstad'
        
        estimate = c(mx, my)
        names(estimate) = c("posterior mean of x", "posterior mean of y")
        
        tstat = (mpost - mu) / se.post
        
      }else{
        kappa1 = kappa + nx + ny
        n1post = nx + n0[1]
        n2post = ny + n0[2]
        S = qchisq(0.5, 1) * sig.med^2
        S1 = S + SSp
        m1post = (nx * mx + n0[1] * m[1]) / n1post
        m2post = (ny * my + n0[2] * m[2]) / n2post
        sigma.sq.B = S1 / kappa1
        mpost = m1post - m2post
        se.post = sqrt(sigma.sq.B * (1/n1post + 1/n2post))
        df = kappa1
        
        prior = dnorm(param.x, m[1] - m[2], sqrt(sum(1/n0)))
        likelihood = dnorm(mx - my, param.x, se.post)
        posterior = dt((param.x - mpost)/se.post, df)
        
        bolstadResult = list(name = name, param.x = param.x, 
                             prior = prior, likelihood = likelihood, posterior = posterior,
                             mean = mpost,
                             var = se.post^2,
                             cdf = function(x)pt((x - mpost) / se.post, df = df),
                             quantileFun = function(probs, ...){se.post * qt(probs, df = df, ...) + mpost})
        class(bolstadResult) = 'Bolstad'
        
        estimate = c(m1post, m2post)
        names(estimate) = c("posterior mean of x", "posterior mean of y")
        
        tstat = (mpost - mu)/se.post
      }
    }else {
      res = bayes.t.gibbs(x, y, nIter, nBurn, sigmaPrior)
      se.post = sd(res$mu.diff)
      
      d = density(res$mu.diff, from = param.x[1], to = param.x[length(param.x)])
      param.x = d$x
      likelihood = dnorm(mx - my, param.x, se.post)
      posterior = d$y
      mpost = mean(res$mu.diff)
      vpost = var(res$mu.diff)
      
      
      bolstadResult = list(name = name, param.x = d$x, 
                           prior = NULL, likelihood = likelihood, posterior = posterior,
                           mean = mpost,
                           var = vpost,
                           cdf = function(x){r  = sintegral(param.x, posterior);
                                             Fx = splinefun(r$x, r$y);
                                             return(Fx(x))},
                           quantileFun = function(probs, ...){quantile(res$mu.diff, probs = probs, ...)})
      class(bolstadResult) = 'Bolstad'
      
      estimate = c(mean(res$mu.x), mean(res$mu.y))
      names(estimate) = c("posterior mean of x", "posterior mean of y")
      
      tstat = mean(res$tstat)
      se.post = sd(res$mu.diff)
      snx = mean(res$sigma.sq.x / nx)
      sny = mean(res$sigma.sq.y / ny)
      df = (snx + sny)^2 / (snx^2 / (nx - 1) + sny^2 / (ny - 1))
    }
  }
  
  if (alternative == "less") {
    pval = pt(tstat, df)
    cint = c(-Inf, tstat + qt(conf.level, df))
  }
  else if (alternative == "greater") {
    pval = pt(tstat, df, lower.tail = FALSE)
    cint = c(tstat - qt(conf.level, df), Inf)
  }
  else {
    pval = 2 * pt(-abs(tstat), df)
    alpha = 1 - conf.level
    cint = qt(1 - alpha/2, df)
    cint = tstat + c(-cint, cint)
  }
  cint = mu + cint * se.post
  names(tstat) = "t"
  names(df) = "df"
  names(mu) = if (paired || !is.null(y)) 
    "difference in means"
  else "mean"
  attr(cint, "conf.level") = conf.level
  rval = list(statistic = tstat, parameter = df, p.value = pval, 
               conf.int = cint, estimate = estimate, null.value = mu, 
               alternative = alternative, method = method, data.name = dname,
              result = bolstadResult)
  class(rval) = "htest"
  return(rval)
}

#' @describeIn bayes.t.test Bayesian t-test
#' @export
bayes.t.test.formula = function(formula, data, subset, na.action, ...){
  ## shamelessly hacked from t.test.formula
  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]), 
                                                                  "term.labels")) != 1L)) 
    stop("'formula' missing or incorrect")
  m = match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame()))) 
    m$data = as.data.frame(data)
  m[[1L]] = quote(stats::model.frame)
  m$... = NULL
  mf = eval(m, parent.frame())
  DNAME = paste(names(mf), collapse = " by ")
  names(mf) = NULL
  response = attr(attr(mf, "terms"), "response")
  g = factor(mf[[-response]])
  if (nlevels(g) != 2L) 
    stop("grouping factor must have exactly 2 levels")
  DATA = setNames(split(mf[[response]], g), c("x", "y"))
  y = do.call("bayes.t.test", c(DATA, list(...)))
  y$data.name = DNAME
  if (length(y$estimate) == 2L) 
    names(y$estimate) = paste("mean in group", levels(g))
  y
}