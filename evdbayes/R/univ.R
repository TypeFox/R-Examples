"prior.quant" <-
function(prob = 10^-(1:3), shape, scale, trendsd = 0)
{
    if(length(prob) != 3 || mode(prob) != "numeric")
        stop("`prob' must be a numeric vector of length three")
    if(min(prob) <= 0 || max(prob) >= 1)
        stop("`prob' must contain values in (0,1)")
    if(any(diff(prob) >= 0))
        stop("`prob' must be a vector of decreasing values")
    if(length(shape) != 3 || mode(shape) != "numeric")
        stop("`shape' must be a numeric vector of length three")
    if(length(scale) != 3 || mode(scale) != "numeric")
        stop("`scale' must be a numeric vector of length three")
    if(any(c(shape, scale) <= 0))
        stop("`shape'  and `scale' must contain positive values")
    if(length(trendsd) != 1 || mode(trendsd) != "numeric")
        stop("`trendsd' must be a numeric vector of length one")
    if(trendsd < 0)
        stop("`trendsd' must be positive")
    structure(list(prior = "dprior.quant", prob = prob, shape = shape,
                   scale = scale, trendsd = trendsd), class = "evprior") 
}

"prior.prob" <-
function(quant, alpha, trendsd = 0)
{
    if(length(quant) != 3 || mode(quant) != "numeric")
        stop("`quant' must be a numeric vector of length three")
    if(any(diff(quant) <= 0))
        stop("`quant' must be a vector of increasing values")
    if(length(alpha) != 4 || mode(alpha) != "numeric")
        stop("`alpha' must be a numeric vector of length four")
    if(any(alpha <= 0))
        stop("`alpha' must contain positive values")
    if(length(trendsd) != 1 || mode(trendsd) != "numeric")
        stop("`trendsd' must be a numeric vector of length one")
    if(trendsd < 0)
        stop("`trendsd' must be positive")
    structure(list(prior = "dprior.prob", quant = quant, alpha = alpha,
              trendsd = trendsd), class = "evprior")
}

"prior.norm" <-
function(mean, cov, trendsd = 0)
{
    if(length(mean) != 3 || mode(mean) != "numeric")
        stop("`mean' must be a numeric vector of length three")
    if(!is.matrix(cov) || any(dim(cov) != 3) || mode(cov) != "numeric")
        stop("`cov' must be a symmetric three by three matrix")
    if(length(trendsd) != 1 || mode(trendsd) != "numeric")
        stop("`trendsd' must be a numeric vector of length one")
    if(trendsd < 0)
        stop("`trendsd' must be positive")
    if(any(abs(cov - t(cov)) > .Machine$double.eps^0.5))
        warning("`cov' may not be symmetric")
    eg <- eigen(cov, symmetric = TRUE, only.values = TRUE)$values
    if(any(eg <= 0))
        warning("`cov' may not be positive definite")
 
    icov <- solve(cov)
    icov <- icov[row(icov) >= col(icov)] 
    structure(list(prior = "dprior.norm", mean = mean, icov = icov,
              trendsd = trendsd), class = "evprior")
}

"prior.loglognorm" <-
  function(mean, cov, trendsd = 0)
{
  if(length(mean) != 3 || mode(mean) != "numeric")
    stop("`mean' must be a numeric vector of length three")
  if(!is.matrix(cov) || any(dim(cov) != 3) || mode(cov) != "numeric")
    stop("`cov' must be a symmetric three by three matrix")
  if(length(trendsd) != 1 || mode(trendsd) != "numeric")
    stop("`trendsd' must be a numeric vector of length one")
  if(trendsd < 0)
    stop("`trendsd' must be positive")
  if(any(abs(cov - t(cov)) > .Machine$double.eps^0.5))
    warning("`cov' may not be symmetric")
  eg <- eigen(cov, symmetric = TRUE, only.values = TRUE)$values
  if(any(eg <= 0))
    warning("`cov' may not be positive definite")
  
  icov <- solve(cov)
  icov <- icov[row(icov) >= col(icov)] 
  structure(list(prior = "dprior.loglognorm", mean = mean, icov = icov,
                 trendsd = trendsd), class = "evprior")
}

# PRIOR DENSITIES

"dprior.quant" <-
## Computes log prior density for (mu,sigma,xi) based on gamma
## distributions for quantile differences corresponding to three
## specified probabilities. Includes optional normal trend for
## location.
function(par, prob = 10^-(1:3), shape, scale, trendsd)
{
    .C("dprior_quant",
        par, prob, shape, scale, trendsd, dns = double(1),
        PACKAGE = "evdbayes")$dns
}

"dprior.prob" <-
## Computes log prior density for (mu,sigma,xi) based on beta
## distributions for probability ratios corresponding to three
## specified quantiles. Includes optional normal trend for location.
function(par, quant, alpha, trendsd)
{
    .C("dprior_prob",
        par, quant, alpha, trendsd, dns = double(1),
        PACKAGE = "evdbayes")$dns
}

"dprior.norm" <-
## Computes log prior density for (mu,sigma,xi) based on a trivariate
## normal distribution for (mu,log(sigma),xi). Includes optional normal
## trend for location.
function(par, mean, icov, trendsd)
{
    .C("dprior_norm",
        par, mean, icov, trendsd, dns = double(1),
        PACKAGE = "evdbayes")$dns
}

"dprior.loglognorm" <-
## Computes log prior density for (loc,scale,shape) based on a trivariate
## normal distribution for (log(loc),log(scale),shape). Includes optional
## normal trend for location.
function(par, mean, icov, trendsd)
{
    .C("dprior_loglognorm",
       par, mean, icov, trendsd, dns = double(1),
       PACKAGE = "evdbayes")$dns
}

# LIKELIHOODS

"gevlik" <-
## Computes log-likelihood of gev model at (mu,sigma,xi). Gumbel is
## computed for small shape.
function(par, data, trend)
{
    nas <- !is.na(data)
    data <- data[nas]
    if(!missing(trend)) trend <- trend[nas]

    gevlik2(par = par, data = data, trend = trend)
}

"gevlik2" <-
## Computes log-likelihood of gev model at (mu,sigma,xi). Gumbel is
## computed for small shape.
function(par, data, trend)
{
    n <- length(data)
    if(missing(trend))
      lik <- .C("gevlik",
          data, n, par, dns = double(1),
          PACKAGE = "evdbayes")$dns
    else
      lik <- .C("gevlikt",
          data, n, par, trend, dns = double(1),
          PACKAGE = "evdbayes")$dns
    lik
}

"gpdlik" <-
## Computes log-likelihood of gpd model at
## (loc,scale,shape). Exponential is computed for small shape.
function(par, data, trend)
{
    nas <- !is.na(data)
    data <- data[nas]
    if(!missing(trend)) trend <- trend[nas]

    gpdlik2(par = par, data = data, trend = trend)
}

"gpdlik2" <-
## Computes log-likelihood of gpd model at
## (loc,scale,shape). Exponential is computed for small shape.
function(par, data, trend)
{
    n <- length(data)
    if(missing(trend))
      lik <- .C("gpdlik",
          data, n, par, dns = double(1),
          PACKAGE = "evdbayes")$dns
    else
      lik <- .C("gpdlikt",
          data, n, par, trend, dns = double(1),
          PACKAGE = "evdbayes")$dns
    lik
}

"pplik" <-
## Computes log-likelihood of Poission process model at (mu,sigma,xi)
## with threshold u and npy observations per block. Gumbel is computed
## for small shape.
function(par, data, thresh, noy, trend, exact = FALSE)
{
    n <- length(data)
    thresh <- rep(thresh, length.out = n)
    nan <- !is.na(data)
    data <- data[nan]
    thresh <- thresh[nan] 
    if(!missing(trend)) trend <- trend[nan]
    
    hd <- (data > thresh)
    data <- data[hd]
    if(length(data) == 0) stop("no data above threshold")
    if(!missing(trend)) htrend <- trend[hd]

    if(!exact) {
        thn <- seq(1, length(thresh), length = length(data))
        thresh <- thresh[thn]
        if(!missing(trend)) trend <- trend[thn]
    }
      
    pplik2(par = par, data = data, thresh = thresh,
           noy = noy, trend = trend, htrend = htrend)
}

"pplik2" <-
## Computes log-likelihood of Poission process model at (mu,sigma,xi)
## with threshold u and npy observations per block. Gumbel is computed
## for small shape.
function(par, data, thresh, noy, trend, htrend)
{   
    nh <- length(data)
    n <- length(thresh)
    if(missing(trend))
      lik <- .C("pplik",
          data, nh, par, thresh, n, noy, dns = double(1),
          PACKAGE = "evdbayes")$dns
    else
      lik <- .C("pplikt",
          data, nh, par, thresh, n, noy, trend, htrend,
                dns = double(1), PACKAGE = "evdbayes")$dns
    lik
}

"oslik" <-
## Computes log-likelihood of gev order statistics model at
## (mu,sigma,xi). Gumbel is computed for small shape.
function(par, data, trend)
{
    nas <- !apply(is.na(data), 1, all)
    data <- data[nas, ,drop = FALSE]
    if(!missing(trend)) trend <- trend[nas]
    thresh <- apply(data, 1, min, na.rm = TRUE)
    rvec <- as.integer(cumsum(rowSums(!is.na(data))))
    data <- t(data)
    data <- data[!is.na(data)]
    
    oslik2(par = par, data = data, trend = trend, thresh = thresh,
           rvec = rvec)
}

"oslik2" <-
## Computes log-likelihood of gev order statistics model at
## (mu,sigma,xi). Gumbel is computed for small shape.
function(par, data, trend, thresh, rvec)
{
    m <- length(thresh)
    n <- length(data)
    if(missing(trend))
      lik <- .C("oslik",
          data, thresh, n, m, par, dns = double(1),
          PACKAGE = "evdbayes")$dns
    else
      lik <- .C("oslikt",
          data, thresh, n, m, rvec, par, trend, dns = double(1),
          PACKAGE = "evdbayes")$dns
    lik
}

# POSTERIOR DENSITIES

"dpost" <-
## Computes log-posterior density at (mu,sigma,xi). Includes optional
## normal trend for location.
function(par, prior, lh, mix = FALSE, pMassProb, normPi0, pMass, ...)
## mix : logical. Is there a point mass for the prior distribution
{
  if ( missing(pMass) && mix )
    stop("pMass should be present.")
  if ( missing(pMassProb) && mix )
    stop("pMassProb should be present.")
  if ( missing(normPi0) && mix )
    stop("normPi0 should be present.")
  
  dprior <- do.call(prior$prior, c(list(par = par), prior[-1]))

  if (mix == TRUE){
    if (par[3] == pMass)
      dprior <- dprior + log(pMassProb) - log(normPi0)
    else
      dprior <- dprior + log(1- pMassProb)
  }
  
  switch(lh,
         gev = dprior + gevlik2(par, ...),
         gpd = dprior + gpdlik2(par, ...),
         pp = dprior + pplik2(par, ...),
         os = dprior + oslik2(par, ...),
         none = dprior)
}


# MAXIMIZE POSTERIOR DENSITY

"mposterior" <-
function(init, prior, lh = c("none","gev","gpd","pp","os"),
         method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),
         lower = -Inf, upper = Inf, control = list(),
         hessian = FALSE, ...)
{
    if (!inherits(prior, "evprior")) 
        stop("`prior' must be a prior distribution")
    if(!prior$trendsd && (length(init) != 3 || mode(init) != "numeric"))
        stop("`init' must be a numeric vector of length three")
    if(prior$trendsd && (length(init) != 4 || mode(init) != "numeric"))
        stop("`init' must be a numeric vector of length four")
    if(init[2] <= 0)
        stop("initial value of scale parameter must be positive")

    lh <- match.arg(lh)
    ndpost <- function(par, prior, lh, ...) {
      if(par[2] <= 0) return(Inf)
      dprior <- do.call(prior$prior, c(list(par = par), prior[-1]))
      switch(lh,
             gev = -dprior - gevlik(par, ...),
             gpd = -dprior - gpdlik(par, ...),
             pp = -dprior - pplik(par, ...),
             os = -dprior - oslik(par, ...),
             none = -dprior)
    }

    inittest <- ndpost(par = init, prior = prior, lh = lh, ...) 
    if(is.infinite(inittest))
        stop("density is zero at initial parameter values")
    
    optim(init, ndpost, method = method, lower = lower,
          upper = upper, control = control, prior = prior,
          lh = lh, ...)
}

# SAMPLER

"gibbs" <-
function(n, init, prior, lh, ..., psd, thin, burn)
{
    if(prior$trendsd) np <- 4
    else np <- 3
    param <- c("mu","sigma","xi")
    if(prior$trendsd) param <- c(param, "mutrend")
    dpst <- function(a) dpost(a, prior, lh, ...)
    
    mc <- .Call("gibbs", as.integer(n), as.integer(np),
                as.integer(thin), as.double(init), as.double(psd),
                quote(dpst(x)), new.env(), PACKAGE = "evdbayes")
    
    naccept <- mc[[2]]
    nexternal <- mc[[3]]
    mc <- matrix(mc[[1]], ncol = np, byrow = TRUE)
    dimnames(mc) <- list(seq(0, n, thin), param)
    mc <- mc[as.numeric(rownames(mc)) >= burn, , drop = FALSE]
    nexternal[np+1] <- sum(nexternal)/np
    naccept[np+1] <- sum(naccept)/np
    ar <- round(c(naccept/n, nexternal/n),2)
    ardn <- list(c("acc.rates","ext.rates"), c(param,"total"))
    ar <- matrix(ar, ncol = np+1, byrow = TRUE, dimnames = ardn)
    structure(mc, ar = ar)
}

"gibbs.mix" <-
function(n, init, prior, lh, ..., psd, thin, burn, pMassProb,
         normPi0, xitilde, pMass, cv)
{
    if(prior$trendsd) np <- 4
    else np <- 3
    param <- c("mu","sigma","xi")
    if(prior$trendsd) param <- c(param, "mutrend")
    dpst <- function(a) dpost(a, prior, lh, mix = TRUE, pMassProb,
                              normPi0, pMass, ...)
    
    mc <- .Call("gibbsmix", as.integer(n), as.integer(np),
                as.integer(thin), as.double(init), as.double(psd),
                quote(dpst(x)), new.env(), as.double(pMassProb),
                as.double(xitilde), as.double(pMass), as.double(cv),
                PACKAGE = "evdbayes")
    
    naccept <- mc[[2]]
    nexternal <- mc[[3]]
    propType <- mc[[4]]
    names(propType) <- c('Type 1', 'Type 2', 'Normal')
    naccType <- mc[[5]]
    naccType <- round(naccType / propType[-3], 2)
    naccType <- matrix(naccType, ncol =2, byrow = TRUE,
                       dimnames=list(c("acc.rates"),
                         c("Type 1", "Type 2")))
    mc <- matrix(mc[[1]], ncol = np, byrow = TRUE)
    dimnames(mc) <- list(seq(0, n, thin), param)
    mc <- mc[as.numeric(rownames(mc)) >= burn, , drop = FALSE]
    nexternal[np+1] <- sum(nexternal)/np
    naccept[np+1] <- sum(naccept)/np
    ar <- round(c(naccept/propType[3], nexternal/propType[3]),2)
    ardn <- list(c("acc.rates","ext.rates"),
                 c(param,"total"))
    ar <- matrix(ar, ncol = np+1, byrow = TRUE, dimnames = ardn)
    structure(mc, ar = ar, propType = propType, arType = naccType)
}

# WRAPPER TO SAMPLER

"posterior" <-
function(n, init, prior, lh = c("none","gev","gpd","pp","os"), ..., psd,
         burn = 0, thin = 1)
{
    if (!inherits(prior, "evprior")) 
        stop("`prior' must be a prior distribution")
    if(!prior$trendsd && (length(init) != 3 || mode(init) != "numeric"))
        stop("`init' must be a numeric vector of length three")
    if(prior$trendsd && (length(init) != 4 || mode(init) != "numeric"))
        stop("`init' must be a numeric vector of length four")
    if(init[2] <= 0)
        stop("initial value of scale parameter must be positive")
    if(!prior$trendsd && (length(psd) != 3 || mode(psd) != "numeric"))
        stop("`psd' must be a numeric vector of length three")
    if(prior$trendsd && (length(psd) != 4 || mode(psd) != "numeric"))
        stop("`psd' must be a numeric vector of length four")
    if(any(psd <= 0))
        stop("`psd'  must contain positive values")
    if(any(psd <= 0))
        stop("`psd'  must contain positive values")
    if(burn > n)
        stop("burn-in period is larger than run length")
    if(thin > n)
        stop("thinning interval is larger than run length")
    
    lh <- match.arg(lh)
    ar <- list(...)
    
    if(lh == "pp") {
        nn <- length(ar$data)
        ar$thresh <- rep(ar$thresh, length.out = nn)
        nan <- !is.na(ar$data)
        ar$data <- ar$data[nan]
        ar$thresh <- ar$thresh[nan]
        if(!is.null(ar$trend)) ar$trend <- ar$trend[nan]
        
        hd <- (ar$data > ar$thresh)
        ar$data <- ar$data[hd]
        if(length(ar$data) == 0) stop("no data above threshold")
        if(!is.null(ar$trend)) ar$htrend <- ar$trend[hd]

        if(is.null(ar$exact)) ar$exact <- FALSE
        if(!ar$exact) {
          thn <- seq(1, length(ar$thresh), length = length(ar$data))
          ar$thresh <- ar$thresh[thn]
          if(!is.null(ar$trend)) ar$trend <- ar$trend[thn]
        }
        ar$exact <- NULL
    }
    
    if(lh == "gev") {
        nas <- !is.na(ar$data)
        ar$data <- ar$data[nas]  
        if(!is.null(ar$trend)) ar$trend <- ar$trend[nas]
    }
    if(lh == "gpd") {
        nas <- !is.na(ar$data)
        ar$data <- ar$data[nas]  
        if(!is.null(ar$trend)) ar$trend <- ar$trend[nas]
    }
    if(lh == "os") {
        nas <- !apply(is.na(ar$data), 1, all)
        ar$data <- ar$data[nas, ,drop = FALSE]  
        if(!is.null(ar$trend)) ar$trend <- ar$trend[nas]
        ar$thresh <- apply(ar$data, 1, min, na.rm = TRUE)
        ar$rvec <- as.integer(cumsum(rowSums(!is.na(ar$data))))
        ar$data <- t(ar$data)
        ar$data <- ar$data[!is.na(ar$data)]
    }

    initar <- c(list(par = init, prior = prior, lh = lh), ar)
    inittest <- do.call("dpost", initar) 
    if(is.infinite(inittest))
        stop("density is zero at initial parameter values")

    ar <- c(list(n = n, init = init, prior = prior, lh = lh), ar,
            list(psd = psd, thin = thin, burn = burn))
    do.call("gibbs", ar)
}


"posterior.mix" <-
function(n, init, prior, lh = c("none","gev","gpd"), ..., psd,
         pMassProb, normPi0, xitilde, pMass = 0, cv,
         burn = 0, thin = 1)
{
    if (!inherits(prior, "evprior")) 
        stop("`prior' must be a prior distribution")
    if(!prior$trendsd && (length(init) != 3 || mode(init) != "numeric"))
        stop("`init' must be a numeric vector of length three")
    if(prior$trendsd && (length(init) != 4 || mode(init) != "numeric"))
        stop("`init' must be a numeric vector of length four")
    if(init[2] <= 0)
        stop("initial value of scale parameter must be positive")
    if(!prior$trendsd && (length(psd) != 3 || mode(psd) != "numeric"))
        stop("`psd' must be a numeric vector of length three")
    if(prior$trendsd && (length(psd) != 4 || mode(psd) != "numeric"))
        stop("`psd' must be a numeric vector of length four")
    if(any(psd <= 0))
        stop("`psd'  must contain positive values")
    if(any(psd <= 0))
        stop("`psd'  must contain positive values")
    if(burn > n)
        stop("burn-in period is larger than run length")
    if(thin > n)
        stop("thinning interval is larger than run length")
    if(pMassProb == 1 && init[3] != pMass)
      stop("initial values should match with the point Mass")
    if(pMassProb == 0 && init[3] == pMass)
      stop("initial values should not match with the point Mass")
    
    lh <- match.arg(lh)
    ar <- list(...)

    if(lh == "gev") {
      cv <- -log(1 - cv)
      nas <- !is.na(ar$data)
      ar$data <- ar$data[nas]  
      if(!is.null(ar$trend)) ar$trend <- ar$trend[nas]
    }
    if(lh == "gpd") {
      cv <- 1 - cv
      nas <- !is.na(ar$data)
      ar$data <- ar$data[nas]  
      if(!is.null(ar$trend)) ar$trend <- ar$trend[nas]
    }

    
    initar <- c(list(par = init, prior = prior, lh = lh), ar,
                list(pMassProb = pMassProb, normPi0 = normPi0,
                     pMass = pMass))
    inittest <- do.call("dpost", initar) 
    if(is.infinite(inittest))
        stop("density is zero at initial parameter values")

    ar <- c(list(n = n, init = init, prior = prior, lh = lh), ar,
            list(psd = psd, thin = thin, burn = burn,
                 pMassProb = pMassProb, normPi0 = normPi0,
                 xitilde = xitilde, pMass = pMass, cv = cv))
    do.call("gibbs.mix", ar)
}

# BETA AND GAMMA INFORMATION

"ibeta" <-
function(mean, var, shape1, shape2)
{
    
    if(missing(shape1) && missing(shape2)) {
        if(min(mean) <= 0 || max(mean) >= 1)
            stop("`mean' must contain values in (0,1)")
        if(min(var) <= 0 || max(var) >= 0.25)
            stop("`var' must contain values in (0,0.25)")
        
        n <- max(length(mean), length(var))
        mean <- rep(mean, length.out=n)
        var <- rep(var, length.out=n)

        tmp <- var^(-1) * mean * (1-mean) - 1
        shape1 <- tmp * mean 
        shape2 <- tmp * (1 - mean)
        nind <- (tmp <= 0)
    }

    if(missing(mean) && missing(var)) {
        if(min(shape1) <= 0)
            stop("`shape1' must contain positive values")
        if(min(shape2) <= 0)
            stop("`shape2' must contain positive values")
        
        n <- max(length(shape1),length(shape2))
        shape1 <- rep(shape1, length.out=n)
        shape2 <- rep(shape2, length.out=n)
        
        mean <- shape1 / (shape1 + shape2)
        var <- mean * (1 - mean) / (shape1 + shape2 + 1)
        nind <- rep(FALSE, n)
    }
    
    mode <- ifelse(shape1 > 1 & shape2 > 1,
                   (shape1 - 1) / (shape1 + shape2 - 2), NA)
    rv <- matrix(c(shape1,shape2,mean,var,mode), ncol = 5,
              dimnames = list(1:n,c("shape1","shape2","mean","var","mode")))
    rv[nind, ] <- NA
    drop(rv)
}

"igamma" <-
function(mean, var, shape, scale)
{
    if(missing(shape) && missing(scale)) {
        if(min(mean) <= 0)
            stop("`mean' must contain positive values")
        if(min(var) <= 0)
            stop("`var' must contain positive values")
        
        n <- max(length(mean),length(var))
        mean <- rep(mean, length.out=n)
        var <- rep(var, length.out=n)

        shape <- mean^2 / var
        scale <- var / mean
    }
    
    if(missing(mean) && missing(var)) {
        if(min(shape) <= 0)
            stop("`shape' must contain positive values")
        if(min(scale) <= 0)
            stop("`scale' must contain positive values")
        
        n <- max(length(shape),length(scale))
        shape <- rep(shape, length.out=n)
        scale <- rep(scale, length.out=n)
        
        mean <- shape * scale
        var <- shape * scale^2
    }
    mode <- ifelse(shape > 1, scale * (shape - 1), NA)
    rv <- matrix(c(shape,scale,mean,var,mode), ncol = 5,
           dimnames = list(1:n, c("shape","scale","mean","var","mode")))
    drop(rv)
}

# COMPUTE POSTERIOR QUANTILES IN UPPER TAIL

"mc.quant" <-
  function(post, p, lh = c("gev", "gpd"))
{
  nc <- length(p)
  nr <- nrow(post)
  dn <- list(rownames(post), p)
  mat <- matrix(0, ncol = nc, nrow = nr, dimnames = dn)
  loc <- post[,"mu"]
  scale <- post[,"sigma"]
  shape <- post[,"xi"]
  
  if (lh == "gev")
    y <- -log(p)
  else
    y <- 1-p

  for(i in 1:nc)
    mat[,i] <- ifelse(shape,
                      loc + scale * (y[i]^(-shape) - 1)/shape,
                      loc - scale * log(y[i]))
    
  drop(mat)
}

# POSTERIOR RETURN LEVEL PLOT

"rl.pst" <-
  function(post, npy, lh = c("gev", "gpd"), ci = 0.9, lty = c(2,1),
           col = c(2,1), xlab = "return period",
           ylab = "return level",  ...)
  {
    if (missing(npy) && lh == "gpd")
      stop("``npy'' should be present with a ``gpd'' likelihood")

    if (lh == "gev")
      npy <- 1
    
    rps <- c(1/npy + 0.001, 10^(seq(0,4,len=20))[-1])
    p.upper <- 1 - 1/(npy * rps)
    mat <- mc.quant(post = post, p = p.upper, lh = lh) 
    mat <- t(apply(mat, 2, quantile, probs = c((1-ci)/2, 0.5, (1+ci)/2)))
    matplot(rps, mat, log = "x", type = "l",
            xlab = xlab, ylab = ylab, lty = lty, col = col, ...)
    invisible(list(x = rps, y = mat))
  }

"rl.pred" <-
  function(post, qlim, npy, lh = c("gev", "gpd"), period = 1, lty = 1,
           col = 1, xlab = "return period",
           ylab = "return level", ...)
{

  if (missing(npy) && lh == "gpd")
    stop("``npy'' should be present with a ``gpd'' likelihood")
  
  if (lh == "gev")
    npy <- 1
    
  np <- length(period)
  p.upper <- matrix(0, nrow = 25, ncol = np)
  qnt <- seq(qlim[1], qlim[2], length = 25)
  
  for(i in 1:25) {
    p <- (qnt[i] - post[,"mu"])/post[,"sigma"]
    p <- ifelse(post[,"xi"],
                exp( - pmax((1 + post[,"xi"] * p),0)^(-1/post[,"xi"])),
                exp(-exp(-p)))
    for(j in 1:np)
      p.upper[i,j] <- 1-mean(p^period[j])
  }

  if (lh == "gpd")
    p <- 1 + log(p)
  
  if(any(p.upper == 1))
    stop("lower q-limit is too small")
  if(any(p.upper == 0))
    stop("upper q-limit is too large")
  matplot(1/(npy * p.upper ), qnt, log = "x", type = "l", lty = lty,
          col = col, xlab = xlab, ylab = ylab, ...)
  invisible(list(x = 1/(npy * p.upper ), y = qnt))
}





