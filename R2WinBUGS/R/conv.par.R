"conv.par" <-
function (x, n.chains, Rupper.keep = TRUE) {
  m <- ncol(x)
  n <- nrow(x)

# We compute the following statistics:
#
#  xdot:  vector of sequence means
#  s2:  vector of sequence sample variances (dividing by n-1)
#  W = mean(s2):  within MS
#  B = n*var(xdot):  between MS.
#  muhat = mean(xdot):  grand mean; unbiased under strong stationarity
#  varW = var(s2)/m:  estimated sampling var of W
#  varB = B^2 * 2/(m+1):  estimated sampling var of B
#  covWB = (n/m)*(cov(s2,xdot^2) - 2*muhat*cov(s^2,xdot)):
#                                               estimated sampling cov(W,B)
#  sig2hat = ((n-1)/n))*W + (1/n)*B:  estimate of sig2; unbiased under
#                                               strong stationarity
#  quantiles:  emipirical quantiles from last half of simulated sequences

  xdot <- apply(x,2,mean)
  muhat <- mean(xdot)
  s2 <- apply(x,2,var)
  W <- mean(s2)
  quantiles <- quantile (as.vector(x), probs=c(.025,.25,.5,.75,.975))

  if ((W > 1.e-8) && (n.chains > 1)) {            # non-degenerate case
  
  B <- n*var(xdot)
  varW <- var(s2)/m
  varB <- B^2 * 2/(m-1)
  covWB <- (n/m)*(var(s2, xdot^2) - 2*muhat*var(s2, xdot))
  sig2hat <- ((n-1)*W + B)/n

# Posterior interval post.range combines all uncertainties
# in a t interval with center muhat, scale sqrt(postvar),
# and postvar.df degrees of freedom.
#
#       postvar = sig2hat + B/(mn):  variance for the posterior interval
#                               The B/(mn) term is there because of the
#                               sampling variance of muhat.
#       varpostvar:  estimated sampling variance of postvar

    postvar <- sig2hat + B/(m*n)
    varpostvar <- max(0, 
        (((n-1)^2) * varW + (1 + 1/m)^2 * varB + 2 * (n-1) * (1 + 1/m) * covWB) / n^2)
    post.df <- min(2*(postvar^2/varpostvar), 1000)

# Estimated potential scale reduction (that would be achieved by
# continuing simulations forever) has two components:  an estimate and
# an approx. 97.5% upper bound.
#
# confshrink = sqrt(postvar/W),
#     multiplied by sqrt(df/(df-2)) as an adjustment for the
###      CHANGED TO sqrt((df+3)/(df+1))
#     width of the t-interval with df degrees of freedom.
#
# postvar/W = (n-1)/n + (1+1/m)(1/n)(B/W); we approximate the sampling dist.
# of (B/W) by an F distribution, with degrees of freedom estimated
# from the approximate chi-squared sampling dists for B and W.  (The
# F approximation assumes that the sampling dists of B and W are independent;
# if they are positively correlated, the approximation is conservative.)

    confshrink.range <- postvar/W
    if(Rupper.keep){
        varlo.df <- 2*(W^2/varW) 
        confshrink.range <- c(confshrink.range, 
            (n-1)/n + (1+1/m)*(1/n)*(B/W) * qf(.975, m-1, varlo.df))
    }
    confshrink.range <- sqrt(confshrink.range * (post.df+3) / (post.df+1))
    
# Calculate effective sample size:  m*n*min(sigma.hat^2/B,1)
# This is a crude measure of sample size because it relies on the between
# variance, B, which can only be estimated with m degrees of freedom.
    
    n.eff <- m*n*min(sig2hat/B,1)
    list(quantiles=quantiles, confshrink=confshrink.range,
         n.eff=n.eff)

  }
  else {      # degenerate case:  all entries in "data matrix" are identical
    list (quantiles=quantiles, confshrink = rep(1, Rupper.keep + 1),
          n.eff=1)

  }
}
