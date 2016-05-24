### functions for computing the many means estimates- Stephen Reid
### returns 
###   i) selected indices
###   ii) selection adjusted point estimates
###   iii) selection adjusted interval estimates
###   iv) selection adjusted p-value of hypothesis testing whether underlying signal is 0

#########################
##### MAIN FUNCTION #####
#########################

#### user-facing function for computing
#### selected set
#### point and interval estimates
#### p-values
#### input:
####    - y = Vector of observations
####    - alpha = Significance level used in CI construction
####    - bh.q = q parameter for BH(q) procedure (default: NULL)
####    - k = Number of largest elements to consider (default: NULL)
####    - sigma = Estimate of standard deviation of one of the components
#### output:
####    * A list (of class "mm") with the following components :
####        - mu.hat = Vector of length length(y) containing the estimated signal size. If a sample element is not selected, then its signal size estimate is 0
####        - selected.set = Indices into the vector y of the sample elements that were selected by our procedure (either BH(q) or top-K)
####        - CIs = Matrix with two columns and number of rows equal to number of elements in selected.set. Provides the post-selection CI bounds for the estimated signal sizes of selected elements. CIs given is rows in the same order as encountered in selected.set
####        - p.vals = Vector of p-values for the test of nullity of the signals of the selected sample elemetns. P-values given in the same order as selected.set

manyMeans <- function(y, alpha=0.1, bh.q=NULL, k=NULL, sigma=1, verbose=FALSE) {
  this.call = match.call()
  if (missing(y) || is.null(y)) stop("y must be specified")
  if (!is.numeric(y)) stop("y must be numeric")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
  if (is.null(bh.q) && is.null(k)) stop("You must either bh.q or k; they cannot both be NULL")
  if (!is.null(bh.q) && (bh.q <= 0 || bh.q >= 1)) stop("bh.q must be between 0 and 1")
  if (!is.null(k) && (k < 1 || k > length(y) || k != round(k))) stop("k must be an integer between 1 and length(y)")
  if (sigma <= 0) stop("sigma must be > 0")
  
  n = length(y)
  if (!is.null(bh.q)) { # use BH selection procedure
    
    if (verbose && !is.null(k)) cat("(Both bh.q and k have been specified; k is being ignored)\n")
    k = NULL
    ci=NULL
    ### find the selected set and threshold
    p.vals = 2*pnorm (abs(y)/sigma, 0, 1, lower.tail=FALSE)
    order.p.vals = order(p.vals)
    sorted.p.vals = p.vals[order.p.vals]
    
    options (warn=-1) # ignore warning if max is over empty set
    last.reject = max(which (sorted.p.vals <= bh.q*(1:n)/n))
    options (warn=0) # reinstitute warnings
    
    if (last.reject == -Inf){ # none rejected
      if (verbose) cat("No sample elements selected.\n")
      out = list(mu.hat=rep(0,n), selected.set=NULL, pv=NULL, ci=NULL, method="BH(q)",
        bh.q=bh.q, k=NULL, threshold=NULL, sigma=sigma, call=this.call)
      class(out) = "manyMeans"
      return(out)
    }
    
    selected.set = order.p.vals[1:n <= last.reject]
    threshold = sigma*qnorm (last.reject/2/n, lower.tail=FALSE)
  }

  else{ # use top-k selection procedure
    
    ### find the selected set and threshold
    if (k == n) { # make no changes - return MLE
      z.alpha = qnorm (alpha/2, lower.tail=FALSE)
      cis = cbind(y - z.alpha*sigma, y + z.alpha*sigma)
      p.vals = 2*pnorm (abs(y), 0, sigma, lower.tail=FALSE)
      
      out = list(mu.hat=y, selected.set=1:n, pv=p.vals, ci=ci, method="top-K",
        bh.q=NULL, k=k, threshold=NULL, sigma=sigma, call=this.call)
      class(out) = "manyMeans"
      return(out)
    }
    
    order.abs.y = order (-abs(y))
    sorted.abs.y = y[order.abs.y]
    
    selected.set = order.abs.y[1:k]
    threshold = abs(sorted.abs.y[k+1])
  }
  
  ### estimate their underlying signal sizes
  mu.hat = sapply (selected.set, function(s){
    uniroot(f=function(mu){tn.mean(mu, -threshold, threshold, sigma=sigma) - y[s]}, lower=-10000*sigma, upper=10000*sigma)$root
  })
  
  ### and CIs
  right.ci = sapply (selected.set, function(s){
    uniroot (f=function(mu){tn.cdf (y[s], mu, -threshold, threshold, sigma=sigma) - (alpha/2)}, lower=-10000*sigma, upper=10000*sigma)$root
  })
  left.ci = sapply (selected.set, function(s){
    uniroot (f=function(mu){tn.cdf (y[s], mu, -threshold, threshold, sigma=sigma) - (1-alpha/2)}, lower=-10000*sigma, upper=10000*sigma)$root
  })
  
  ### and p-values
  p.vals = sapply (selected.set, function(s){tn.cdf (y[s], 0, -threshold, threshold, sigma=sigma)})
  p.vals = 2*pmin(p.vals, 1-p.vals)
  
  ### arrange
  order.selected.set = order (selected.set)
  selected.set = selected.set[order.selected.set]
  mu.hat = mu.hat[order.selected.set]
  left.ci = left.ci[order.selected.set]
  right.ci = right.ci[order.selected.set]
  p.vals = p.vals[order.selected.set]
  
  mu.hat.final = rep(0, n)
  mu.hat.final[selected.set] = mu.hat
  
  out = list(mu.hat=mu.hat.final, selected.set=selected.set, pv=p.vals, ci=cbind(left.ci,right.ci),
    method=ifelse(is.null(bh.q), "top-K", "BH(q)"), sigma=sigma, bh.q=bh.q, k=k, threshold=threshold,
    call=this.call)
  class(out) = "manyMeans"
  return(out)
}

#### prints a pretty data frame summarising the information of an object of the mm class
#### columns for index, signal size estimate, left and right CI bounds and p values
#### only for those sample elements selected by the selection procedure associated with the mmObj
print.manyMeans <- function(x, ...){
  cat("\nCall:\n")
  dput(x$call)

  cat(sprintf("\nStandard deviation of noise sigma = %0.3f\n\n",
              x$sigma))
  
  tab = cbind(x$selected.set,x$mu.hat[x$selected.set],x$pv,x$ci)
  tab = round(tab,3)
  colnames(tab) = c("SelInd","MuHat","P-value","LowConfPt","UpConfPt")
  rownames(tab) = rep("",nrow(tab))
  print(tab)
}

###############################
##### AUXILIARY FUNCTIONS #####
###############################

#### function returning the cumulative distribution function value
#### of a truncated Gaussian RV, truncated to interval (-Inf, a) \union (b, Inf)
#### with underlying Gaussian having mean parameter mu and standard deviation sigma
#### at value y
tn.cdf = function(y, mu, a, b, sigma=1){
  ## denominator
  d_right = pnorm (b, mu, sigma, lower.tail=FALSE, log.p=TRUE)
  d_left = pnorm (a, mu, sigma, lower.tail=TRUE, log.p=TRUE)
  d_max = max(d_right, d_left)
  d_log = d_max + log(exp(d_left - d_max) + exp(d_right - d_max))
  
  
  # numerator
  if (y > a & y < b){
    n_log = d_left
    return (exp(n_log-d_log))
  }else{
    if (y > b){
      # b and y
      n_y_tilde = pnorm (y, mu, sigma, lower.tail=FALSE, log.p=TRUE)
      n_b_tilde = pnorm (b, mu, sigma, lower.tail=FALSE, log.p=TRUE)
      n_yb = n_b_tilde + log(1 - exp(n_y_tilde-n_b_tilde))
      
      # a
      n_a = d_left
      
      # combine
      return(exp(n_yb-d_log) + exp(n_a-d_log))
    }else{
      n_log = pnorm (y, mu, sigma, lower.tail=TRUE, log.p=TRUE)
      return (exp(n_log-d_log))
    }
  }
}

##### function for computing the mean of an N(mu, 1) RV
##### truncated to be on the interval (-Inf, a) \union (b, Inf)
tn.mean = function(mu, a, b, sigma=1){
  # denominator
  d_left = pnorm (a, mu, sigma, lower.tail=TRUE, log.p=TRUE)
  d_right = pnorm (b, mu, sigma, lower.tail=FALSE, log.p=TRUE)
  d_max = max(d_left, d_right)
  d_log = d_max + log(exp(d_left - d_max) + exp(d_right - d_max))
  
  # numerator
  n_left = dnorm (b, mu, sigma, log=TRUE)
  n_right = dnorm (a, mu, sigma, log=TRUE)
  
  if (n_left > n_right){
    mu + exp(n_left + log(1 - exp(n_right-n_left)) - d_log)
  }else{
    mu - exp(n_right + log(1 - exp(n_left-n_right)) - d_log)
  }
}
