Metropolis = function(loglikelihood, sigma, m1, niter, gen, logproposal, logprior = function(x) 0, burn = 0, save_int = 10, verbose = TRUE, ...){
  # logprior: function(x) given prior opinions on model likelihood
  # loglikelihood: function(x, ...) giving likelihood of model (similar to misfit)
  # sigma: standard deviations of model elements used for calculating proposal distribution--exp(-0.5*sum((x-y)^2/sigma^2))
  # m1: starting model
  # niter: number of iterations to run
  # burn: length of burn-in period
  # save_int: number of iterations between saving models
  # ...: additional parameters to feed to likelihood

  if(length(sigma) != length(m1)){
    stop('length(sigma) != length(m1)')
  }

  # posterior distribution function
  logposterior = function(x, ...){
    p = logprior(x)
    l = loglikelihood(x, ...)
    l + p
  }
  
  if(missing(logproposal)){
    logproposal = function(x, y, sigma){
      -0.5 * sum((x - y)^2/sigma^2)
    }
  }

  # test initial model to make sure it has nonzero likelihood
  if(is.na(logposterior(m1, ...)) || logposterior(m1, ...) == -Inf){
    stop('m1 is not allowed: posterior distribution value is either 0 or NaN')
  }
  # initialize everything
    
  best = list(l = -Inf)
  acceptance = rep(0, 100)
  nsave = sum((1:niter >= burn) & (1:niter/save_int == floor(1:niter/save_int)))
  n = length(sigma)
  modmat = matrix(0, nsave, n)
  modmat_nosave = matrix(0, niter, n)
  l = rep(0, nsave)
  a = rep(0, nsave)
  q1 = logposterior(m1, ...) 
  j = 0
  
  # iterate niter times
  for(i in 1:niter){
    q0 = q1
    m2 = gen(m1, sigma)

    q2 = logposterior(m2, ...) # note that s, alpha, are logs

    if(is.na(q2)){
      print('Error: Metropolis crashed on proposed model:')
      print(m2)
      stop('logposterior of model is NaN')
    }

    r1 = logproposal(m1, m2, sigma) # p(jumping from m1 to m2)
    r2 = logproposal(m2, m1, sigma)
    alpha = min(0, (q2 + r2) - (q1 + r1))
    # if q(m2) >= q(m1), definitely jump.  otherwise, possibly jump.
    if(log(runif(1)) < alpha){
      action = 'accept'
      m1 = m2
      q1 = q2
      acceptance = c(acceptance[2:100], 1)
    }else{
      action = 'reject'
      acceptance = c(acceptance[2:100], 0)
    }
    modmat_nosave[i, ] = m1

    # print updates to screen
    if(i/1000 == floor(i/1000)){
      uncentered_cor = function(x, y)sum(x*y)/sqrt(sum(x^2)*sum(y^2)+1e-100)
      if(verbose){
          print(paste(i, 'of', niter, ';',  signif(q0, 4), ';', signif(q2, 4), ';', signif(best$l, 4), ';', action, ';', sum(acceptance)/min(100, i)))#, ';', mean(apply(modmat_nosave[max(i-1000,1):max(i-1,2),], 1, uncentered_cor, m1))))
      }
    }    
    # save result
    if(i >= burn && i/save_int == floor(i/save_int)){
      j = j + 1
      l[j] = q1
      modmat[j,] = m1
      a[j] = sum(acceptance)/min(i, 100)
    }
    # update the best result even if this isn't a 'save' iteration
    if(q1 > best$l){
      best$l = q1
      best$mod = m1
    }

  }
  return(list(m = modmat, l = l, a = a, best = best))
}
