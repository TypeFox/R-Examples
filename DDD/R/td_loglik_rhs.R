td_loglik_rhs = function(t,x,pars)
{
  # number of probabilities
  lp = pars[length(pars)]
  
  # number of sigmas
  lrs = length(x) - lp 
  
  la = pars[1]
  mu = pars[2]
  K = pars[3]
  
  nn = -1:lp  
  lnn = length(nn)
  
  p = c(0,x[1:lp],0)
  sigdiv = x[(lp + 1):(lp + lrs)]
  
  lavec = pmax(rep(0,lnn),la - (la - mu)/K * nn)
  muvec = mu * rep(1,lnn)  
  dp = lavec[(2:(lp+1))-1] * nn[(2:(lp + 1))-1] * p[(2:(lp + 1))-1] + muvec[(2:(lp+1))+1] * nn[(2:(lp+1))+1] * p[(2:(lp+1))+1] - (lavec[(2:(lp+1))] + muvec[(2:(lp+1))]) * nn[(2:(lp+1))] * p[(2:(lp+1))]
  
  mutd = rep(mu,lrs)
  En = sum((0:(lp - 1)) * x[1:lp] )
  dsigdiv = mutd / En
  
  return(list(c(dp,dsigdiv)))  
}