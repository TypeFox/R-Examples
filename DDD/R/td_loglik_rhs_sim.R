td_loglik_rhs_sim = function (t, x, pars) 
{
  lp = pars[length(pars)]
  lrs = length(x) - lp
  la = pars[1]
  mu = pars[2]
  K = pars[3]
  nn = -1:lp
  lnn = length(nn)
  p = c(0, x[1:lp], 0)
  lavec = pmax(rep(0, lnn), la - (la - mu)/K * nn)
  muvec = mu * rep(1, lnn)
  dp = lavec[(2:(lp + 1)) - 1] * nn[(2:(lp + 1)) - 1] * p[(2:(lp + 1)) - 1] + muvec[(2:(lp + 1)) + 1] * nn[(2:(lp + 1)) +  1] * p[(2:(lp + 1)) + 1] - (lavec[(2:(lp + 1))] + muvec[(2:(lp + 1))]) * nn[(2:(lp + 1))] * p[(2:(lp + 1))]
  return(list(dp))
}