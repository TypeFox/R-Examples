asht <- 
function (X, f, theta, cmat, alternative, conf.level, args) 
{
  esti <- theta(X = X, f = as.factor(f))
  out <- waldci(cmat = cmat, estp = esti$estimate, varp = esti$varest, varcor = esti$varest, alternative = alternative, conf.level = conf.level)
  return(out)
}
