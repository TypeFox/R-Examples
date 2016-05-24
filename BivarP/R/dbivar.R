dibivar = function(x, y, par, afa, rodina, fam){
  # require(copula)
  nx = length(x)
  ny = length(y)
  dxy = array(0,c(ny, nx))
  colnames(dxy) = x
  rownames(dxy) = y
  if(rodina[1]=="weibull"){
    ux = pweibull(x, shape=par[1], scale=par[2])
    dux = dweibull(x, shape=par[1], scale=par[2])
  } else if(rodina[1]=="gamma"){
    ux = pgamma(x, shape=par[1], scale=par[2])
    dux = dgamma(x, shape=par[1], scale=par[2])
  } else if(rodina[1]=="lnorm"){
    ux = plnorm(x, meanlog=par[1], sdlog=par[2])
    dux = dlnorm(x, meanlog=par[1], sdlog=par[2])
  } else {
    # rodina[1]=="norm
    ux = pnorm(x, mean=par[1], sd=par[2])
    dux = dnorm(x, mean=par[1], sd=par[2])
  }
  if(rodina[2]=="weibull"){
    vy = pweibull(y, shape=par[3], scale=par[4])
    dvy = dweibull(y, shape=par[3], scale=par[4])
  } else if(rodina[2]=="gamma"){
    vy = pgamma(y, shape=par[3], scale=par[4])
    dvy = dgamma(y, shape=par[3], scale=par[4])
  } else if(rodina[2]=="lnorm"){
    vy = plnorm(y, meanlog=par[3], sdlog=par[4])
    dvy = dlnorm(y, meanlog=par[3], sdlog=par[4])
  } else {
    # rodina[2]=="norm
    vy = pnorm(y, mean=par[3], sd=par[4])
    dvy = dnorm(y, mean=par[3], sd=par[4])
  }
  #oooooooooooooooooooooooooooooooooo
  kop = archmCopula(fam, afa, dim = 2)
  cxy = array(0,c(ny, nx))
  for(i in 1:nx){
    yx = cbind(vy,rep(ux[i],ny))
    cxy[, i] = dCopula(yx, kop)
  }
  #oooooooooooooooooooooooooooooooooo
  for(i in 1:nx){
    for(j in 1:ny){
      dxy[j, i] = cxy[j, i]*dvy[j]*dux[i]
    }
  }
  return(dxy)
}
#=======================================
pibivar = function(x, y, par, afa, rodina, fam){
#  require(copula)
  nx = length(x)
  ny = length(y)
  pxy = array(0,c(ny, nx))
  colnames(pxy) = x
  rownames(pxy) = y
  if(rodina[1]=="weibull") ux = pweibull(x, shape=par[1], scale=par[2]) else {
    if(rodina[1]=="gamma") ux = pgamma(x, shape=par[1], scale=par[2]) else {
      if(rodina[1]=="lnorm") ux = plnorm(x, meanlog=par[1], sdlog=par[2]) else {
        # rodina[1]=="norm"
        ux = pnorm(x, mean=par[1], sd=par[2])
      }}}
  #=============================
  if(rodina[2]=="weibull") vy = pweibull(y, shape=par[3], scale=par[4]) else {
    if(rodina[2]=="gamma") vy = pgamma(y, shape=par[3], scale=par[4]) else {
      if(rodina[2]=="lnorm") vy = plnorm(y, meanlog=par[3], sdlog=par[4]) else {
        # rodina[2]=="norm"
        vy = pnorm(y, mean=par[3], sd=par[4])
      }}}
  #oooooooooooooooooooooooooooooooooo
  kop = archmCopula(fam, afa, 2)
  pxy = array(0,c(ny, nx))
  for(i in 1:nx){
    yx = cbind(vy,rep(ux[i],ny))
    pxy[, i] = pCopula(yx, kop)
  }
  #oooooooooooooooooooooooooooooooooo
  return(pxy)
}
#=========================================
subivar = function(x, y, par, afa, rodina, fam){
#  require(copula)
  nx = length(x)
  ny = length(y)
  sxy = array(0,c(ny, nx))
  colnames(sxy) = x
  rownames(sxy) = y
  if(rodina[1]=="weibull") ux = pweibull(x, shape=par[1], scale=par[2], lower.tail=FALSE) else {
    if(rodina[1]=="gamma") ux = pgamma(x, shape=par[1], scale=par[2], lower.tail=FALSE) else {
      if(rodina[1]=="lnorm") ux = plnorm(x, meanlog=par[1], sdlog=par[2], lower.tail=FALSE) else {
        # rodina[1]=="norm"
        ux = pnorm(x, mean=par[1], sd=par[2], lower.tail=FALSE)
      }}}
  if(rodina[2]=="weibull") vy = pweibull(y, shape=par[3], scale=par[4], lower.tail=FALSE) else {
    if(rodina[2]=="gamma") vy = pgamma(y, shape=par[3], scale=par[4], lower.tail=FALSE) else {
      if(rodina[2]=="lnorm") vy = plnorm(y, meanlog=par[3], sdlog=par[4], lower.tail=FALSE) else {
        # rodina[2]=="norm"
        vy = pnorm(y, mean=par[3], sd=par[4], lower.tail=FALSE)
      }}}
  #oooooooooooooooooooooooooooooooooo
  kop = archmCopula(fam, afa, 2)
  sxy = array(0,c(ny, nx))
  for(i in 1:nx){
    yx = cbind(vy,rep(ux[i],ny))
    sxy[, i] = pCopula(yx, kop)
  }
  #oooooooooooooooooooooooooooooooooo
  return(sxy)
}
