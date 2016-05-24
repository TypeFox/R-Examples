sYamKrige = function(newCor,cmat,Obs,depVar="value",nmax,model,inv,ikri){
# Function for predicting at one location
  if (!missing(ikri) && ikri %% 100 == 0) print(ikri)
  
  c0dist = spDistsN1(coordinates(Obs),newCor)
  clen = length(c0dist)
  if (nmax < clen) clen = nmax
  iobs = order(c0dist)[1:clen]
  c0dist = c0dist[iobs]
  if (!inv) {
    cmat = cmat[iobs,iobs]
    dl = dim(cmat)[1]
    cmat = rbind(cmat,rep(1,dl))
    cmat = cbind(cmat,rep(1,dl+1) )
    dl = dl+1
    cmat[dl,dl] = 0
    cinv = solve(cmat)
  } else cinv = cmat
  c0arr = variogramLine(model,dist_vector = c0dist)$gamma
  c0arr[(clen+1)] = 1

  ww = cinv %*% c0arr
  if (min(ww[1:clen]) < 0) {
    ww[1:clen] = ww[1:clen]-min(ww[1:clen])
    ww = ww/sum(ww[1:clen])
  }
  var1.pred = ww[1:clen] %*% as.matrix(Obs[iobs,depVar]@data)
  var1.ok = t(ww) %*% c0arr
  var1.var = t(ww[1:clen]) %*% ((as.matrix(Obs[iobs,depVar]@data)-as.numeric(var1.pred))^2)
  return(c(var1.pred=var1.pred,var1.ok=var1.ok,var1.var = var1.var))
}



yamamotoKrige = function(formula,Obs, newPoints,model,nsim = 0,nmax = 20) {
  depVar=as.character(formula[[2]])
 
  if (nsim >0) {
    mSim = newPoints
    for (i in 1:nsim) {
      cat(paste("Conditional simulation ",i,"\n"))
#  cSim = condSim(Obs,newPoints,isim=i,...)
      cSim = condSimYama(Obs,newPoints,isim=i,model = model,depVar=depVar,nmax = nmax)
      if (i ==1) mSim = cSim  else   mSim = cbind(mSim,cSim@data)
    } 
    names(mSim)[-(1:2)] = mapply(FUN = function(i) paste("sim",i,sep=""),seq(1:nsim))
    coordinates(mSim) = as.formula(paste("~",names(mSim)[1],"+",names(mSim)[2]))
    return(mSim)
  }
  cObs = coordinates(Obs)
  if (is(newPoints,"Spatial")) {
    cNew = coordinates(newPoints)
  } else {
    cNew = newPoints
  }
  if (is.null(dim(cNew))) cNew = as.data.frame(cNew)
#  if (dim(cNew)[2] == 1) cNew = t(cNew)
  dvar = as.matrix(dist(cObs))
  cmat = variogramLine(model, dist_vector = dvar)
  dl = dim(cmat)[1]
  cmat = rbind(cmat,rep(1,dl))
  cmat = cbind(cmat,rep(1,dl+1) )
  dl = dl+1
  cmat[dl,dl] = 0
  if (nmax >= dim(Obs)[1]) {
    cmat = solve(cmat)
    inv = TRUE
  } else inv = FALSE
  ikri = c(1:dl) 
  preds = t(apply(cNew,MARGIN = 1,FUN = sYamKrige,model = model, Obs = Obs, depVar=depVar,
    cmat = cmat,nmax = nmax,inv = inv,ikri=ikri))
  preds = as.data.frame(cbind(preds,cNew))
  coordinates(preds) = as.formula(paste("~",dimnames(cNew)[2][[1]][1],"+",dimnames(cNew)[2][[1]][2]))
#  print(preds)
  return(preds)
}



condSimYama = function(Obs,newPoints,isim=1,model,depVar="value",nmax = 25) {
  if (length(names(Obs)) == 1) depVar = names(Obs)
  ccObs = coordinates(Obs)
  nObs = dim(ccObs)[1]
#  dataObs = matrix(nrow = nObs,ncol = ssim)
#  dataObs = t(mapply(Obs@data[[val]],FUN = function(X,ssim) c(rep(X,ssim)),MoreArgs = list(ssim)))
  sim = Obs@data[[depVar]]
  coords = coordinates(newPoints)
  nPred = dim(coords)[1]
  iord = sample(c(1:nPred))
#  nPred = 100

  cObs = matrix(ncol = 2,nrow=(nObs+nPred))
  cObs[1:nObs,] = ccObs
  iObs = nObs
  for (inew in 1:nPred) {
    if ((iObs-nObs+1) %% 100 == 0) cat(paste("Simulation number",isim,"Simulating point",iObs-nObs+1,"of",nPred,"\n"))
    cNew = coords[iord[inew],]
    c0dist = spDistsN1(cObs[1:iObs,],cNew)
    clen = length(c0dist)
    if (nmax < iObs) {
      jord = order(c0dist)
      c0dist = c0dist[jord[1:nmax]]
      rObs = cObs[jord[1:nmax],]
      dObs = sim[jord[1:nmax]]
    } else {
      rObs = cObs[1:nObs]
      dObs = sim
    }
    c0arr = variogramLine(model,dist_vector = c0dist)$gamma
    rvar = as.matrix(dist(rObs))
    cmat = variogramLine(model,dist_vector = rvar)
    dl = dim(cmat)[1]
    clen = dl
    cmat = rbind(cmat,rep(1,dl))
    cmat = cbind(cmat,rep(1,dl+1) )
    dl = dl+1
    cmat[dl,dl] = 0
    cinv = solve(cmat)
    c0arr[dl] = 1

    ww = cinv %*% c0arr
    if (min(ww[1:clen]) < 0) {
      ww[1:clen] = ww[1:clen]-min(ww[1:clen])
      ww = ww/sum(ww[1:clen])
    }
    var1.pred = ww[1:clen] %*% as.matrix(dObs)
    pred.var = t(ww) %*% c0arr
    var1.var = t(ww[1:clen]) %*% ((as.matrix(dObs)-as.numeric(var1.pred))^2)
    iObs = iObs + 1
    cObs[iObs,] = cNew
#    cObs = rbind(cObs,cNew)
    sim[iObs] = rnorm(1,var1.pred,sqrt(var1.var))

  }

  sim = sim[(nObs+1):(nObs+nPred)]
  sim = sim[order(iord)]

  preds = as.data.frame(cbind(coords,sim = sim))
  coordinates(preds) = as.formula(paste("~",dimnames(as.data.frame(coords))[2][[1]][1],"+",dimnames(as.data.frame(coords))[2][[1]][2]))
#  print(preds)
  return(preds)
}

