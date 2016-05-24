VarInflation <-
function(dta,Blist,maxnbfactors,dig) {
  m = ncol(dta)
  n = nrow(dta)
  vecrho = round(seq(10^(-dig),1,10^(-dig)),digits=dig)
  vecdt = unlist(lapply(vecrho,Dt))
  sampled = sample(1:m,min(500,m))
  sampsize = length(sampled)
  cordata = t(dta[,sampled])%*%dta[,sampled]/(n-1)
  sdt = rep(0,maxnbfactors+1)
  names(sdt) = paste(0:maxnbfactors,"factors")
  for (i in 1:(maxnbfactors+1)) {
    #    print(paste("Calculating criterion for the model with",i-1,"factors"))
    B = matrix(Blist[[i]][sampled,],nrow=sampsize)
    sdb = sqrt(1-apply(B^2,1,sum))   ################################### NaNs are generated !!!
    matrho = cordata - B%*%t(B)
    matrho = sweep(matrho,2,FUN="/",STATS=sdb)
    matrho = sweep(matrho,1,FUN="/",STATS=sdb)
    rho = matrho[col(matrho)>row(matrho)]
    rho[abs(rho)>=1] = 1
    veccor = sort(round(abs(rho),digits=dig))
    duplic = duplicated(veccor)
    vduplic = sort(unique(veccor[duplic]))
    vunic = setdiff(unique(veccor),vduplic)
    dtunic = vecdt[is.element(vecrho,vunic)]
    dtduplic = vecdt[is.element(vecrho,vduplic)]
    vmatch = match(vecrho,veccor,0)
    nboccur = diff(c(vmatch[vmatch>0],length(veccor)+1))
    nboccur = nboccur[nboccur>1]
    sdt[i] = sum(dtunic)+crossprod(nboccur,dtduplic)  }
  return(sdt) }
