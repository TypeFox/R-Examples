# Check if zmin and zmax can actually be given through the ... notation


wEst = function(zSa,acdfl,cP,pThresh,...) {
  p = acdfFindSimp(acdfl,zSa,...)
#  cat(paste("wEst",zSa,cP,a,p,"\n"))
  1./(1.+exp(-cP*(p-pThresh)))
}


wEstA = function(zArr,acdf,cP,pThresh,...) {
  wArr = c(1:length(zArr))
  nacdf = dim(acdf)[1]
  zMin = acdf[1,1]
  zmax = acdf[nacdf,1]
  zInc = acdf[2,1]-acdf[1,1]
  pArr = mapply(wEs,zArr,MoreArgs=list(zMin = zMin,zInc = zInc,nacdf = nacdf,acdf = acdf))
  wArr =   1./(1.+exp(-cP*(pArr-pThresh)))
  return(wArr) 
}

wEs = function(zm,zMin,zInc,nacdf,acdf){
    id =  as.integer(((zm-zMin)/zInc)+2)
    if (id >nacdf) id = nacdf
    if (id <1) id = 1
    acdf[id,2]
}


#        iwqselRes = iwqsel(zPred@data,acdf,thresh,pThresh)
  
iwqsel = function(zPred,acdf,thresh,pThresh,fbate,maxit = 500, cpAddLim = 0.0001, debug.level=0,...) {
  nPred = dim(zPred)[1]
  if (!missing("fbate")) {
    itte = max(which(fbate$cP !=0))
    cP = fbate$cP[itte]
    cPadd = fbate$cPadd[itte]
    ate = fbate$ate[itte]
    late = fbate$ate[itte-1]
    bate = min(abs(fbate$ate))
  } else {
    itte = 0
    cPadd = (pThresh-0.5)/0.1
    if (abs(cPadd) < 0.2) cPadd = sign(cPadd)*0.2
    late = cPadd
    ate = cPadd
    cP = -cPadd/2
#    cPadd = 2
    bate = 1
    sate = 1
    fbate = data.frame(itte = itte,cP = 0,cPadd = 0,ate = 0)
  }
  if (debug.level != 1) cat(paste("\n after",itte,round(cP,4),round(cPadd,4),round(ate,4),round(late,4),
                  round(sate,4),as.integer(nPred/50)))
  while ((itte < maxit | abs(ate) > (1/nPred+0.00001)) & abs(cPadd) > cpAddLim) {
    itte = itte + 1
    cPadd = ifelse(sign(ate) == sign(late),cPadd,-cPadd/2)
    cP = cP+cPadd
    zEst = array(c(1:2))
    ieq = which(fbate$cP == cP)
    if (length(ieq) > 0) {
      late = ate
      ate = fbate$ate[ieq[1]]
    } else {
      t0 = proc.time()[3]
      for (iaa in 1:nPred) {
         Zloc = as.numeric(zPred[iaa,])
#      wEstArr = mapply(wEst,zSa = zPred[iaa,],
#            MoreArgs=list(acdfl = acdf,cP = cP,pThresh=pThresh,...))
        wEstArr = wEstA(Zloc,acdf,cP = cP,pThresh=pThresh,...)
        if (is.list(wEstArr)) wEstArr = unlist(wEstArr)
        lzEst =  sum(wEstArr*Zloc)/sum(wEstArr)
        if (is.na(lzEst)) stop
        zEst[iaa] = lzEst
        if ((iaa %% 50) == 0) cat(paste(1))
#      if ((iaa %% 50) == 0) cat(paste(itte,round(cP,4),round(cPadd,4),"bate: ",
#            round(bate,4),round(ate,4),round(late,4),"iaa: ",iaa,nPred,
#          " lzEst ",round(lzEst,5),round(wEstArr[1],5),pThresh,"\n"))
      }
#    print(round(zEst[zEst>thresh],5))
      t1 = proc.time()[3]
      late = ate
      ate = -(quantile(zEst,pThresh)[[1]]-thresh)[[1]]
      sate = mean(I(zEst<thresh))-pThresh
    }
    if (abs(ate) < abs(bate)) {
      bate = ate
      bitte = itte
    }
    fbate[itte,] = fbate[1,]
    fbate$cP[itte] = cP
    fbate$ate[itte] = ate
    fbate$cPadd[itte] = cPadd
#    if ((abs(ate) < 0.0002 & abs(cPadd) < 0.05) | abs(cPadd) < 0.0001) itte = 1000
    if (debug.level != 1) cat(paste("\n after",itte,round(cP,4),round(cPadd,4),round(ate,4),round(late,4),
# Removed timing, to avoid check-errors
#            round(sate,4),"time: ", round(t1-t0,2),as.integer(nPred/50)))
            round(sate,4),"time: ",as.integer(nPred/50)))
  }
  return(list(fbate = fbate,zEst=zEst))
}

