ScanCBS2dEM = function(AT, BT, AN, BN, initialTau=NULL, minGridSize=10, statistic="binomial2d", grid.size="auto", takeN=5, maxNCut=100, minStat=0, alpha=0.05, verbose=FALSE, timing=TRUE, error=1e-5, maxIter=1000, pOri=c(0.49,0.51)){
  # error, maxIter, pOri: parameters for GetP
  # initialTau=NULL; minGridSize=10; statistic="binomial2d"; grid.size="auto"; takeN=5; maxNCut=100; minStat=0; alpha=0.05; verbose=FALSE; timing=TRUE; error=1e-5; maxIter=1000; pOri=c(0.49,0.51)

  timeCBSTotal = proc.time()[3]
  timeCBSPreProcess = proc.time()[3]
  p = as.numeric(.Call("GetP", as.numeric(AT), as.numeric(BT), as.numeric(AN), as.numeric(BN), as.numeric(error), as.numeric(maxIter), as.numeric(pOri), PACKAGE="falcon"))
  likh = as.numeric(.Call("LikH", as.numeric(AT), as.numeric(BT), as.numeric(AN), as.numeric(BN), as.numeric(p), PACKAGE="falcon"))
  nL = length(AT)
  tauHat = c(1,nL+1)

  if (is.null(initialTau)){
    mBIC = ScanBIC2dEM(AT, BT, AN, BN, tauHat, error, maxIter, pOri)
  }else{
    mBIC = ScanBIC2dEM(AT, BT, AN, BN, sort(unique(c(tauHat, initialTau))), error, maxIter, pOri)
  }
  curBIC = mBIC
  if (verbose) print(paste("curBIC:", curBIC))
  if (grid.size=="auto"){
    grid.size = getAutoGridSize(nL)
    if (grid.size[length(grid.size)] >= nL){
      grid.size = grid.size[-length(grid.size)]
    }
  }else if (nL/max(grid.size)>100){
    print("Grid is too fine; may result in slow computation")
  }
  grid.size = sort(grid.size, decreasing=TRUE)
  nGridSize = length(grid.size)
  while(grid.size[nGridSize] < minGridSize){
    grid.size = grid.size[-nGridSize]
    nGridSize = length(grid.size)
  }

  min.grid.size.T2 = 2*min(grid.size)
  if (verbose) print(grid.size)

  timeCBSPreProcess = proc.time()[3] - timeCBSPreProcess

  timeIGSTotal = proc.time()[3]
  timeIGSBreakDown = matrix(0, nrow=nGridSize, ncol=3)
  colnames(timeIGSBreakDown) = c("newComp", "refineComp", "misc")
  rownames(timeIGSBreakDown) = grid.size

  maxSRes = ScanIterateGrid2dEM(AT, BT, AN, BN, statistic, grid.size, nGridSize, minGridSize, timeIGSBreakDown, takeN, verbose, timing, error, maxIter, pOri)
  if(verbose) print(maxSRes)
  maxS = matrix(c(maxSRes$cptsRet, 1, nL+1), ncol=5)
  timeIGSBreakDown = maxSRes$timeIGSBreakDown
  statHat = matrix(nrow=0, ncol=6)
  cutStat = minStat + 1
  nCut = 1

  timeIGSTotal = proc.time()[3] - timeIGSTotal
  timeCBSBreakDown = c()

  while((nCut < maxNCut) && (abs(cutStat) > minStat)){
    timeCBSnew = proc.time()[3]
    if (verbose){
      print(dim(maxS))
      print(maxS)
      print(maxS[(max(c(maxS[,2]-maxS[,1], maxS[,1]-maxS[,4], maxS[,5]-maxS[,2])) >= min.grid.size.T2),])
    }
    maxS = matrix(maxS[(max(c(maxS[,2]-maxS[,1], maxS[,1]-maxS[,4], maxS[,5]-maxS[,2])) >= min.grid.size.T2),], ncol=5)
    if (verbose) print(maxS)
    if (length(maxS)==0) break
    maxS.ind = which.max(abs(maxS[,3]))
    maxS.cut = maxS[maxS.ind,]
    if (is.null(initialTau)){
      mBIC = ScanBIC2dEM(AT, BT, AN, BN, unique(sort(c(tauHat, maxS.cut[1:2]))), error, maxIter, pOri)
    }else{
      mBIC = ScanBIC2dEM(AT, BT, AN, BN, unique(sort(c(tauHat, maxS.cut[1:2], initialTau))), error, maxIter, pOri)
    }
    # if (mBIC < curBIC && nCut != 1) break
    if (mBIC < curBIC) break
    tauHat = unique(sort(c(tauHat, maxS.cut[1:2])))
    if (verbose) print(tauHat)
    statHat = rbind(statHat, c(maxS.cut, mBIC))
    if (verbose) print(c(maxS.cut, mBIC))

    cutStat = maxS.cut[3]
    maxS.cut.L = maxS.cut[4]
    maxS.cut.C1 = maxS.cut[1]
    maxS.cut.C2 = maxS.cut[2]
    maxS.cut.R = maxS.cut[5]

    grid.size.L = grid.size[2*grid.size <= maxS.cut.C1 - maxS.cut.L]
    grid.size.M = grid.size[2*grid.size <= maxS.cut.C2 - maxS.cut.C1]
    grid.size.R = grid.size[2*grid.size <= maxS.cut.R - maxS.cut.C2]

    timeIGSNew = proc.time()[3]

    if (length(grid.size.L) > 0){
      ids = maxS.cut.L:(maxS.cut.C1-1)
      maxS.LRes = ScanIterateGrid2dEM(AT[ids], BT[ids], AN[ids], BN[ids], statistic, grid.size.L, nGridSize, minGridSize, timeIGSBreakDown, takeN, verbose, timing, error, maxIter, pOri)
      maxS.L = maxS.LRes$cptsRet
      timeIGSBreakDown = maxS.LRes$timeIGSBreakDown
      maxS.L[c(1,2)] = maxS.L[c(1,2)] + maxS.cut.L - 1
      maxS.L = c(maxS.L, maxS.cut.L, maxS.cut.C1)
      maxS = rbind(maxS, maxS.L)
    }
    if (length(grid.size.M) > 0){
      ids = maxS.cut.C1:(maxS.cut.C2-1)
      maxS.MRes = ScanIterateGrid2dEM(AT[ids], BT[ids], AN[ids], BN[ids], statistic, grid.size.M, nGridSize, minGridSize, timeIGSBreakDown, takeN, verbose, timing, error, maxIter, pOri)
      maxS.M = maxS.MRes$cptsRet
      timeIGSBreakDown = maxS.MRes$timeIGSBreakDown
      maxS.M[c(1,2)] = maxS.M[c(1,2)] + maxS.cut.C1 -1
      maxS.M = c(maxS.M, maxS.cut.C1, maxS.cut.C2)
      maxS = rbind(maxS, maxS.M)
    }
    if (length(grid.size.R) > 0){
      ids = maxS.cut.C2:(maxS.cut.R-1)
      maxS.RRes = ScanIterateGrid2dEM(AT[ids], BT[ids], AN[ids], BN[ids], statistic, grid.size.R, nGridSize, minGridSize, timeIGSBreakDown, takeN, verbose, timing, error, maxIter, pOri)
      maxS.R = maxS.RRes$cptsRet
      timeIGSBreakDown = maxS.RRes$timeIGSBreakDown
      maxS.R[c(1,2)] = maxS.R[c(1,2)] + maxS.cut.C2 -1
      maxS.R = c(maxS.R, maxS.cut.C2, maxS.cut.R)
      maxS = rbind(maxS, maxS.R)
    }
    if(verbose) {
      print("maxS after cut")
      print(maxS)
    }
    maxS = matrix(maxS[-maxS.ind,], ncol=5)
    nCut = nCut+1
    curBIC = mBIC
    
    timeIGSTotal = timeIGSTotal + proc.time()[3] - timeIGSNew
    timeCBSBreakDown = c(timeCBSBreakDown, proc.time()[3]-timeCBSnew)
  }
  if (nCut>1){
    tauHatInd = sort(tauHat)
    statHat = cbind(matrix(statHat[,1:2],nrow=nrow(statHat)), statHat)
    colnames(statHat) = c("cptL", "cptR", "cptLReadInd","cptRReadInd", "stat", "parentL", "parentR", "BIC")
  }  
  timeCBSTotal = proc.time()[3] - timeCBSTotal
  return(tauHat)
}



ScanIterateGrid2dEM = function(AT, BT, AN, BN, statistic, grid.size, nGridSize, minGridSize, timeIGSBreakDown, takeN, verbose, timing, error, maxIter, pOri){
  ## Return:
  ## Do a variable window scan statistic of the case/control processes
  timeIGSnew = proc.time()[3]
  nL = length(AT)
  p = as.numeric(.Call("GetP", as.numeric(AT), as.numeric(BT), as.numeric(AN), as.numeric(BN), as.numeric(error), as.numeric(maxIter), as.numeric(pOri), PACKAGE="falcon"))
  timeIGSBreakDown[,3] = proc.time()[3] - timeIGSnew + timeIGSBreakDown[,3]

  ## Initialize the bag of change point calls
  cpts = matrix(nrow=0, ncol=3)
  if (grid.size[length(grid.size)] > 2*minGridSize) grid.size = c(grid.size, minGridSize)
  # grid.size = c(grid.size, 1) # too fine
  
  ## Variable Windows by Grid Refinement
  for (g in 1:length(grid.size)){
    timeRow = nGridSize - length(grid.size) + g
    timeIGSnew = proc.time()[3]
    ## 1. Construct the current grid and define max.win
    grid.cur = as.numeric(seq(1, nL+1, grid.size[g]))
    if (grid.cur[length(grid.cur)] != nL+1) grid.cur = c(grid.cur, nL+1)
    if (grid.cur[length(grid.cur)] - grid.cur[length(grid.cur)-1] == 1) grid.cur = grid.cur[-(length(grid.cur)-1)]
    if (g==1){
      max.win = length(grid.cur)-1
    }else{
      max.win = 2*floor(grid.size[g-1]/grid.size[g])
    }
    max.win = as.numeric(max.win)
    timeIGSBreakDown[timeRow,3] = timeIGSBreakDown[timeRow,3] + proc.time()[3] - timeIGSnew
    if (verbose) print(paste("PassedCumSum", g))

    ## 2. Refine existing change-points in bag
    timeIGSnew = proc.time()[3]
    if (nrow(cpts)>0){
      # grid.LR = cpts %/% grid.size[g]
      for (r in 1:nrow(cpts)){
        idL = which(abs(grid.cur-cpts[r,1]) <= (grid.size[g]*max.win/2))
        idR = which(abs(grid.cur-cpts[r,2]) <= (grid.size[g]*max.win/2))
        refineRes =.Call("ScanStatRefineCompBinom2dEMC", as.numeric(AT), as.numeric(BT), as.numeric(AN), as.numeric(BN), as.numeric(error), as.numeric(maxIter), as.numeric(pOri), as.numeric(grid.cur), as.numeric(idL-1), as.numeric(idR-1), PACKAGE="falcon")
        cpts[r,] = refineRes[which.max(abs(refineRes[,3])),]
      }
    }
    if (verbose) print(paste("PassedRefineScan",g))
    timeIGSBreakDown[timeRow,2] = timeIGSBreakDown[timeRow,2] + proc.time()[3] - timeIGSnew

    ## 3. New Scan with the current grid size
    if (grid.size[g] > 4){
      timeIGSnew = proc.time()[3]
      newRes = .Call("ScanStatNewCompBinom2dEMC", as.numeric(AT), as.numeric(BT), as.numeric(AN), as.numeric(BN), as.numeric(error), as.numeric(maxIter), as.numeric(pOri), as.numeric(grid.cur), as.numeric(max.win), PACKAGE="falcon")
      if (verbose) print(dim(newRes))
      if (takeN/nrow(newRes) > 0.2){
        cpts = rbind(cpts, newRes[order(abs(newRes[,3]), decreasing=TRUE)[1:min(c(takeN, nrow(newRes)))],])
      }else{
        abs3NewRes = abs(newRes[,3])
        for (k in 1:takeN){
          maxIndex = which.max(abs3NewRes)
          abs3NewRes = abs3NewRes[-maxIndex]
          cpts = rbind(cpts, newRes[maxIndex,])
          newRes = newRes[-maxIndex,]
        }
      }
      if (verbose) print(paste("PassedNewScan", g))
      timeIGSBreakDown[timeRow,1] = timeIGSBreakDown[timeRow,1] + proc.time()[3] - timeIGSnew
    }
  }
  cptsRet = cpts[which.max(abs(cpts[,3])),]
  return(list(cptsRet=cptsRet, timeIGSBreakDown=timeIGSBreakDown))
}
    

ScanBIC2dEM = function(AT, BT, AN, BN, tauHat, error, maxIter, pOri){
  nCpts = length(tauHat)-1
  p = as.numeric(.Call("GetP", as.numeric(AT), as.numeric(BT), as.numeric(AN), as.numeric(BN), as.numeric(error), as.numeric(maxIter), as.numeric(pOri), PACKAGE="falcon"))
  likh = as.numeric(.Call("LikH", as.numeric(AT), as.numeric(BT), as.numeric(AN), as.numeric(BN), as.numeric(p), PACKAGE="falcon"))
  mBIC = likh[2]/2 - likh[1]
  for (i in 1:nCpts){
    ids = tauHat[i]:(tauHat[i+1]-1)
    p = as.numeric(.Call("GetP", as.numeric(AT[ids]), as.numeric(BT[ids]), as.numeric(AN[ids]), as.numeric(BN[ids]), as.numeric(error), as.numeric(maxIter), as.numeric(pOri), PACKAGE="falcon"))
    likh = as.numeric(.Call("LikH", as.numeric(AT[ids]), as.numeric(BT[ids]), as.numeric(AN[ids]), as.numeric(BN[ids]), as.numeric(p), PACKAGE="falcon"))
    if ((!is.na(likh[1])) && (!is.na(likh[2]))){
      mBIC = mBIC + likh[1] - likh[2]/2
    }else if (!is.na(likh[1])){
      mBIC = mBIC + likh[1]
      cat("Hassian matrix singular!", tauHat[i], tauHat[i+1], "\n")
    }else{
      cat("log-likelihood NA!", tauHat[i], tauHat[i+1], "\n")
    }
  }
  mBIC = mBIC - nCpts*log(length(AT))
  return(mBIC)
  return(0)
}
 
