
#Level 2 test: check if functions work properly with different sets of
#parameters ("vgLargeParam" in data folder are called)#
### dpqr group - dvg, ddvg, pvg, qvg, rvg, vgCalcRange, vgBreaks

test.vgL2dpqr <- function () {
  for (i in 1:nrow(testParam)) {    
    param <- testParam[i,]
    dvgReturn1 <- dvg(x = param[1], param = param)
    dvgReturn2 <- dvg(x = 1, param = param)
    ddvgReturn1 <- ddvg(x = param[1], param = param)
    ddvgReturn2 <- ddvg(x = 1, param = param)
    pvgReturn <- pvg(q = 0, param = param)
    qvgReturn <- qvg(p = 0.5, param = param)
    rvgReturn <- rvg(n =100, param = param)  
    vgCalcRangeReturn <- vgCalcRange(param = param)
    vgBreaksReturn <-vgBreaks(param = param)   
         
    checkTrue(is(dvgReturn1, "numeric"), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(dvgReturn1),1, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    
    checkEquals(length(ddvgReturn1),1, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkTrue(is(dvgReturn2, "numeric"), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(dvgReturn2),1, msg = paste(" param= ", param[1], param[2], param[3], 
      param[4]))
    checkTrue(is(ddvgReturn2, "numeric"), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(ddvgReturn2),1, msg = paste(" param= ", param[1], param[2],
      param[3], param[4]))
    checkTrue(is(pvgReturn, "numeric"), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(pvgReturn),1, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkTrue(is(qvgReturn, "numeric"), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(qvgReturn),1, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkTrue(is(rvgReturn, "numeric"), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(rvgReturn),100)
    checkTrue(is(vgCalcRangeReturn, "numeric"), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(vgCalcRangeReturn), 2, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkTrue(is(vgBreaksReturn, "list"), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(vgBreaksReturn), 7, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
  }
}

### moments group - vgMode, vgMean, vgVar, vgSkew, vgKurt, vgMom
test.vgL2moments <- function () {    
  for (i in 1:nrow(testParam)) {
    param <- testParam[i,]
    vgModeReturn <- vgMode(param = param)
    vgMeanReturn <- vgMean(param = param)
    vgVarReturn <- vgVar(param = param)
    vgSkewReturn <- vgSkew(param = param)
    vgKurtReturn <- vgKurt(param = param)
    vgMomReturnRaw8 <- vgMom (order = 8, param = param)
    vgMomReturnRaw7 <- vgMom (order = 7, param = param)
    vgMomReturnCen8 <- vgMom (order = 8, param = param, momType = "central") 
    vgMomReturnCen7 <- vgMom (order = 7, param = param, momType = "central") 
    vgMomReturnMu8 <- vgMom (order = 8, param = param, momType = "mu") 
    vgMomReturnMu7 <- vgMom (order = 7, param = param, momType = "mu") 
    vgMomReturnAbo8 <- vgMom (order = 8, param = param, about = 1)
    vgMomReturnAbo7 <- vgMom (order = 7, param = param, about = 1)
  
    checkTrue(is(vgModeReturn, "numeric"), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(vgModeReturn),1, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))  
    checkTrue(is(vgMeanReturn, "numeric"), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(vgMeanReturn),1, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))  
    checkTrue(is(vgMeanReturn, "numeric"), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(vgMeanReturn),1, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4])) 
    checkTrue(is(vgVarReturn, "numeric"), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(vgVarReturn),1, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4])) 
    checkTrue(vgVarReturn > 0, msg = paste(" param= ", param[1], param[2], param[3], 
      param[4]))
    checkTrue(is(vgSkewReturn, "numeric"), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(vgSkewReturn),1, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4])) 
    checkTrue(is(vgKurtReturn, "numeric"), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(vgKurtReturn),1, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4])) 
    checkTrue(is(vgMomReturnRaw8, "numeric"), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(vgMomReturnRaw8),1, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4])) 
    checkTrue(is(vgMomReturnRaw7, "numeric"), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(vgMomReturnRaw7),1, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4])) 
    checkTrue(is(vgMomReturnCen8, "numeric"), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(vgMomReturnCen8),1, msg = paste(" param= ", param[1], param[2],
      param[3], param[4])) 
    checkTrue(is(vgMomReturnCen7, "numeric"), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(vgMomReturnCen7),1, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4])) 
    checkTrue(is(vgMomReturnMu8, "numeric"), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(vgMomReturnMu8),1, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4])) 
    checkTrue(is(vgMomReturnMu7, "numeric"), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(vgMomReturnMu7),1, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4])) 
    checkTrue(is(vgMomReturnAbo8, "numeric"), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(vgMomReturnAbo8),1, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4])) 
    checkTrue(is(vgMomReturnAbo7, "numeric"), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(vgMomReturnAbo7),1, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4])) 
  }
}


### fitting group - vgFitStartMoM, vgFitStart, vgFit
test.vgL2fitting <- function () {
  for (i in 1:nrow(testParam)) {
    param <- testParam[i,]
    vgFitStartMoMReturn <- vgFitStartMoM (rvg(100, param = param))
    vgFitStartReturn <- vgFitStart(rvg(100, param = param))$paramStart
    vgFitNMReturn <- vgFit(rvg(100, param = param))
    vgFitBFGSReturn <- vgFit(rvg(100, param = param), method = "BFGS")
    vgFitnlmReturn <- vgFit(rvg(100, param = param), method = "nlm")
    
    checkTrue(is(vgFitStartMoMReturn, "numeric"), msg = paste(" param= ", param[1], 
      param[2], param[3], param[4]))
    checkEquals(length(vgFitStartMoMReturn ),4, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4])) 
    checkTrue(is(vgFitStartReturn, "numeric"), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(vgFitStartReturn ),4, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4])) 
    checkTrue(is.list(vgFitNMReturn), msg = paste(" param= ", param[1], param[2], param[3], 
      param[4]))
    checkEquals(length(vgFitNMReturn),14, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4])) 
    checkTrue(vgFitNMReturn$param[2] > 0, msg = paste(" param= ", param[1], param[2],
      param[3], param[4]))
    checkTrue(vgFitNMReturn$param[4] > 0, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkTrue(is.list(vgFitBFGSReturn), msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkEquals(length(vgFitBFGSReturn),14, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkTrue(vgFitBFGSReturn$param[2] > 0, msg = paste(" param= ", param[1], param[2],
      param[3], param[4]))
    checkTrue(vgFitBFGSReturn$param[4] > 0, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4])) 
    checkTrue(is.list(vgFitnlmReturn), msg = paste(" param= ", param[1], param[2], param[3], 
      param[4]))
    checkEquals(length(vgFitnlmReturn),14, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4])) 
    checkTrue(vgFitnlmReturn$param[2] > 0, msg = paste(" param= ", param[1], param[2], 
      param[3], param[4]))
    checkTrue(vgFitnlmReturn$param[4] > 0, msg = paste(" param= ", param[1], param[2],  
      param[3], param[4]))
  }
}




