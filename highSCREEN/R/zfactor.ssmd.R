zfactor.ssmd = function(dat, pos.cont, neg.cont, MainPlate, triplicate){
  neg.cont.before = dat[[paste("scorebefore", triplicate, sep="")]][grep(neg.cont, as.character(dat[["welltype"]]))]
  pos.cont.before = dat[[paste("scorebefore", triplicate, sep="")]][grep(pos.cont, as.character(dat[["welltype"]]))]

  neg.cont.after = dat[[paste("scoreafter", triplicate, sep="")]][grep(neg.cont, as.character(dat[["welltype"]]))]
  pos.cont.after = dat[[paste("scoreafter", triplicate, sep="")]][grep(pos.cont, as.character(dat[["welltype"]]))]

  mean.neg.cont.before = mean(neg.cont.before)
  var.neg.cont.before = var(neg.cont.before)
  mean.pos.cont.before = mean(pos.cont.before)
  var.pos.cont.before = var(pos.cont.before)  
  cov.neg.cont.pos.cont.before = cov(neg.cont.before, pos.cont.before) 

  mean.neg.cont.after = mean(neg.cont.after)
  var.neg.cont.after = var(neg.cont.after)
  mean.pos.cont.after = mean(pos.cont.after)
  var.pos.cont.after = var(pos.cont.after)
  cov.neg.cont.pos.cont.after = cov(neg.cont.after, pos.cont.after) 
  
  SSMD_Before = (mean.neg.cont.before-mean.pos.cont.before) / sqrt(var.neg.cont.before+var.pos.cont.before-2*cov.neg.cont.pos.cont.before)
  SSMD_After = (mean.neg.cont.after-mean.pos.cont.after) / sqrt(var.neg.cont.after+var.pos.cont.after-2*cov.neg.cont.pos.cont.after)

  ZFactor_Before = 1-3*(sqrt(var.neg.cont.before)+sqrt(var.pos.cont.before))/abs(mean.pos.cont.before-mean.neg.cont.before)
  ZFactor_After = 1-3*(sqrt(var.neg.cont.after)+sqrt(var.pos.cont.after))/abs(mean.pos.cont.after-mean.neg.cont.after)

  res = data.frame(MainPlate, triplicate, ZFactor_Before, ZFactor_After, SSMD_Before, SSMD_After)
  return(res)
}
  
