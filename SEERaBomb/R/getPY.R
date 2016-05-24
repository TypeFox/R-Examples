getPY=function(surv,bin,binS,brks) {
  binLengths=diff(brks)
  LL=getBinInfo(bin,binS)["LL"]
  UL=getBinInfo(bin,binS)["UL"]
  indx=getBinInfo(bin,binS)["index"]
  ifelse(surv>=UL,binLengths[indx],ifelse(surv>LL,surv-LL,0))
}
