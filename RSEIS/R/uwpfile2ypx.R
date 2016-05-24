`uwpfile2ypx` <-
function(P)
{
  P = cleanpickfile(P)
  YPX=P$STAS
  n = length(P$STAS$name)
  YPX$yr = rep(P$LOC$yr, length=n)
  YPX$mo = rep(P$LOC$mo, length=n)
  YPX$dom = rep(P$LOC$dom, length=n)
  YPX$jd = rep(P$LOC$jd, length=n)
  YPX$hr = rep(P$LOC$hr, length=n)
  YPX$mi = rep(P$LOC$mi, length=n)
  YPX$sec = P$STAS$sec
  
  pcol=rep("springgreen4", length(YPX$sec))
  phas = YPX$phase
  pcol[phas=="P"] = "violetred2"
  pcol[phas=="S"] = "deepskyblue4"
  
  YPX$col = pcol
  YPX$onoff = rep(0, length(YPX$sec))

  return(YPX)

}

