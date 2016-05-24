xtract.trace<-function(GH, sel=1, WIN=c(0,1) )
{
  if(missing(sel)) sel = 1

  y1= GH$JSTR[[sel]]
  ex = seq(from=0, length=length(y1), by=GH$dt[sel] )


  if(missing(WIN)) {

    WIN = range(ex)

  }
  y  = y1[ex>=WIN[1] & ex<=WIN[2] ]
  
GH$info$sec[sel] = GH$info$sec[sel]+WIN[1] 
  ZD = Zdate(GH$info)

  
  attr(y, "dt")<- GH$dt[sel]
  attr(y, "start")<- ZD[sel]
  attr(y, "name")<- GH$STNS[sel]
  attr(y, "comp")<- GH$COMPS[sel]

  
  return(y)

}
