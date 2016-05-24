WINGH<-function(GH, sel=1 , WIN = c(0,1) )
{
###  GH is a standard RSEIS format list
  if(missing(sel)) sel = 1:length(GH$STNS)

### first copy the GH structure completely
  NEWH = GH

  
  NEWH$JSTR = vector(mode="list")
  NEWH$info = vector(mode="list")
  NEWH$dt = vector(mode="numeric")
  
  for (ii in 1:length(sel) ) {

    
    i = sel[ii]
    tim = GH$dt[i] * seq(from = 0, to = length(GH$JSTR[[i]]) -
      1)
    tflag = tim >= WIN[1] & tim <= WIN[2]
    amp = GH$JSTR[[i]][tflag]
    n1 = length(amp)
    NEWH$JSTR[[ii]] = amp
    NEWH$dt[ii] = GH$dt[i]
    RDATE = recdate(GH$info$jd[i], GH$info$hr[i], GH$info$mi[i],
      GH$info$sec[i] + GH$info$t1[i] + GH$info$msec[i]/1000 +
      WIN[1] - GH$info$off[i], GH$info$yr[i])
    GDOM = getmoday(RDATE$jd, RDATE$yr)
    NEWH$info$fn[ii] = GH$info$fn[i]
    NEWH$info$name[ii] = GH$info$fn[i]
    NEWH$info$yr[ii] = RDATE$yr
    NEWH$info$jd[ii] = RDATE$jd
    NEWH$info$mo[ii] = GDOM$mo
    NEWH$info$dom[ii] = GDOM$dom
    NEWH$info$hr[ii] = RDATE$hr
    NEWH$info$mi[ii] = RDATE$mi
    NEWH$info$sec[ii] = RDATE$sec
    NEWH$info$msec[ii] = 0
    NEWH$info$t1[ii] = 0
    NEWH$info$off[ii] = 0
    NEWH$info$t2[ii] = NEWH$info$t1[ii] + n1 * GH$dt[i]
    NEWH$info$n1[ii] = n1
    NEWH$info$n2[ii] = n1
    NEWH$info$n3[ii] = n1
    NEWH$info$n[ii] = n1
    
  }

NEWH$STNS = GH$STNS[sel]
NEWH$COMPS = GH$COMPS[sel]
NEWH$pcol  =   GH$pcol[sel]

NEWH$KNOTES  =  GH$KNOTES[sel]

NEWH$OCOMPS  =  GH$OCOMPS[sel]
NEWH$ftime  =  GH$ftime[sel]
  nn = length(NEWH$STNS)
   NEWH$nn = nn
  NEWH$ok  =  1:nn
  NEWH$ex =
    NEWH$dt[1]*seq(from=0,to=length(NEWH$JSTR[[1]])-1)

  return(NEWH)
}

