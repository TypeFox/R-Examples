CHOP.SEISN = function (GH, sel = 1:4, WIN = NULL)
{
    if (missing(sel)) {
        sel = 1:length(GH$dt)
    }
    if (is.logical(sel)) {
        sel = which(sel)
    }
    NEWH = GH
    if (missing(WIN)) {
        WIN = NULL
        zloc = ZOOM.SEISN(GH, sel, WIN = WIN)
    }
    else {
        zloc = WIN
    }
    if (is.list(zloc) == FALSE) {
        zloc = list(x = zloc)
    }
    NEWH$JSTR = list()
    NEWH$dt = NULL
    NEWH$info$yr = NULL
    NEWH$info$jd = NULL
    NEWH$info$mo = NULL
    NEWH$info$dom = NULL
    NEWH$info$hr = NULL
    NEWH$info$mi = NULL
    NEWH$info$sec = NULL
    NEWH$info$msec = NULL
    NEWH$info$t1 = NULL
    NEWH$info$off = NULL
    NEWH$info$t2 = NULL
    NEWH$info$n1 = NULL
    NEWH$info$n2 = NULL
    NEWH$info$n3 = NULL
    NEWH$info$n = NULL
    NEWH$info$gain = NULL
    NEWH$info$scalefac = NULL

    NEWH$STNS = GH$STNS[sel]
    NEWH$COMPS = GH$COMPS[sel]
    NEWH$OCOMPS = GH$OCOMPS[sel]
    NEWH$KNOTES = GH$KNOTES[sel]
    NEWH$pcol = GH$pcol[sel]
    NEWH$ok = GH$ok[sel]
    NEWH$ftime = GH$ftime[sel]
   
    for (i in sel) {
        ii = which(i==sel)
        tim = GH$dt[i] * seq(from = 0, to = length(GH$JSTR[[i]]) -
            1)
        tflag = tim >= zloc$x[1] & tim <= zloc$x[2]
        amp = GH$JSTR[[i]][tflag]
        n1 = length(amp)
        NEWH$JSTR[[ii]] = amp
        NEWH$dt[ii] = GH$dt[i]
        RDATE = recdate(GH$info$jd[i], GH$info$hr[i], GH$info$mi[i],
            GH$info$sec[i] + GH$info$t1[i] + GH$info$msec[i]/1000 +
                zloc$x[1] - GH$info$off[i], GH$info$yr[i])
        GDOM = getmoday(RDATE$jd, RDATE$yr)
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
       
 NEWH$info$gain[ii] =  GH$info$gain[i]

 	  NEWH$info$scalefac[ii] = GH$info$scalefac[i]
          
       
    }
     NEWH$wintim = NEWH$info$jd[1] + NEWH$info$hr[1]/24 + NEWH$info$mi[1]/(24 * 60) +
          (NEWH$info$sec[1] + NEWH$info$msec[1]/1000 + NEWH$info$t1[1] - NEWH$info$off[1])/(24 * 3600)

        
    NEWH$ex = NEWH$dt[1] * seq(from = 0, to = length(NEWH$JSTR[[1]]) -
        1)
    invisible(NEWH)
}
