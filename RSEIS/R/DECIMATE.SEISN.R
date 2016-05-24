`DECIMATE.SEISN` <-
function(TH, sel=1:length(TH$JSTR), dec=5 , type="LP", proto="C1" , fl=2,  fh=10 )
  {
    if(missing(sel)) { sel = 1:length(TH$JSTR) }
    if(missing(fh)) { fh=10 }
    if(missing(dec)) { dec=5 }
    
    for(i in 1:length(sel))
      {
        ii = sel[i]

        xamp = TH$JSTR[[ii]]
        mn = mean(xamp, na.rm=TRUE)
        nacop = is.na(xamp)
        xamp[nacop] = mn
        
        dt = TH$dt[ii]
        
        ynew = butfilt(xamp, fl=fl, fh=fh, deltat=dt, type=type, proto=proto )
        
        ynew[nacop] = NA
        TH$JSTR[[ii]] = ynew[seq(from=1, to=length(ynew), by=dec)]
        
        TH$dt[ii] = TH$dt[ii]*dec
        
        TH$info$dt[ii]  =TH$dt[ii]
        TH$info$n[ii] =TH$info$n1[ii] = TH$info$n2[ii] = length(TH$JSTR[[ii]])
        TH$info$t2[ii] = TH$info$dt[ii] * TH$info$n[ii]
      }
    return(TH)
  }

