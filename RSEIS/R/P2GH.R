P2GH<-function(P1)
  {
    ##########  P1 = swig(GH)
    ########  select a small window on a trace, then click XTR
    
    DATTIM = P1$TIMEpick

    GG = list()
    N = length(P1$y)
    DATTIM$msec=0
    DATTIM$dt=P1$dt
    DATTIM$t1=0
    DATTIM$t2=N*P1$dt
    DATTIM$off=0

    GG[[1]] = list(fn="temp", sta=P1$STNS,  comp=P1$COMPS, dt=P1$dt, DATTIM=DATTIM, N=1, units="volts" , coords="temp" ,  amp=P1$y )

    LH = prepSEIS(GG)

    invisible(LH)

  }

