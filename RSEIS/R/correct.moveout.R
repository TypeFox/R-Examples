correct.moveout<-function(GH, sel=1, tims=0)
  {
    ###  SHIFT traces by a given number of seconds
    ###  shift selected traces by a given number of seconds
    ###  shift by some moveout calculated previously
    #########   does moveout correction on set of traces in RSEIS
    #######   beware: this routine changes the data in the RSEIS structure
    
    if(missing(sel)) { sel = 1:length(GH$JSTR) }
    if(missing(tims)) { rep(0,length(GH$JSTR))  }


    
    for(i in 1:length(sel))
      {
        j = sel[i]
        if(tims[i]==0) next
        nsamps = round(abs(tims[i])/GH$dt[j])
        iz = sign(tims[i])
        a = GH$JSTR[[j]]
        n = length(a)
        if(iz<0)
          {
            b = c(a[(nsamps+1):n], rep(NA, length=nsamps))
          }
        else
          {
            b = c(rep(NA, length=nsamps),   a[1:(n-nsamps-1)] )

          }
        GH$JSTR[[j]] = b
        

      }

    invisible(GH)


  }
