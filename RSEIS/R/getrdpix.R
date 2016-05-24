getrdpix<-function(zloc,zenclick,sel,NH   )
  {
############   zloc are the clicks on the screen
############   zenclick is the number of clicks
############   sel is the selected index of traces
############   NH is the rseis structure
    
    kix = legitpix(sel, zloc, zenclick)

    ##  print(kix$ypick)
    
    if(length(kix$ypick)<1) {

      ## print("No legit picks")
      return(NULL)


    }
    
    ypick =  kix$ypick
    ppick =  kix$ppick
    
    ipick = sel[ypick]
    
    asec = NH$info$sec[ipick]+NH$info$msec[ipick]/1000+NH$info$t1[ipick]-NH$info$off[ipick]+ppick
    
    rd = recdate( NH$info$jd[ipick], NH$info$hr[ipick], NH$info$mi[ipick], asec, yr=NH$info$yr[ipick])
    
    
    rd$stn =  NH$STNS[ipick]
    rd$comp = NH$COMPS[ipick]
    invisible(rd) 
  }
