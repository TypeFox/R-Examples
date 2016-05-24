save.wpix<-function(KOUT, fn= 'wpix.out')
  {
    if(missing(fn)) { fn= 'wpix.out' }


        fac = file.access(fn, mode=2)
        if(fac!=0)    fac =  file.create(fn)
    ####      if(fac!=0)
    ####        {
    ####          print("failed to save wpix.  Permission to write file denied")
    ####          return(0)

    ####        }

    notisna = !is.na(KOUT$WPX$yr) & KOUT$WPX$onoff>0

    
    d = data.frame(yr=KOUT$WPX$yr[notisna],
      mo=KOUT$WPX$mo[notisna],
      dom=KOUT$WPX$dom[notisna],
      jd=KOUT$WPX$jd[notisna],
      hr=KOUT$WPX$hr[notisna],
      mi=KOUT$WPX$mi[notisna],
      sec=KOUT$WPX$sec[notisna],
      dur=KOUT$WPX$res[notisna])
    
    write.table(file=fn, d, row.names = FALSE, col.names = FALSE, append=TRUE)


  }


