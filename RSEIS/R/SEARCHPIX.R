`SEARCHPIX` <-
function(KPX, IPX, tol=0.5)
  {
    ####  IPX is the set of pix in memory
    ####  KPX is the user locator pix
    
    t1 =yeardate(IPX$yr, IPX$jd, IPX$hr,IPX$mi,  IPX$sec )
    t2 =yeardate(KPX$yr, KPX$jd, KPX$hr,KPX$mi,  KPX$sec )
    
    wn = which(  abs(t2-t1) < tol)
    
    return(wn)
    
  }

