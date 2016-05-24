pickseis24<-function(pjj, DB, usta, ucomp )
{
  w = list(hr=c(1,1) )

  
origyr=attr(DB, "origyr")
  
  while(length(w$hr)>1)
    {
      dev.set(2)
      w = winseis24(pjj)
      print(w)

      eday = EPOCHday(w$yr, jd =w$jd, origyr = origyr)


      at1 = eday$jday + (w$hr[1]) /24 
      at2 = eday$jday + (w$hr[2])/24

      GH = Mine.seis(at1, at2, DB, usta, ucomp)

      dev.set(3)
      swig(GH)

      dev.set(2)



    }

}

