`T12.pix` <-
function(A)
  {
    if(is.null(A$dur)) A$dur=A$res

   origyr = attr(A, "origyr")
    if(is.null(origyr)) { origyr = min(A$yr, na.rm = TRUE)  }
    
    eday = EPOCHday(A$yr,A$jd, origyr = origyr)
    
    A$t1 = eday$jday  + A$hr/24 + A$mi/(1440) + A$sec/(86400) 
    A$t2 = A$t1 + A$dur/86400


    attr(A, "origyr")<-origyr 
    
    invisible(A)
  }

