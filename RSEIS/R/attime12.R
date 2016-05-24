attime12<-function(t1, t2=t1, origyr=1972, pre=0, post=0)
  {
    ####   get two times suitable for use in extracting
    ####   data from a DB list and file system
    if(missing(t2)) t2 = t1
    if(missing(origyr)) origyr=1972
    if(missing(pre)) pre=0
    if(missing(post)) post=0
    
    

    eday1 = EPOCHday(t1$yr[1] , jd = t1$jd[1], origyr = origyr[1])
    eday2 = EPOCHday(t2$yr[1] , jd = t2$jd[1], origyr = origyr[1])

    at1 = eday1$jday +t1$hr[1]/24+ t1$mi/(24*60)+(t1$sec-pre)/(24*3600)
    at2 = eday2$jday +t2$hr[1]/24+ t2$mi/(24*60)+(t2$sec+post)/(24*3600)


    return(c(at1, at2)) 


  }


