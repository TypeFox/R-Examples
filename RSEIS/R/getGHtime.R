getGHtime<-function(GH, wi=1, pix=NULL)
  {
    if(missing(wi)) { wi =1 } 
    if(missing(pix)) { pix=NULL }

    yr = GH$info$yr[wi]
    jd = GH$info$jd[wi]
    hr = GH$info$hr[wi]
    mi = GH$info$mi[wi]
    sec = GH$info$sec[wi]

    ORGtime = list(yr=yr, jd=jd, hr=hr, mi=mi, sec=sec)
    
    if(!is.null(pix))
      {
        spix=secdifL(ORGtime, pix)


      }
    else
      {

        spix= NULL
      }

    
    return(list(yr=yr, jd=jd, hr=hr, mi=mi, sec=sec, spix=spix))


  }
