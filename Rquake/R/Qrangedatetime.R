Qrangedatetime<-function(D)
{
  ##   given a list of dates, get the range and difference

  if(is.list(D))
    {

      D$yr[is.na(D$yr)] = D$yr[!is.na(D$yr)][1]
      D$jd[is.na(D$jd)] = D$jd[!is.na(D$jd)][1]
      D$hr[is.na(D$hr)] = D$hr[!is.na(D$hr)][1]
      D$mi[is.na(D$mi)] = D$mi[!is.na(D$mi)][1]
      D$sec[is.na(D$sec)] = D$sec[!is.na(D$sec)][1]

      
      NN = lapply(D, "length")
      
      n = max(c(NN$yr, NN$jd, NN$hr, NN$mi, NN$sec), na.rm=TRUE)

      RDtmes = D$yr+D$jd/366+D$hr/(366*24)+D$mi/(366*24*60)+D$sec/(366*24*3600)

      w1 = which.min(RDtmes)
      w2 = which.max(RDtmes)

      if(NN$yr<n) { D$yr = rep(D$yr[1],times=n) }
      if(NN$jd<n) { D$jd = rep(D$jd[1],times=n) }
      if(NN$hr<n) { D$hr = rep(D$hr[1],times=n) }
      if(NN$mi<n) { D$mi = rep(D$mi[1],times=n) }
      if(NN$sec<n) { D$sec = rep(D$sec[1],times=n) }

      T1 = list(jd=D$jd[w1], hr=D$hr[w1], mi=D$mi[w1], sec=D$sec[w1], yr=D$yr[w1])
      T2 = list(jd=D$jd[w2], hr=D$hr[w2], mi=D$mi[w2], sec=D$sec[w2],  yr=D$yr[w2])

      ydif = RSEIS::YRsecdif(T1$jd, T1$hr, T1$mi, T1$sec,  T2$jd, T2$hr, T2$mi, T2$sec,
        D$yr[w1] ,  D$yr[w2])
      return(list(min=T1, max=T2, dif=ydif))
    }
  else
    {
      return(NULL)
    }
}

