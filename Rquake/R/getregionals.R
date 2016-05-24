getregionals <-
function(KAT, Mlat, Mlon, rad=1000 ,  t1=1, t2=2)
  {
###  given a katalog a lat and lon
    ###   radius and a time window, extract the events that comply
    ###   get the katalog from function prepPDE
    d = GEOmap::distaz(Mlat, Mlon, KAT$lat, KAT$lon)

    if(!is.null(rad)  )
      {
        wrad = (d$dist<rad)
      }
    else
      {
        wrad = rep(TRUE, length(KAT$jsec))

      }

   if(!is.null(t1) & !is.null(t2) )
     { wt =  KAT$jsec>t1 & KAT$jsec<t2 }
    else
      {
        wt = rep(TRUE, length(KAT$jsec))
      }

    wsel = which(wrad & wt)
    return(wsel)
  }
