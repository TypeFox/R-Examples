prepFOCS<-function(CMTSOL)
  {


    PHAIT = list(
      Paz =vector(),
      Pdip =vector(),
      Taz =vector(),
      Tdip =vector(),
      h = vector(),
      v = vector(),
      fcols = vector(),
      LATS =vector(),
      LONS = vector(),
      IFcol = vector(), yr=vector(), JDHM=vector() )

    for(i in 1:length(CMTSOL$str1) )
      {
        ##   M = FS2[[i]]$M

        MEC =  SDRfoc(CMTSOL$str1[i], CMTSOL$dip1[i], CMTSOL$rake1[i] , u=FALSE, ALIM=c(-1,-1, +1, +1), PLOT=FALSE)

        MEC$UP = FALSE
        PHAIT$IFcol[i] = foc.icolor(MEC$rake1)
        PHAIT$fcols[i] = foc.color(PHAIT$IFcol[i],1)

        az1 = MEC$az1
        dip1 = MEC$dip1
        az2 = MEC$az2
        dip2 = MEC$dip2
        BBB = Bfocvec(az1, dip1,  az2,  dip2)
        V = ternfoc.point(BBB$Bdip, MEC$P$dip, MEC$T$dip )
        ## print(BBB)
        ## print(V)

        PHAIT$Paz[i] = MEC$P$az
        PHAIT$Pdip[i] = MEC$P$dip
        PHAIT$Taz[i] = MEC$T$az
        PHAIT$Tdip[i] = MEC$T$dip

        PHAIT$h[i] = V$h
        PHAIT$v[i] = V$v
        PHAIT$LATS[i] = CMTSOL$lat[i]
        PHAIT$LONS[i] = CMTSOL$lon[i]
        PHAIT$yr[i] = CMTSOL$yr[i]
        PHAIT$JDHMS[i] = CMTSOL$jd[i]

      }

    
 invisible(PHAIT)

    
  }

