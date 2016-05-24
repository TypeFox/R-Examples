GETpsTT<-function(phase, eqz=6, staz=0, delx=1, dely=1,  deltadis=6 , vel)
  {
    ########   get travel times and derivatives for multiple eq-stations
    neqns = length(deltadis)

    if(length(staz)<neqns) { staz=rep(staz, neqns) }
    
    ROWZ = matrix(rep(0, neqns*3), ncol=3, nrow = neqns)

    PredictedTT = rep(NA, times=neqns)

    Pphase = phase=="P"
    Sphase = phase=="S"

    if(any(Pphase))
      {
        wp = which(Pphase)

        udist = deltadis[wp]
        np = length(udist)
        
        Hp = RSEIS::many.time1D(udist , eqz , staz[wp], length(vel$zp)  , vel$zp ,  vel$vp)
        
         dtdxp = rep(0, np)
        dtdyp =  dtdxp

        
        useP = wp[udist>0]
        dtdxp[udist>0] = Hp$dtdr[udist>0] *delx[useP ]/udist[udist>0]
        dtdyp[udist>0] = Hp$dtdr[udist>0] *dely[useP ]/udist[udist>0]
        
       
        dtdzp = Hp$dtdz
        
        ROWZ[wp, ] = cbind(dtdxp, dtdyp, dtdzp)
        PredictedTT[wp] = Hp$tt

      }


    if(any(Sphase))
      {
        ws = which(Sphase)

        udist = deltadis[ws]
        ns = length(udist)
        
        Hs = RSEIS::many.time1D(udist , eqz , staz[ws], length(vel$zs)  , vel$zs ,  vel$vs)
         dtdxs = rep(0, ns)
        dtdys =  dtdxs

        useS =  ws[udist>0]

        dtdxs[udist>0] = Hs$dtdr[udist>0] *delx[useS ]/udist[udist>0]
        dtdys[udist>0] = Hs$dtdr[udist>0] *dely[useS ]/udist[udist>0]
        
        dtdzs = Hs$dtdz
        
        ROWZ[ws, ] = cbind(dtdxs, dtdys, dtdzs)
        PredictedTT[ws] = Hs$tt

      }
    
    return(list(TT=PredictedTT, Derivs=ROWZ))

  }


