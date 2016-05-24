Klocate<-function(Ldat,
                  sol=c(0,0,0,0),
                  vel=defaultVEL(6),
                  distwt=20,
                  errtol=c(0.01, 0.01, 0.01),
                  maxit=20 ,
                  Lambda=1,
                  guessdepth = 6,
                  APLOT = FALSE ,
                  stas=list(name="", lat=NA, lon=NA, z=NA)  )
  {

   
    
###   this location program has removed the projections.
###   distances are calculated between station and hypocenter
    ### and used in the travel time calculation

    
    ##   errtol = vector, 3 (ex,ey,ez) tolerances for stopping iterations

###  maxit = maximum number of non-linear iterations
    
    ##  distwt =  a single constant used to weight
    ##         observations based on their distance from the event.
    
    ##   vel = defaultVEL(kind = 6)
###  find first arrival

 dtor = pi/180
   rtod = 180/pi  
     rearth = 6371.0;

    
    BPLOT = FALSE
    if(missing(APLOT) ) { APLOT = FALSE }
    
    if(missing(distwt)) distwt=20
    if(missing(errtol)) errtol=c(0.01, 0.01, 0.01)
    if(missing(vel)) { vel=defaultVEL(6) }

    if(missing(Lambda))
      {
        lambdareg = 1
      }
    else
      {
        lambdareg = Lambda
      }

    if(missing(maxit)) maxit=20
    
    if(missing(guessdepth) )  guessdepth = 6
    
  
    if(missing(stas)) { stas=list(name="", lat=NA, lon=NA, z=NA) }
    
###     A1 = Gfirstguess(Ldat, type="median")
    
###  establish coordinate system, convert to XY
###    proj = GEOmap::setPROJ(type = 2, LAT0 = A1[1], LON0 = A1[2])

###    XY = GEOmap::GLOB.XY(Ldat$lat, Ldat$lon, proj)
###    EQ = GEOmap::GLOB.XY(sol[1], sol[2] , proj)


    ##########  must have station lat-lon-z values 
    if(is.null(Ldat$lat) | is.null(Ldat$lon))
      {

        if(!is.na(stas$lat[1]))
          {
            Ldat = latlonz2wpx(Ldat, stas)
      }
        else
          {
            print("ERROR: No LAT-LON-Z for stations")
            return(NULL)

          }

      }
######   must have an error estimate
    if(any(Ldat$err<=0))
      {
        Ldat$err[Ldat$err<=0] = 0.1

      }

    ######  if no P or S phases, assume all are first arrivals (P)
    if(  !any(Ldat$phase=="P") &  !any(Ldat$phase=="S") )
      {
        Ldat$phase=rep("P", length(Ldat$phase))

      }

  if(missing(sol) )
      {
        A1 = Gfirstguess(Ldat, type="first")
        sol = A1
        sol[3] = guessdepth
        sol[4] = 0 
      }
    if(is.null(sol) )
      {
        A1 = Gfirstguess(Ldat, type="first")
        sol = A1
        sol[3] = guessdepth
        sol[4] = 0 
      }


 ###    print(data.frame(Ldat))
    
    EQ = list(lat=sol[1], lon=sol[2], z=sol[3] ,t=0)
    
    if(APLOT)
      {
        plot(  RPMG::fmod(Ldat$lon, 360) , Ldat$lat, pch=6)
        points(RPMG::fmod( EQ$lon, 360) , EQ$lat, pch=8, col='red')
      }
    
    for(Kiters in 1:maxit)
      {
        
        
###  Set up Equations
        neqns = length(Ldat$sec)
        
        ROWZ = matrix(ncol=4, nrow = neqns)
        
        RHS = rep(0, length=neqns)

       DEL =  GEOmap::distaz(EQ$lat, EQ$lon, Ldat$lat, Ldat$lon)

        deltadis = DEL$dist
        daz =   DEL$az*dtor

        cosAZ = cos(daz)
        sinAZ = sin(daz)
        cosLAT = cos(Ldat$lat * dtor)
        
############   set up weights for each equation
        wts = rep(1, neqns)

        
        temp = Ldat$err
        dwt =  (1.0/(1. + ((deltadis^2)/(distwt^2))))
#########   this is the Lquake weighting scheme:
        wts = dwt/sqrt(temp^3);
        
        ##  W = diag(wts)

        ###################    go through the data and accumulate the matrix

        for(j in 1:neqns)
          {

            if(Ldat$phase[j]!="P" & Ldat$phase[j]!="S") next


           if(Ldat$phase[j]=="P")
              {
                TAV =   RSEIS::travel.time1D(deltadis[j], EQ$z , 0, length(vel$zp) , vel$zp , vel$vp)
              }
            if(Ldat$phase[j]=="S")
              {
                TAV =   RSEIS::travel.time1D(deltadis[j], EQ$z , 0, length(vel$zs) , vel$zs , vel$vs)
              }

            dtdz = TAV$dtdz

            if(is.nan(dtdz)) { dtdz =  0 }
            

            if(deltadis[j] == 0 )
              { 
                dtdx = 0
                dtdy = 0
              }
            else
              {
            	dtdx  = -TAV$dtdr * cosAZ[j] * rearth * dtor;
		dtdy =  -TAV$dtdr * sinAZ[j] * cosLAT[j] * rearth * dtor;
                
              }

##########   apply weights and accumulate 
            ROWZ[j, ] = wts[j] * c(1, dtdx, dtdy,   dtdz) 
            RHS[j]  =   wts[j] *  ( Ldat$sec[j] - EQ$t - TAV$tt ) 
          }

      #####   print(ROWZ)
     ###   print(RHS)

########  establish regularization....
        
###########  SVD solution
        S1 = svd(ROWZ)

     ###   print(S1$d)
        
        LAM = diag(S1$d/(S1$d^2+lambdareg^2) )

        Gdagger = S1$v %*% LAM %*% t(S1$u)
        Ssol = Gdagger  %*% RHS

     ###    print( Ssol )   
##########   get perturbation solution

####  dLL = GEOmap::XY.GLOB(EQ$x+Ssol[1], EQ$y+Ssol[2], proj   )

####  A1 = c(dLL$lat, dLL$lon, sol[3]+Ssol[3], sol[4]+Ssol[3])

        newlocX = EQ$lat+Ssol[2,]
        newlocY = EQ$lon +Ssol[3,]
        newlocZ = EQ$z + Ssol[4,]
        newlocT = EQ$t +Ssol[1,]

        if(APLOT)
          {
        points( newlocX ,   newlocY, col='purple' , pch=3) 
        arrows(EQ$x,EQ$y,    newlocX ,  newlocY, col='green', length=0.02 )
      }

        EQ$lat = newlocX
        EQ$lon = newlocY
        EQ$z = newlocZ
        EQ$t = newlocT
        if(abs(Ssol[4])<errtol[3] & abs(Ssol[2])<errtol[1] & abs(Ssol[3])<errtol[2] )
          {
            ##  print(paste("Iterations=",Kiters))

           
            break
          }

      }

    ###  calculate the error bars

  ##   covB = t(S1$v) %*% S1$v

    ####  covariance of the data:
    covD = diag(Ldat$err^2)

    ###  covariance of the model parameters
    covB =   Gdagger %*% covD  %*% t(Gdagger)

    #######  extract the diagonals (variances) and get sqrt
    dels = sqrt(diag(covB))

    return(list(lat=EQ$lat, lon=EQ$lon, z=EQ$z, t=EQ$t, QUAL=0, ITS=Kiters, cov=covB ))
  }
