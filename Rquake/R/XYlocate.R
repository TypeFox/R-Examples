XYlocate<-function(Ldat,EQ,vel, 
                   maxITER = 10,
                   distwt = 10,
                   lambdareg =100,
                   FIXZ = FALSE,
                   REG = TRUE,
                   WTS = TRUE,
                   STOPPING = TRUE,
                   RESMAX = c(.4,.5),
                   tolx = 0.005,
                   toly = 0.005,
                   tolz = 0.01, PLOT=FALSE)
  {

###  the Ldat list must have members: x,y,err,sec, cor

####   the EQ list must have x,y,z,t
    
###  Earthquake location program for local earthquakes.
    ##  this code includes a projection to cartesian coordinates.
##########  EARTHQUAKE Location Program
    
    guesses = vector(mode='list')
    kguess = 0 
    NUMrows = length(Ldat$x)
    
     if(is.null(Ldat$cor)) { Ldat$cor =rep(0, NUMrows) }
      Ldat$err[Ldat$err<=0]  =  0.05


      kguess = kguess+1
            guesses[[kguess]] = list(x=EQ$x, y=EQ$y, z=EQ$z, t=EQ$t)
    

    for(K in 1:maxITER)
      {
#########   Now build up the partial derivative matrix:
        
        delx = EQ$x-Ldat$x
        dely = EQ$y-Ldat$y
        deltadis =sqrt( (delx)^2 +  (dely)^2)

        

        ROWZ = matrix(ncol=4, nrow = NUMrows)

        PredictedTT = rep(NA, length=NUMrows)
       
        wts = rep(1, length=NUMrows)

        if(WTS)
          {
            wts = DistWeightXY(Ldat$x, Ldat$y, EQ$x, EQ$y, Ldat$err, distwt)
            
          }
        
        G1 = GETpsTT(Ldat$phase, eqz=EQ$z, staz=0, delx=delx, dely=dely,  deltadis=deltadis , vel)

        kindex = Rowz2Keep(Ldat, EQ,  G1,  RESMAX)

        PredictedTT = EQ$t + G1$TT[kindex]
        Derivs = G1$Derivs[kindex, ]

        Observed = Ldat$sec[kindex]

        wheights = wts[kindex]

        cors = Ldat$cor[kindex]

        RHS  =  wheights  *  ( Observed - PredictedTT - cors )
        
        neqns = length(RHS)

        if(neqns<2)
          {

            cat(paste("############## BIG Problems: bad ROWZ"), sep="\n")
         #####     cat(paste(RESMAX), sep="\n")
            
         #####     cat(paste(kindex), sep="\n")
            
          #####  print(data.frame(Ldat))
           #####  testTT(Ldat,EQ, stas , vel)
          #####  print(ROWZ)
             return(list(EQ=NULL, its=NULL, rms=NULL, wrms=NULL))
 

          }

        
        ROWZ =  wheights  * cbind(rep(1, neqns), Derivs )

        

        ##############################
        
        if(FIXZ)
          {
            ROWZ = ROWZ[,1:3]
          }
        
        ##############  get the inverse

        if(any(!is.finite(ROWZ)))
          {
            cat(paste("############## BIG Problems: bad ROWZ"), sep="\n")
           ##### print(data.frame(Ldat))
           #####  testTT(Ldat,EQ, stas , vel)
           ##### print(ROWZ)
             return(list(EQ=NULL, its=NULL, rms=NULL, wrms=NULL))

          }

        
        S = svd(ROWZ)



        
        if(REG)
          {
            LAM = diag(S$d/(S$d^2+lambdareg^2) )
            h2 = S$v %*% LAM %*% t(S$u) %*% RHS
            
          }  else {
            h2 = S$v %*% diag(1/S$d) %*% t(S$u) %*% RHS
          }

        if(FIXZ)
          {
            h2 = c(h2, 0)
          }
       

###  update the solution:
        OLD = list(x=EQ$x, y=EQ$y)
        if(PLOT)
          {
            points(EQ$x, EQ$y, col='red' , cex=0.6, pch=8)
          }
        
        EQ$t = EQ$t+h2[1]
        EQ$x = EQ$x+h2[2]
        EQ$y = EQ$y+h2[3]
        EQ$z = EQ$z+h2[4]

        if(STOPPING)
          {
            if(abs(h2[2])<tolx & abs(h2[3])<toly & abs(h2[4])<tolz)
              {

                break
              }
          }

  ####
        kguess = kguess+1
        guesses[[kguess]] = list(x=EQ$x, y=EQ$y, z=EQ$z, t=EQ$t) 
        if(PLOT)
          {
            arrows(OLD$x, OLD$y, EQ$x, EQ$y, length=0.1, col='blue')
            
            points(EQ$x, EQ$y, col='red' , cex=0.6, pch=8)
            
          }
        
        
      }
#############  END LOOP
    

    wts = rep(1, length=NUMrows)

    if(WTS)
      {
        wts = DistWeightXY(Ldat$x, Ldat$y, EQ$x, EQ$y, Ldat$err, distwt)
            
      }
        

    delx = EQ$x-Ldat$x
    dely = EQ$y-Ldat$y
    deltadis =sqrt( (delx)^2 +  (dely)^2)
    
    G1 = GETpsTT(Ldat$phase, eqz=EQ$z, staz=0, delx=delx, dely=dely,  deltadis=deltadis , vel)
    
    kindex = Rowz2Keep(Ldat, EQ,  G1,  RESMAX)
    
    PredictedTT = EQ$t + G1$TT[kindex]
    
    
    Observed = Ldat$sec[kindex]
    
    wheights = wts[kindex]
    
    cors = Ldat$cor[kindex]

    resids = (Observed-(PredictedTT-cors))
    rms = sqrt(mean(resids^2))
    wrms = sqrt(mean((wheights*resids)^2))


    guesses= matrix(unlist(guesses), ncol=4, byrow=TRUE)
    
    return(list(EQ=EQ, its=K, rms=rms, wrms=wrms, used=kindex, guesses=guesses  ))


  }
