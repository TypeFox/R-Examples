symshot1<-function(PLOT=FALSE)
  {

    if(missing(PLOT)) { PLOIT = FALSE }

    
    nly=6
    ioff=3
    xdistance = 6.0
    iskip=2
    xspace= .1
    dt=.004
    maxtrc = 60
    velair = .350
    ampscal = 1.0
    ampnois = .2

    ampnois = 0.01
    ampreflex = 0.1
    amprefrac = 0.1
    ampair = 0.1

   rickfreq  = 25
   ricklen   = 35
    
    tlen = 2.4 

    tracelen=600
    wavelen=50


    GRLoc=list()
    GRLoc$x=c(6.01078415186,2.77768081931)
    GRLoc$y=c(-1.51478902016,-2.38634599907)

    GRLslope1 = GRLoc$y[1]/  GRLoc$x[1]
    GRLslope2 = GRLoc$y[2]/  GRLoc$x[2]
    

 ######    abline(0, GRLoc$y[1]/  GRLoc$x[1], col='blue')
 ######    abline(0, GRLoc$y[2]/  GRLoc$x[2], col='blue')
   
    
    
##########     m is the velocity model

    m=matrix(c(0.0, 3.0,
      0.3, 4.0,
      1.0, 6.5,
      2.0, 7.0,  
      4.5, 8.5,
      7.5, 9.5,
      50.0, 10), ncol=2, byrow=TRUE)


    y = m[1:6,1]
    v = m[1:6 ,2]

    x1 = rep(0, length(y))
    x2 = rep(1, length(y))


    x  = seq(from=ioff*xspace, to=xdistance, by=xspace)
    
    t2tim = seq(from=0, to=tlen, by=dt)
    
    NLEN = length(t2tim)
    Ntrace = length(x)

    tair = x/velair

    klay = 3

    h = diff(y)
    xcrit = vector()

    trefrac = matrix(ncol=length(x), nrow=klay)
    treflex = matrix(ncol=length(x), nrow=3)

    xcrit = vector(length=nly)
    top = y
    vel = v
    xc = 0
    for( i in 1:(nly-1) )
      {
        
        xc = xc +2*h[i]* sqrt(   (v[i+1]/v[i])^2 - 1 )
        xcrit[i] = xc
      }
############  Calculate the theoretical arrival times based on geometry
############
############   refractions

    for( n in 1:klay)
      {
        vn = v[n]
        i = seq(from=1, to=n, by=1)
        thetai = asin(v[i]/vn)
        tn = sum(2*h[i]*cos(thetai)/v[i])
        trefrac[n,] = x/v[n] + tn 
        if(n>1) trefrac[n, x<xcrit[n-1] ] =NA
      }

############  calculate reflections

    k = 0
    for(n in 3:(nly-1))
      {
        k = k+1
        Delt0 = h[1:n]/v[1:length(h[1:n])]
        vrms = sqrt( sum( v[1:n]^2*Delt0)/sum(Delt0) )
        t0 = sum(Delt0[1:n])
        tt = sqrt((x^2)/vrms^2 + t0^2)
        ## print(c(n, vrms, t0))
        ## print(range(tt))
        ##points(x, tt, col=n)
        treflex[k,] = tt
      }

############
    
###############   set up wavelets:

 ############   ricker wavelet, shifted so it will be centered on the spike

    wavelet = genrick(rickfreq,dt,ricklen)
    klem = 11
### nwave = -1 + 2 * (wavelet - min(wavelet) )/diff(range(wavelet))
    nwave = RPMG::RESCALE(wavelet, 0, 1, wavelet[1], max(wavelet))




    shkip = length(nwave)/2
    ones = rep(1, klem)
    jspred = applytaper(ones, p = 0.5)
    ##  plot(jspred)
    nspred = 0.5+length(jspred)/2

    reach = seq(from =1, by=shkip, length=klem)


###############   ground roll:
#################  this is a sinusoidal signal
    
    grlen = floor(.6/dt)
    fgr = 10
    tape = applytaper( rep(1, grlen), p = 0.2)
    tgr = seq(from=0, by=dt, length=grlen)
    siggr = tape*sin(2*pi*fgr*tgr)
   ##  plot(tgr, siggr, type='l')

#######################################################################
#######################################################################

    x1 = ampnois*runif(NLEN, -1, 1)

    KL = length(nwave)

###
    smograms = matrix(ncol=Ntrace , nrow=NLEN)

    #################################################   Air wave
    AIRrec = matrix(ncol=Ntrace , nrow=NLEN)

    for(i in 1:Ntrace)
      {

        x1 = rep(0, times=(NLEN))
        ##  x1 = ampnois*runif(NLEN, -1, 1)

        ##  air
        iair = round( tair[i]/dt )

        if(iair>0 & iair<NLEN)
          {
            x1[iair] = x1[iair]+ampair

          }
        cx1 = x1

        AIRrec[,i ] = cx1

      }

    AIRrec = sigconv(AIRrec, nwave)

#########################################################   ground roll
    n = 1
    GRrec = matrix(ncol=Ntrace , nrow=NLEN)

    for(i in 1:Ntrace)
      {
        x1 = rep(0, times=(NLEN))
        zim  = round( trefrac[n,i]/dt )

        if(zim>0 & zim<NLEN & !is.na(zim))
          {
            x1[zim] = x1[zim]+amprefrac
            
          }
        cx1 = x1

        grlen = floor(.6/dt)
        fgr = 10
        tape = applytaper( rep(1, grlen), p = 0.2)
        tgr = seq(from=0, by=dt, length=grlen)
        siggr = tape*sin(2*pi*fgr*tgr)

        GRrec[,i ] = cx1


      }

    GRrec = sigconv(GRrec, siggr)


    ############################################   refractions
    REFR = matrix(ncol=Ntrace , nrow=NLEN)

    for(i in 1:Ntrace)
      {
        x1 = rep(0, times=(NLEN))


        for(n in 2:klay)
          {
            zim  = round( trefrac[n,i]/dt )

            if(zim>0 & zim<NLEN & !is.na(zim))
              {
                x1[zim] = x1[zim]+amprefrac
                
              }
          }


        cx1 = x1

        REFR[,i ] = cx1


      }

    REFR  = sigconv(REFR , nwave)


##################################   reflections

    REFL = matrix(ncol=Ntrace , nrow=NLEN)

    for(i in 1:Ntrace)
      {
        x1 = rep(0, times=(NLEN))


        for(n in 1:klay)
          {
            zim  = round( treflex[n,i]/dt )

            if(zim>0 & zim<NLEN & !is.na(zim))
              {
                x1[zim] = x1[zim]+ampreflex
                
              }
          }


        cx1 = x1

        REFL[,i ] = cx1


      }

    REFL = sigconv(REFL , nwave)

    smograms = REFL + REFR +GRrec +AIRrec

    if(PLOT==TRUE) wiggleimage(smograms  , dt=(-dt), dx=x)

    THEORY = list(trefrac=trefrac, treflex=treflex, tair=tair, velair=velair, mod=m)
    dx=x[2]-x[1]

      invisible(list( smograms = smograms,   dt=dt,  x=x,  dx=xspace , REFL=REFL,  REFR=REFR,  GRrec=GRrec  , AIRrec=AIRrec ,   THEORY = THEORY))  

    }
