`MTMgabor` <-
function(a, dt=0, ppoint=95, numf=1024, Ns=0, Nov=0, fl=0, fh=10 )
  {
    ###  Nfft=1024;Ns=250;Nov=240;fl=0; fh=10
    if(missing(dt)) { dt=1;}
    if(missing(numf)) { numf=1024}
    if(missing(Ns)) { Ns=250;}
    if(missing(Nov)) { Nov=240;}
    if(missing(fl)) { fl=0;}
    if(missing(fh)) { fh=1/(2*dt);}
     if(missing(ppoint)) {  ppoint=95 }

    Ns = floor(Ns)
    Nfft = numf
    NT = length(a);
    nyquistf = 1/(2*dt)
    
    if(Nov<1)
      {
        Nov = floor(Ns - 0.1*Ns);
      }
    
    
    kcol =floor( (NT-floor(Nov) )/(Ns-floor(Nov)))
    if(kcol<Ns)
      {
        Ns = kcol
        Nov = floor(Ns-0.1*Ns)
         kcol =floor( (NT-floor(Nov) )/(Ns-floor(Nov)))
      }
    
    min1 = Nfft%%2;
    if(min1 == 0)
      {
        ## /* even */
        krow = (Nfft/2);
      } else {
        ##  /*  odd */
        krow = (Nfft+1)/2;
      }
    
    skiplen = Ns - Nov;
    
    df = 1.0/(Nfft*dt);

       krow = numf
    numfreqs=numf
    
 len2 = 2*next2(numf)

  #####  print(paste(sep=' ', "evolfft kcol=", kcol, "krow=", krow, "Ns=", Ns, "Nov=", Nov,  "numf=" , numf,  "len2=",   len2 ))
    if(kcol<1)
      {
        print(paste(sep=' ', "error in evolMTM kcol=", kcol, "krow=", krow, "NT", NT, "Ns", Ns, "Nov", Nov))
        return()
      }
          
    DMAT = matrix(rep(0,krow*kcol), ncol=kcol, nrow=krow)
    HIMAT = matrix(rep(0,krow*kcol), ncol=kcol, nrow=krow)
    DOFMAT = matrix(rep(0,krow*kcol), ncol=kcol, nrow=krow)
    FVMAT = matrix(rep(0,krow*kcol), ncol=kcol, nrow=krow)

    m = 1:(kcol)
    ibeg=((m-1)*skiplen)+1;
     iend = ibeg+Ns-1;
    ##   ppoints  =  c(50.0, 90.0, 95.0, 99.0, 99.5, 99.9)

###  print(paste(sep=" ", m, ibeg, iend, NT))


    ##########################################   check

    if(FALSE)
      {
    i = 1
        tem = a[ibeg[i]:iend[i]]
        tem = tem-mean(tem, na.rm=TRUE)
        
        ##  DOLEESPARK(tem,dt, tappercent=0.1)


        ######################
        Mspec =   mtapspec(tem,dt, klen=len2,  MTP=list(kind=2,nwin=5, npi=5,inorm=1)  )

        f=Mspec$freq
    print(paste(sep=" ", "len f=" , length(f)))

  }
##########################################   end check

    
    for( i in m)
      { 
###  print(paste(sep=" ", m, ibeg, iend, NT))
        tem = a[ibeg[i]:iend[i]]
        tem = tem-mean(tem, na.rm=TRUE)
        
        ##  DOLEESPARK(tem,dt, tappercent=0.1)


        
        Mspec =   mtapspec(tem,dt, klen=len2,  MTP=list(kind=2,nwin=5, npi=5,inorm=1)  )

        f=Mspec$freq

        ##############   mtapspec sends back one extra f value - discard this?
        findex  = 1:(length(f)-1)

        amp = Mspec$spec[findex]
        

        kdof = 2*Mspec$mtm$nwin-2
        
        myf = qf(ppoint/100, 2, kdof)
        fvals = Mspec$Fv[findex]
        w1 = which(fvals>myf)

########    plot(f, log(amp), type='l')
        HIMAT[w1  ,i] = 1
        FVMAT[,i] = fvals
        DOFMAT[,i] = Mspec$dof[findex]
        DMAT[,i] = amp
      }
    
   #### DFFT = DMAT

   DSPEC = DMAT
     # col=heat.colors(50)


    
    x = (ibeg+Ns/2)*dt
    
    freqs = f[findex]

    y = f[findex]
       
     ##   pdat = log(t(DMAT))
       
  ##  image(x, y, pdat, col=rainbow(100))
  ##  image(x, y, pdat, col=rainbow(100), ylim=c(0, 20) )


    ##     TEV = evolfft(a,dt, Nfft=4096, Ns=250 , Nov=240,  fl=0, fh=nyquistf  )
  ##       plotevol(TEV, log=1, fl=0, fh=nyquistf, col=rainbow(100))

  ####  print("FVMAT")
   ####   print(dim(FVMAT))
  ####    print("DOFMAT")
  ####    print(dim(DOFMAT))
  ####    print("HIMAT")
  ####    print(dim(HIMAT))
  ####    print("DSPEC")
  ####    print(dim(DSPEC))

       
    RET = list(sig=a, dt=dt, numfreqs=numfreqs, wpars=list(Nfft=numfreqs,  Ns=Ns, Nov=Nov, fl=fl, fh=fh), DSPEC=DSPEC, HIMAT=HIMAT, DOFMAT=DOFMAT, FVMAT=FVMAT, kdof=kdof,  ppoint=ppoint, freqs=y, tims=x)

    ## plotevol(RET)
    
    invisible(RET)

  }

