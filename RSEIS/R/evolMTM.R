`evolMTM` <-
function(a, dt=0, numf=1024, Ns=0, Nov=0, fl=0, fh=10 )
  {
    ###  Nfft=1024;Ns=250;Nov=240;fl=0; fh=10
    if(missing(dt)) { dt=1;}
    if(missing(numf)) { numf=1024}
    if(missing(Ns)) { Ns=250;}
    if(missing(Nov)) { Nov=240;}
    if(missing(fl)) { fl=0;}
    if(missing(fh)) { fh=1/(2*dt);}
    

    Nfft = numf
    NT = length(a);
    nyquistf = 1/(2*dt)
    if(Nov<1)
      {
        Nov = floor(Ns - 0.1*Ns);
      }

     Ns = floor(Ns)
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

       krow = numf+1
    numfreqs=numf
    
 len2 = 2*next2(numf)

   ###  print(paste(sep=' ', "evolfft kcol=", kcol, "krow=", krow, "Ns", Ns, "Nov", Nov))
    if(kcol<1)
      {
        print(paste(sep=' ', "error in evolMTM kcol=", kcol, "krow=", krow, "NT", NT, "Ns", Ns, "Nov", Nov))
        return()
      }
          
    DMAT = matrix(rep(0,krow*kcol), ncol=kcol, nrow=krow)

    m = 1:(kcol)
    ibeg=((m-1)*skiplen)+1;
     iend = ibeg+Ns-1;
    
    for( i in m)
      { 
       ###  print(paste(sep=" ", m, ibeg, iend, NT))
        tem = a[ibeg[i]:iend[i]]
        tem = tem-mean(tem, na.rm=TRUE)

  Mspec =   mtapspec(tem,dt, klen=len2,  MTP=list(kind=2,nwin=5, npi=5,inorm=1)  )

  f=Mspec$freq

  amp = Mspec$spec[1:length(f)]

  DMAT[,i] = amp
      }
    
    DFFT = DMAT

   DSPEC = DMAT
     # col=heat.colors(50)


    
    x = (ibeg+Ns/2)*dt
    
    freqs = f

    y = f 
       
     ##   pdat = log(t(DMAT))
       
  ##  image(x, y, pdat, col=rainbow(100))
  ##  image(x, y, pdat, col=rainbow(100), ylim=c(0, 20) )


    ##     TEV = evolfft(a,dt, Nfft=4096, Ns=250 , Nov=240,  fl=0, fh=nyquistf  )
  ##       plotevol(TEV, log=1, fl=0, fh=nyquistf, col=rainbow(100))


       
    RET = list(sig=a, dt=dt, numfreqs=numfreqs, wpars=list(Nfft=numfreqs,  Ns=Ns, Nov=Nov, fl=fl, fh=fh), DSPEC=DSPEC, freqs=y, tims=x)

    ## plotevol(RET)
    
    invisible(RET)

  }

