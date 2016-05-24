DOsgram<-function(Xamp, DT = 0.008, multi = 1, scale.def =0, TWIN=2, TSKIP=.2, PCTTAP=0.05, PLOT=TRUE)
{
  if(missing(multi)) { multi= 1 }
  if(missing(scale.def)) { scale.def =0 }
  if(missing(TWIN)) { TWIN=2 }
  if(missing( TSKIP)) {TSKIP=.2}
  if(missing(PCTTAP)) {PCTTAP=0.05}
  if(missing(PLOT)) { PLOT = TRUE }
  STAMP = NULL
  
  gridon = FALSE
  TPALS = c("rainbow", "topo.colors", "terrain.colors", "heat.colors", "tomo.colors")
  pal = RPMG::Gcols(plow=5, phi=0,  N=100, pal=TPALS[1])


  fl=0
  fh=0.25*(1/DT)

  flshow = 0

###  fh is 60% the fh  frequency
  fhshow = round(0.6*fh)

###  here need to choose a length for the roving window.
###  I used to use 2 second windows by default, but here I think we need to
###  be more creative - how about 2% of the total length of the record?
  ##  tsecs = DT*(length(Xamp)*.02)
### No, now lets do 256 per wind approximately


  N = length(Xamp)
  NWIN = ceiling( TWIN/DT )
  Nskip = ceiling((TSKIP)/DT )

  NOV = NWIN-Nskip
  NS = NWIN

  Nhalf = floor(NS/2)

  Kcol  = seq(from=Nhalf, to=N-Nhalf, by=Nskip)

##  print(paste(sep=" ", "Len Kcol=", length(Kcol) ))
  
  if(FALSE)
    {
      tsecs = DT*256


      TWOSEC = tsecs*(1/DT)

      NS = floor(multi*TWOSEC)
      NOV = floor(multi*(TWOSEC-.2*TWOSEC))

    }

  Nfft=4096

  DEV = evolfft(Xamp, DT , Nfft=Nfft, Ns=NS , Nov=NOV,  fl=fl, fh=fh, pcttap=PCTTAP  )


  if(PLOT)  PE = plotevol(DEV, log=scale.def, fl=flshow, fh=fhshow, col=pal, ygrid=gridon, STAMP=STAMP)




  
  invisible(DEV)

}

