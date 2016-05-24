`longfft` <-
function(DB, DAYS=c(233,234), HRS=1:24 , sta="KR1", comp=c("V", "I"), NPP=6 , CSCALE=FALSE, pal=rainbow(100), PS=FALSE , kind = 1,  Iendian=1, BIGLONG=FALSE )
{
  ###  longfft  =     Long Evolutive FFT
     if (missing(kind)) {
      kind = 1
    }
      if(missing(Iendian)) { Iendian=1 }
     if(missing(BIGLONG)) { BIGLONG=FALSE}

 if(missing(CSCALE)) { CSCALE=TRUE  }
 if(missing(NPP)) { NPP=6  }  ###########  number of plots per page
 if(missing(PS)) { PS=FALSE  }
 if(missing(HRS)) { HRS=1:24 }
 if(missing(pal)) {

   pal = list(p1=RPMG::Gcols(plow=5, phi=0,  N=100, pal="topo.colors", mingray=0.8),
     p2=RPMG::Gcols(plow=5, phi=0,  N=100, pal="rainbow", mingray=0.8))
 }
     
### if(missing(Use)) { Use = seq(from=1, to=length(gnam), by=2) }
### dev.set(2)

 if(is.list(pal))
   {
     pal1 = pal[[1]]
     pal2 = pal[[2]]
   }
 else
   {
     pal1 = pal
     pal2 = pal

   }
def.par = par(no.readonly = TRUE) 
  
 
 if(PS==FALSE)
   {
   ### graphics.off()
  
    par(mfrow=c(NPP,1))
    par(mai=c(0.3,0.35,0.3,0.1))
   }
 

  fl=0
  fh=15
  fhshow =  15
  flshow = 0   
  
  for(i in 1:length(DAYS) )
    {

      theday = DAYS[i]

      KK = 1
      kplot = 0

      for(j in 1:length(HRS))
        {

          print(paste(sep=' ',"##########################################",KK, (KK %% (NPP/2))) )

          
          hr = HRS[j]
          
         
          
          
          at1 = theday+hr/24
          at2 = at1+1/24

          print(paste(sep=' ',"########", theday, hr, at1, at2, sta, comp))
          KH = Mine.seis(at1, at2, DB, sta, comp, kind = kind, Iendian=Iendian, BIGLONG=BIGLONG )

          if(is.null(KH)) { next }


           if( identical(kplot, 0) )
            {
              if(PS==TRUE)
                {
                  kplot = longpstart(NPP=NPP, asta=sta, acomp=comp, theday=theday, hr=hr)
                }

            }
          
          
          ##  swig(KH)

          if(length(KH$JSTR)>0)
            {
              if(length(KH$JSTR)==2)
                {
                  if(kplot>NPP-2)
                    {
                      kplot=longreset(NPP, PS)
                      kplot = longpstart(NPP=NPP, asta=sta[1], acomp="ALL", theday=theday, hr=hr)
                    }
                }
              ipick = 1
              famp = KH$JSTR[[ipick]]
              DT=KH$dt[ipick]
###  plot(famp)
              
              temp =  famp 
              Xamp =   temp-mean(temp, na.rm =TRUE)

####   SPECT.drive(Xamp, DT=KH$dt[ipick], NEW=FALSE)
          NS = round((0.1*60)/DT)

          NOV = round(NS/10)
          NFFT = next2(NS*2)

           print(paste(sep=' ',"########", length(Xamp), NS, NOV, NFFT, fl, fh))

          DEV = evolfft(Xamp, DT , Nfft=NFFT, Ns=NS , Nov=NOV,  fl=fl, fh=fh  )

          plotevol(DEV, log=0, fl=flshow, fh=fhshow, col=pal1, ygrid=FALSE, AXE=c(1,2) , CSCALE=CSCALE)
         ## image(t(DEV$DSPEC[1:(DEV$numfreqs/2),], col=pal1))
          anam = KH$info$fn[ipick]

          snam = unlist( strsplit(split="/", anam))

          inam = snam[length(snam)]
           mtext(paste(sep=' ',"Day=",theday, "Hour=", hr, KH$STNS[ipick], KH$COMPS[ipick] ), side=3, at=0, adj=0, col=4)
          kplot = kplot +1
          
        }

          if(length(KH$JSTR)>1)
            {
          ipick = 2
          famp = KH$JSTR[[ipick]]
          DT=KH$dt[ipick]
###  plot(famp)

          temp =  famp 
          Xamp =   temp-mean(temp, na.rm=TRUE)

####   SPECT.drive(Xamp, DT=KH$dt[ipick], NEW=FALSE)
          NS = round((0.1*60)/DT)

          NOV = round(NS/10)
          NFFT = NS*4

          DEV = evolfft(Xamp, DT , Nfft=NFFT, Ns=NS , Nov=NOV,  fl=fl, fh=fh  )
          ##  dev.set(dev.next())

          plotevol(DEV, log=0, fl=flshow, fh=fhshow, col=pal2, ygrid=FALSE, AXE=c(1,2), CSCALE=CSCALE)
              ##  image(t(DEV$DSPEC[1:(DEV$numfreqs/2),], col=pal2))

          ## 

          anam = KH$info$fn[ipick]

          snam = unlist( strsplit(split="/", anam))

          inam = snam[length(snam)]
          ## title(main=paste(sep='    ',inam, KH$ftime[ipick]))
          mtext(paste(sep=' ',"Day=",theday, "Hour=", hr, KH$STNS[ipick], KH$COMPS[ipick] ), side=3, at=0, adj=0, col=3)
          kplot = kplot+1
          ###  to get thse to stay together on a page need to add some logic...
          
        }

          
          ##  dev.set(dev.next())

   
         if(kplot>=NPP)
            {
            kplot = longreset(NPP, PS)
            }

          KK = KK + 1
        }
    }

 if(PS==TRUE)
   {
     if(identical(names(dev.cur()), "postscript")) {  dev.off() }
   }
 par(def.par)
}

