`getphaselag2` <-
function(y1,y2, DT=0.008, frange=c(0,20), PLOT=FALSE, PLOT1=FALSE,  PLOT2=FALSE)
  {

############ dyn.load("/home/lees/Progs/Rc/MTAPSRC.so")
        if(missing(DT)) { DT=1 }
    if(missing(PLOT)) { PLOT=FALSE }

     if(missing(PLOT1)) { PLOT1=FALSE }

       
    if(missing(PLOT2)) { PLOT2=FALSE }
    if(missing(frange)) { frange=c(0,20) }

       ##   PLOT1 = FALSE
     ##   PLOT2 = FALSE

  n1 = length(y1)
  n2 = length(y2)
  n = max(c(n1,n2))

    k1 = next2(length(y1))
    k2  = next2(length(y2))
    KLEN = max(c(k1, k2))

  ts1 = c(y1-mean(y1), rep(0,n-n1))
  ts2 = c(y2-mean(y2), rep(0,n-n2))

    MTS1 =   mtapspec(y1 , DT , klen=KLEN,  MTP=list(kind=2,nwin=5, npi=3,inorm=0)  )
    MTS2 =   mtapspec(y2 , DT , klen=KLEN,  MTP=list(kind=2,nwin=5, npi=3,inorm=0)  )

#####names(MTS1)
#####dim(MTS1$Rspec)
    AA = MTS1$Rspec*MTS1$Rspec+MTS1$Ispec*MTS1$Ispec
    BB = MTS2$Rspec*MTS2$Rspec+MTS2$Ispec*MTS2$Ispec
    ABr = MTS1$Rspec*MTS2$Rspec+MTS1$Ispec*MTS2$Ispec
    ABi = MTS1$Rspec*MTS2$Ispec-MTS1$Ispec*MTS2$Rspec

    aa = apply(AA, 1, sum)
    bb = apply(BB, 1, sum)
    abre  = apply(ABr, 1, sum)
    abim  = apply(ABi, 1, sum)


    phase = atan2(abim,abre)

    abre = abre / (sqrt(aa) * sqrt(bb));
    abim = abim / (sqrt(aa) * sqrt(bb));

    cohere = sqrt((abre)^2 + (abim)^2);


###matplot(cbind(y1,y2), type='l', lty=1)

    nphase = LocalUnwrap(phase)

        
###plot(MTS1$freq, nphase)

    wgt = (cohere)^2 / (1.0- (cohere)^2 ) ;


###plot(MTS1$freq, cohere)


    flag =MTS1$freq>=frange[1] &  MTS1$freq<=frange[2]
    MOD = lm(nphase[flag] ~ MTS1$freq[flag] ,  weights=wgt[flag])

  
    phaselag = MOD$coefficients[2]/(2.0*pi);


if(identical(PLOT, TRUE))
  {
    par(mfrow=c(3,1))

    
    plot(MTS1$freq[flag], nphase[flag],  main=paste(sep=' ', "thelag=",phaselag) )
 plot(MTS1$freq[flag], phase[flag],  main=paste(sep=' ', "Unwrapped" ) )
 plot(MTS1$freq[flag], nphase[flag],  main=paste(sep=' ', "Unwrapped" ) )

    
    abline(MOD)
    ex = seq(from=0, length=length(ts1), by=DT)
    
    lagsamples  = floor(phaselag/DT)
    k = which.max(ts2)
    
    mlag = lagsamples
    gex = k-lagsamples
    
    plot(ex, ts1, type='l')
    abline(v=ex[gex], col=rgb(1,0,0), lty=2)
    plot(ex, ts2, type='l')
    
         if(gex>0 & gex<=n)
        {
          segments(ex[k], ts2[k], ex[gex], ts2[k], col=rgb(1,0,0))
          abline(v=ex[gex], col=rgb(1,0,0), lty=2)
        }
      else
        {
          segments(ex[n/2], ts2[k], ex[(n/2)+mlag/DT], ts2[k], col=rgb(1,0,0))
        }
  

  }
        

if(identical(PLOT1, TRUE))
  {
    par(mfrow=c(4,1))

   plot(MTS1$freq, nphase,  main=paste(sep=' ', "thelag=",phaselag) )
     abline(MOD)
    


    ##  plot.ts(cbind(z1,z2),  plot.type="multiple" , title=paste(sep=' ', "thelag=",phaselag) )
 
    lagplot(ts1, DT, phaselag, PLOT=TRUE )
    ex = seq(from=0, length=length(ts2), by=DT)
    
    plot(ex, ts2, type='l')
    

    }




  if(identical(PLOT2,TRUE))
    {
      screens(2)
      dev.set(dev.next())
      
      EX = DT*(1:n)
     plot(DT*(1:n), ts1, xlab="time, s", type = 'l')
     plot(DT*(1:n), ts2, xlab="time, s", type = 'l')
      k = which.max(ts2)
     xc = xcor2(ts1, ts2, DT, PLOT = FALSE)

      thelag = xc$mlag
      mlag = thelag 
      gex = k+mlag/DT
      if(gex>0 & gex<=n)
        {
          segments(EX[k], ts2[k], EX[gex], ts2[k], col=rgb(1,0,0))
          abline(v=EX[gex], col=rgb(1,0,0), lty=2)
        }
      else
        {
          segments(EX[n/2], ts2[k], EX[(n/2)+mlag/DT], ts2[k], col=rgb(1,0,0))
        }
      plot(xc)
      points(xc$lag[which.max(xc$acf)], max(xc$acf), col=2)
    }






    
    return(phaselag)
  }

