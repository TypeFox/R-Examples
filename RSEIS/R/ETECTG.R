`ETECTG` <-
function(GH, sel=sel, FRWD=8,  BKWD=8,  sbef=1, saft=6, DFRWD=.5,  DBKWD=.5, thresh=2,  Tthresh2=7, stretch=1000, flo=0.1, fhi=5.0, PLOT = FALSE, Kmin=7, perc=0.05, kind=1,  DOARAIC = FALSE )
{ 
#######    do automatic picking on 3 component traces
#######  fnames4, fnames5, fnames6=names of files to read in vertical north east
#######  FRWD=forward window in seconds
#######  BKWD=backward window in seconds
#######  sbef=before window in seconds
#######  saft=back window in seconds
#######  thresh=threshhold
#######  Tthresh2=window threshhold
#######  stretch=stretching factor to multiply prior i=to rat curve
#######  flo=low freq cutoff for band pass
#######  fhi=high freq for  bandpass 
  if(missing(FRWD)) { FRWD=8 }
  if(missing(BKWD)) { BKWD=8 }
  if(missing(DFRWD)) { DFRWD=.5 }
  if(missing(DBKWD)) { DBKWD=.5 }

  if(missing(sbef)) { sbef = 1 }
  if(missing(saft)) { saft = 6}
  if(missing(thresh)) { thresh=2 }
  if(missing(Tthresh2)) { Tthresh2= (sbef+saft) }
  if(missing(stretch)) { stretch=1000 }
  if(missing(flo)) { flo = .1 }
  if(missing(fhi)) { fhi=5.0 }
  if(missing(PLOT)) {  PLOT = TRUE }
  if(missing(Kmin)) { Kmin = 7 }
  if(missing(kind)) { kind=1 }
  if(missing(sel)) { 1:length(GH$JSTR) }
  if(missing(perc)) { perc=0.05 }
  if(missing(DOARAIC)) { DOARAIC = FALSE }

  N1 = length(GH$JSTR)
  
  PTIMES = as.list(1:N1)


  ####  main loop

  for(M in 1:length(sel) )
    {
      M3 = sel[M]
      
      nam3 = paste("ay",M, sep="")
      
      assign(nam3, GH$JSTR[[M3]] )
      g = get(nam3)
      g[is.na(g)] = mean(g, na.rm=TRUE)
      
      nam4 = paste("PP",M, sep="")
      
      deltat = GH$info$dt[M3]
pik1 = pickit(g , deltat=deltat,  FRWD=FRWD,  BKWD=BKWD,sbef=sbef, saft=saft,
                          thresh=thresh, Tthresh2 =  Tthresh2, flo=flo, fhi=fhi, stretch=stretch, Kmin=Kmin)

      
      assign(nam4, pickit(g , deltat=deltat,  FRWD=FRWD,  BKWD=BKWD,sbef=sbef, saft=saft,
                          thresh=thresh, Tthresh2 =  Tthresh2, flo=flo, fhi=fhi, stretch=stretch, Kmin=Kmin))



      
      v = get(nam4)
      v$STNS = GH$STNS[M3]
      v$COMPS = GH$COMPS[M3]

      assign(nam4, v)

      if(PLOT==TRUE)
        {
###  plotting
          par(mfrow=c(2,1))
      
          plot.ts(v$ay)
          winmark(v$a1,v$a2,col=4)	
          plot.ts(v$RAT)

          winmark(v$a1,v$a2,col=4)
          abline(h=v$thresh, col=3)
          locator(1)
        }
      
###     ls(pat='PP')
###	rm(list=ls(pat='PP'))

    }

  ###   end initial detection


  
  for(M in 1:length(sel) )
    {
      nam4 = paste("PP",M, sep="")
      v = get(nam4)
      print(paste(" ", "ETECTG",  v$STNS, v$COMPS, length(v$RAT)) )
      nam5 = paste("tt",M, sep="")
      assign(nam5, rep(0, length(v$x)))
      tt  = get(nam5)
      for(j in 1:length(v$a1))
        {
          tt[v$x>=v$a1[j]&v$x<=v$a2[j]] = 1
        }
      assign(nam5, tt)
###     ls(pat='tt')
###	rm(list=ls(pat='tt'))

      
    }	
      ###############  create a series of ones for hits on the STLT algor

     ###############
      ## weights: might weigh certain stations more? 

  K1 = length(sel)
      
    ###   ALLP = tt4+tt5+tt6
      nam5 = paste("tt",1, sep="")
      ALLP  = get(nam5)
      if(K1>1)
        {
          for(M in 2:K1)
            {
              M3 = M
              nam5 = paste("tt",M3, sep="")
              tt = get(nam5)
              ALLP = ALLP+tt
            }
        }
  
   ## plot.ts(ALLP)
  
                                        #  weight the vertical more than the horizontals
      ## ALLP = 2.0*tt4+ tt5+ tt6
      
      JJ = Thresh.J(ALLP, K1-0.5)
###   JJ = Thresh.J(ALLP,1.5)
  NNALLP = length(ALLP)
      NJ = length(JJ$J)
      if(NJ<1)
        {

          detpix=0
          next;
        }


  for(M in 1:length(sel) )
    {
      K = sel[M]
      deltat = GH$info$dt[K]
      namdp = paste("DP",M, sep="")
      assign( namdp , rep(0,length(JJ$J)))
      detpix = get(namdp)
      namtp = paste("TP",M, sep="")
      
      nam4 = paste("PP",M, sep="")
      PP = get(nam4)

      ifrwd  =  round(DFRWD/PP$deltat)
      ibkwd  =  round(DBKWD/PP$deltat)    

      zwin = max(c(ifrwd, ibkwd))
###  for each PICK-WINDOW get a detailed pick on each component

      for(j in 1:NJ)
        {
          ## print(paste(sep=" ", "****************** sub win=", i, j, "of ", NJ, BIGN))
          ##  b1 = (JJ$J[j]-sbef/deltat)
          ##  b2 = (JJ$L[j]+saft/deltat)

          B = getb1b2(JJ$J[j], JJ$L[j], zwin, 375,  NNALLP )
          b1 = B[1]
          b2 = B[2]

          if(is.na(b1) | is.na(b2) ) { next }
          
          z4= PP$fy[b1:b2]
          
          xz = PP$x[b1:b2]

          if(length(z4)<375)
            {
              
              print(paste(sep=" " , "Very short window",M, j, length(z4), b1, b2)   )
              detpix[j] = 0
              next
            }
####  here do the detailed picking to find a good first arrival
####     RATP = ratcurve(z4, dt=PP4$deltat, fwlen =  75,  bwlen  =200, PLOT=TRUE)
          print("ETECTG Going to PSTLTcurve")
          RATP = PSTLTcurve(z4, dt=PP$deltat, fwlen=ifrwd,  bwlen=ibkwd, perc=perc, stretch=1000 , MED=77, PLOT=FALSE)
          if(is.na(RATP[[1]]) )
            {
              next;
            }
          ## locator(1)
          T1=TFIN=RATP$eye;

          if(DOARAIC==TRUE)
            {
              Nz4 = length(z4)
              Mar = 8
              O1=2; 
              O2=0.2; 
              WW=2
              aout = rep(0, Nz4)
              deltat = PP$deltat[1]
              plot.ts(z4)
              locator(1)

              ary = .C("CALL_ARAIC", PACKAGE = "RSEIS",
                as.double(z4), as.integer(Nz4),as.double(deltat), as.integer(Mar),
                as.integer(T1), as.double(O1), as.double(O2), as.double(WW), as.double(aout)) 

              kaic = ary[[9]]
              kaic[kaic==0]=NA
              Taic =TFIN=  which.min(kaic)

              xkaic = 1:length(kaic)
              plot(xkaic,kaic, type='l')
              lm1 = lm(kaic ~ cbind(xkaic, xkaic^2, xkaic^3, xkaic^4))
              lines(xkaic[!is.na(kaic)], lm1$fitted.values, col=2)

              vline(c(RATP$eye, RATP$ind, RATP$mix), per=-1,  LAB=c("eye", "ind", "mix"), COL=c(2,3,5))
              vline(which.min(kaic), COL=rgb(.4,.8,1) )
            }
          
          detpix[j] = xz[1]-1+TFIN
          
        }


      
      assign( namdp , detpix)
      detpix=get(namdp)
      PTIMES = recdate(GH$info$jd[K], GH$info$hr[K], GH$info$mi[K],
        GH$info$sec[K]+GH$info$msec[K]/1000+GH$info$t1[K]+detpix*GH$info$dt[K])
      PTIMES$yr = rep(GH$info$yr[K], length=length(detpix))
      PTIMES$STAID = list(stn=GH$STNS[K] , comp=GH$COMPS[K])
      
      assign( namtp , PTIMES)

      
    }


###     ls(pat='^TP')
###	rm(list=ls(pat='^TP'))

  
  ###     ls(pat='^DP')
###	rm(list=ls(pat='^DP'))

  
  if(PLOT==TRUE)
    {
      print(paste("going to plotting", K1))
      
      ##  PLOT.SEISN(GH, sel = sel, notes=GH$KNOTES[sel]);  gin = locator(2)
      ##  gin = locator(2)


###  swig(GH, sel=sel)
###  y = swig(GH, sel=sel, WIN=gin)

    ###  gin = locator(2)
      gin = NULL

      one()
      PLOT.SEISN(GH, sel = sel, notes=GH$KNOTES[sel], WIN=gin)

      
      for(M in 1:length(sel) )
        {
          K = sel[M]
          namtp = paste("TP",M, sep="")
          tp = get(namtp)

          ypos = (length(sel)-M+0.5)/length(sel)

          zloc = list(x=rep(NA, length(tp$jd)),  y=rep(ypos, length(tp$jd))  )
          
          zloc$x = secdif(GH$info$jd[K], GH$info$hr[K], GH$info$mi[K], GH$info$sec[K]+GH$info$msec[K]/1000+GH$info$t1[K],
            tp$jd, tp$hr, tp$mi, tp$sec)
          
          
          PPIX(zloc, YN=length(sel), col=4, lab='P')
        }

      

    }

      PPALL = as.list(1:K1)
      PPTIM = as.list(1:K1)

      for(M in 1:K1)
        {
          M3 = sel[M]
          nam2 = paste("PP",M, sep="")
          P = get(nam2)
          PPALL[[M]] = P

          namtp = paste("TP",M, sep="")
          ptims = get( namtp)
          PPTIM[[M]] = ptims

        }
  
    return(list(sel=sel, JJ=JJ, PPTIM=PPTIM, PP=PPALL))
  
}

