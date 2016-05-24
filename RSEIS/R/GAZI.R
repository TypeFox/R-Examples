`GAZI` <-
function(ADAT, dt=1, ex=seq(0, 100), comp=c(4,5,6), sta="ZZZ", az=0, len=50, shift=10, prev=1, fileid="", picks=NA, labs=NA )
{
##X##    particle motion estimator
##X## INPUT:
##X##     ADAT = N by 3 matrix of three component seismic data
##X##     dt = deltat  sample interval (s)
##X##     ex = time values for X axis
##X##     comp  = components
##X##     sta = station name
##X##     az = azimuth from source to reciever
##X##     len = length of sliding window (samples)
##X##     shift = shift of each window (samples)
##X##     prev = pre-event window for calculation of noise average
##X##     fileid = character string identifying the event (yy:mm:hh:mi:se, for example)
##X##     picks = vector of times relative to the start for lines to be drawn
##X##     labs = vector of labels for the picks lines
  

       if(missing(dt)) { dt=1 }
       if(missing(ex)) { ex=seq(from=0, length=dim(ADAT)[1], by=dt) }
       if(missing(len)) { len=50 }
       if(missing(shift)) { shift=10 }
       if(missing(comp)) { comp=c(4,5,6) }
       if(missing(sta)) {  sta="ZZZ"}
       if(missing(az)) {  az=0 }
       if(missing(prev)) {  prev=1 }
       if(missing(fileid)) { fileid=" " }
       if(missing(picks)) { picks=NA  }
       if(missing(labs)) {  labs=NA  }




	opar=par(no.readonly = TRUE)
	alen=length(ADAT[,1])
	dt=dt
	ex = ex

	dat = ADAT

	ascd = ADAT
	
  	aex=rep(0,alen)
 	aaz=rep(0,alen)
  	ai=rep(0,alen)
  	rateig=rep(0,alen)

	winlen=len*dt
 	winn=len
  	winhalf=winn/2
	k = winn/2
   	wincen=ex[k]-ex[1]	
   	wina=(wincen/dt)-winhalf
    	winb=(wincen/dt)+winhalf
  
  	k=len/2
  	j=1
  

  # 	xtics=pretty(seq(from=min(ex), to=max(ex), N=10))
  	xtics=pretty(ex, n=10)


	mintic=min(xtics)
	maxtic=max(xtics)

  	# for each trace, find pre-event DC offset and remove that from the whole trace
	# do not remove the mean again below, that would be wrong

  	ax=1:length(ex)

	#   here we determine a limit on X

  	flagax = ax[ex<prev]
  	tem=dat[ flagax ,]
  	mns=apply(tem,2,mean)
  	dtem=sweep(dat, 2, mns)

  while(k<(alen-len/2))
    {
      wincen=ex[k]-ex[1]
      wina=round((wincen/dt)-winhalf+1)
      winb=round((wincen/dt)+winhalf+1)
      
      winb=min(winb,alen)
      #    print(c(wina, winb))

      tem=dtem[wina:winb,]
					# need to remove the mean value from each column (we did this above)
					#    NO:  tem=sweep(tem, 2, apply(tem,2,mean))


      
      
     ###### covtem=t(tem) %*% (tem)

      covtem = var(tem)

      
      eg=eigen(covtem, symmetric = TRUE )
					# Be=winn*diag(1,nrow=3) + matrix(c(-1,1,1,1,-1,1,1,1,-1),nrow=3)*covtem
					# Beg=eigen(Be)
      
					# Kappa<-log(Beg$values[1]/Beg$values[2])/log(Beg$values[2]/Beg$values[3])
      
      
      aex[j]=ex[k]
      ## rateig[j]=sqrt( eg$values[2]^2 + eg$values[3]^2 ) / eg$values[1]


      ##  Joydeep recommends using the following measure of rectilinearity
      ## jepsen and kennett, 1990, bssa, 80b, #6, 2032-2052.

      rateig[j]= 1 - ((eg$values[2]+eg$values[3])/(2*eg$values[1]))


      
					#  rateig[j]=Kappa

      #   careful here: be sure the azimuth below is calculated in the N-E-Down coordinate system
#   1=Z   2=N   3=E
#  this means that the real azimuth is 90-alpha  where alpha is the counter-clockwise
#  coordinate angle derived below
      RAD2DEG = 180/pi
      alpha=RAD2DEG*atan2(eg$vectors[2,1], eg$vectors[3,1])

     ###    az=90-alpha
      az=alpha


      inci=RAD2DEG*atan2(eg$vectors[1,1], sqrt(eg$vectors[2,1]^2+eg$vectors[3,1]^2))


### convert angles so that they are orientations as shiftnd not simply directions
###  this is because the direction is irrelevant and -10deg=170deg orientation

     ### if(az<0) az=az+180

      aaz[j]=az
      if(inci<0)inci=abs(inci)
      ai[j]=inci
      
      k=k+shift
      j=j+1
    }
       
  jall=j-1
  
#  dev.set(which=2)
#########   old: par(mfrow=c(6, 1) )
       
       
  par(mfrow=c(6, 1) )
  par(mai=c(0.1, .5, 0.1, 0.5) )
  for(i in 1:3)
    { 
      plot(ex,dat[,i], axes=FALSE, xlab="",ylab="", type="n")
      lines(ex,dat[,i],type="l")
      axis(1,tck=.03,at=xtics,labels=FALSE)
      # axis(2, las=1)
      axis(3,tck=.03,at=xtics,labels=FALSE)
      box()
      locy=0.8*max(ascd[,i])

      tcomp = fixcompname(comp[i])
      text(ex[1], locy,paste(sta,tcomp,sep=" : ") ,cex=.8, adj=0)


      if(!is.na(picks)) { PLTpicks(picks, labs) }
      
       ###  plot.ps(ain)    ###  plots the P and S lines on the graph
        ###   plot.t1t2(ain)
      
      letter.it(i,2)
      
    }	
  i=3
  locy=0.8*max(ascd[,i])
					#  locy=0.95*min(dat[,i])
					# text(max(ex), locy, paste(ain$fil, ain$pfil,ain$id, ain$sec,sep=" : ") , cex=.8,  adj=1, col=3)
  
#######  NOW plot New Stuff  ############################
##  this switches to the other opened window

##    dev.set(which=3)
##  par(mfrow=c(3, 1) )
##   par(mai=c(0.0, .5, 0.1, 0.5) )
  par(mai=c(0.1, .5, 0.1, 0.5) )



#####    INC ANGLE
  
  plot(aex[0:jall],ai[0:jall],xlim=range(ex),ylim=c(0,90),type="n",  axes=FALSE, xlab="",ylab="IncAng, deg")
        lines(aex[0:jall],ai[0:jall],type="l")
					#   abline(h=c(0))
					#   axis(2, las=1)
 #  axis(2, at=c(-60, -30,0,30 , 60), tck=1, las=1, lty=2, lwd=0.5)
axis(2, at=seq(0,90, by=10), tck=1, las=1, lty=2, lwd=0.5)
  axis(3,at=xtics,tck=.03,labels=FALSE)

   if(!is.na(picks)) { PLTpicks(picks, labs) }
       
  ### plot.ps(ain)
  ### plot.t1t2(ain)
  
   letter.it(4,2)
  box()

  incpar<-par()

  figinc=par("fig")
  
####  RATIO
  

  plot(aex[0:jall],rateig[0:jall],xlim=range(ex),type='n',  axes=FALSE, xlab="",ylab="RatEig")
        lines(aex[0:jall],rateig[0:jall], type='l')
  locy=0.8*max(rateig[0:jall])
  
  axis(2, las=1)
  axis(1,at=xtics,tck=.03, las=1,   mgp=c(.1,.1,0))
  axis(3,at=xtics,tck=.03,labels=FALSE)
  mtext( paste(fileid) , line=0.1)
     if(!is.na(picks)) { PLTpicks(picks, labs) }
   ### plot.ps(ain)
   ### plot.t1t2(ain)
  box()
   letter.it(5,2)
#####   Azimuth
   par(mai=c(0.2, .5, 0.15, 0.5) )

   azims =     RPMG::fmod(aaz[1:jall], 180)

  plot(aex[1:jall], azims ,xlim=range(ex),ylim=c(0,180),axes=FALSE, xlab="Time, s",ylab="Az, deg")

 # axis(2, at=c(-150,-100, -50,0,50, 100, 150), tck=1, las=1, lty=2, lwd=0.5)
  axis(2, at=seq(0,180, by=20), tck=1, las=1, lty=2, lwd=0.5)
  axis(3,at=xtics,tck=.03,labels=FALSE)
   axis(1,at=xtics,tck=.03, las=1,   mgp=c(.1,.1,0))

  AZ=  RPMG::fmod(az, 180)

###    if(AZ>180) AZ<-AZ-360
					# abline(h=c(0))
  abline(h=c(AZ),lty=4, col=2)
    #   locy=0.8*max(aaz[0:jall])
  locy=165
  text(max(xtics), locy, paste("AZ=",format.default(AZ, digits=3)) ,adj=1, cex=1.2, col=2)
  box()
  
   if(!is.na(picks)) { PLTpicks(picks, labs) }
  ### plot.ps(ain)
 ###   plot.t1t2(ain)

 #  plot.medbars(ain,aex[1:jall], aaz[1:jall]   )

   figaz=par("fig")
   usraz=par("usr")
  
     locy=0.9*min(aaz[0:jall])
 	#   locy=-165
  m=max(aex[0:jall])
  segments( m-winlen, locy, m , locy, lwd=3)
  
   letter.it(6,2)
  azpar<-par()
  # dev.prev()

  #  invisible(par(opar))
  #
	newpar=par(no.readonly = TRUE)
       par(opar)
  invisible(list(aex=aex[1:jall], rateig=rateig[1:jall], aaz=aaz[1:jall], ai=ai[1:jall], figaz=figaz, azpar=azpar, incpar=incpar, par=newpar  )  )	
}

