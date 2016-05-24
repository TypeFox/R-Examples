#dyn.load("Cbaselinefit.so")
#dyn.load("Fbaselinefit.so")

### default Werte: ###
#
### tau = 2.5 ### tau = 2.0 für mehr Peaks ###
#
### gam = 1 ###
#
### scl.factor = 1.2 ###
#
### maxwdth = 5 ###

"baselinefit" <- function(data,tau=2.5,gam=1, scl.factor=1.2, maxwdth=5){
  #Abfrage äquidistante Winkelwerte für 2Theta
  if(any(diff(data[,1], 1, 2) > 10^-10)) 
    stop("2 theta values have to be equidistant")
  #Abfrage ob NAs vorhanden
  if(any(is.na(data))) {
    stop("data contains NAs")
  }
  #Abfrage ob negative Werte vorhanden
  if(any(data[,2] < 0)) {
    stop("negative counts are not permitted")
  }
  #Tau muss im Interall [2,4] liegen
  if(!(tau >= 2 && tau <= 4)) {
    stop("tau must be between [2,4]")
  }  
  # Taut String Anpassung #
  # out=TRUE: Ausreisser werden eleminiert, dann muessen im folgenden pmg$y als Daten verwendet werden
  pmg <- fnpregm(data[,2],tau=tau,p=1.25,scl.factor=scl.factor, out=TRUE)
  # Abfrage ob Peak vorhanden
  if(pmg$pks == 0)
    stop("no peaks are found in data")
  # Bestimmung der Extremwerte der TS #
  expmg <- extrempmreg.fn(pmg)

  # Anpassung eines weighted smoothing Splines #
  spl <- wsspoiss(data[,1],pmg$y,tau=tau,sqfn=pmg$scl)
	if (spl$mem!=0) {
		print("Not enough memory for spline approximation!")
		return(0)
		}

  # Erweiterung der Extremwertbereiche aus TS mit spl #
  # gam wird verwendet in:  gr <- gamma*median(abs(diffwss))
  exb <- exber.maxwdth(expmg,spl,gamma=gam,xval=data[,1],xex=data[expmg$pkloc,1],maxwdth=maxwdth)
  
  if(length(exb$indl)<=2)
    indextr <- indextremw(exb)
  else
    indextr <- indextremw.c(exb)

# Anpassung der Basislinie mit WSS #
  bs <- wsspoiss(data[-indextr,1],pmg$y[-indextr],sqfn=pmg$scl[-indextr],thresh=spl$thresh)
	if (bs$mem!=0) {
		print("Not enough memory for spline approximation!")
		return(0)
		}
 
  # lineare Zwischenstuecke unter Peaks #
  basiserg <- basiserg(data[,1],pmg$y,bs$reg, exb$indl,exb$indr,indextr)

  npks <- pkcnt(exb$indlsep,exb$indrsep,pmg)

  pks <- rep(0,le=length(data[,1]))
  pks[indextr] <- data[indextr,2]-basiserg$basisl[indextr]
  
  list(pmg=pmg, spl=spl, baseline=basiserg,npks=npks, indlsep=exb$indlsep,indrsep=exb$indrsep, indextr=indextr,bs=bs,pks=pks,exb=exb, x=data[,1], y=data[,2])
}



"fnpregm" <- function(y,scl=0,tau=3,p=3,pks=-1,mult=TRUE,triang=FALSE,lcl=TRUE,run=FALSE,bmm=FALSE,fig=FALSE,verb=FALSE,eps=-1,dmax=0,out=FALSE,scl.factor=1)
{
	### INPUTS ###
	# data				y
        #	
	# scale to allow for            scl
	# heteroskedastic noise	
        #	
	# number of peaks               pks (0,1,2,3,...)
	# if pks=-1 the number
	# of peaks is determined
	# automatically. Otherwise
	# a function with pks peaks
	# is returned
        #
	# Multiresolution scheme	mult (TRUE, FALSE)
	# for interavls (TRUE) or
	# all intervals (FALSE)		
        #
	# Local or global squeezing	lcl  (TRUE,FALSE)
        #
	# Create plot of scale		fig  (TRUE,FALSE)
        #	
	# Controls size of		tau  (0 < tau ) 
	# thresholding	
        #
	#############	


  n<-length(y)	
  mny<-mean(y)
  y<-y-mny
  lscl<-length(scl)
  if(lscl==1){
    scl<-1.048358*median(abs(diff(y)))
    scl<-0*(1:n)+scl
  }
  else{
    scl<-scl
  }
  tmp <- .Fortran("npreg",
                  as.double(y),
                  as.double(scl),
                  double(n+1),
                  double(n+1),
                  double(n+1),
                  double(n+1),
                  double(n+1),
                  double(n+1),
                  as.double(p),
                  as.double(tau),
                  integer(n+1),
                  integer(n+1),
                  integer(n+1),
                  as.integer(n),
                  integer(1),
                  as.integer(pks),
                  as.logical(mult),
                  as.logical(lcl),
                  as.logical(triang),
                  as.logical(run),
                  as.logical(verb),
                  as.logical(bmm),
                  as.double(eps),
                  as.double(dmax),
                  PACKAGE="diffractometry"
                  )
  
  fn<-tmp[[6]][1:n]
  res<-y-fn
  tmps<-fcnst(res)
  if(lscl==1){
    scl<-pmax(scl.factor*tmps$fn,sqrt(fn+mny))
  }
  if(out){
    y<-frunmed(y,win=7,sig=scl)
  }
  
  tmp <- .Fortran("npreg",
                  as.double(y),
                  as.double(scl),
                  double(n+1),
                  double(n+1),
                  double(n+1),
                  double(n+1),
                  double(n+1),
                  double(n+1),
                  as.double(p),
                  as.double(tau),
                  integer(n+1),
                  integer(n+1),
                  integer(n+1),
                  as.integer(n),
                  integer(1),
                  as.integer(pks),
                  as.logical(mult),
                  as.logical(lcl),
                  as.logical(triang),
                  as.logical(run),
                  as.logical(verb),
                  as.logical(bmm),
                  as.double(eps),
                  as.double(dmax),
                  PACKAGE="diffractometry"
                  )	

  ###OUTPUTS###

  ### regression function ###
  fn<- tmp[[6]][1:n]+mny
  y<-y+mny
        
  ### number of local extremes ###
  pks<- tmp[[16]]

  ### indices of local extremes ###
  pkloc<-0
  if(pks>0){
    pkloc<-tmp[[11]][1:pks]
  }
        
  ### taut string with boundaries###
  ss<-tmp[[3]]   ### centre 
  ll<-tmp[[4]]   ### lower bound
  uu<-tmp[[5]]   ### upper bound
  nknts<-tmp[[15]]
  knts<-abs(tmp[[13]][1:nknts])
  str<-tmp[[8]][1:nknts]  ### taut string
  
  ###plot regression function###
  if(fig){
    plot(y,col="grey",ylab="Regression function")
    lines(fn,col="red",lwd=2)
  }
  list(y=y,fn=fn,pks=pks,pkloc=pkloc,scl=scl,ll=ll,uu=uu,ss=ss,knts=knts,str=str,tau=tmp[[11]],dmax=tmp[[24]])
}


"frunmed" <- function(x,win,fac=2.5,sig){
  n<-length(x)
  xx<-x
  k<-(win-1)/2
  i<-1
  while(i<=n){
    l1<-max(1,i-k)
    l2<-min(n,i+k)
    xm<-median(x[l1:l2])
    if(abs(x[i]-xm)>fac*sig[i]){xx[i]<-xm+fac*sig[i]*sign(x[i]-xm)}
    i<-i+1
  }
  return(xx)
}


"fcnst" <- function(x,alpha=0.5,fig=FALSE)
{
###INPUTS###
#
# x	     = data
#
# pks        = if >0 then this gives number of peaks	
#
# sqzf	     = squeeze factor for tube
#
# mult       = if 0 then multiresolution intervals, if 1 all intervals
#
# lcl        = if 0 local squeezing, otherwise global squeezing
#
# tau	     = threshold parameter
#
###
  n<-length(x)
  medx<- median(abs(x))		
  y<-abs(x)/medx
  y<-pmax(10^(-1),y)
  qxtrml<-double(n)
  qxtrmu<-double(n)
#
#	calculate cutoff value
#
  if(alpha==0) {
    alpha<--(1-1/(n*sqrt(pi*log(n))))
  }	
### fortran subroutine###
#	
  tmp <- .Fortran("fcnst",
                  as.double(y),
                  double(n),
                  double(n+1),
                  double(n+1),
                  as.double(qxtrml),
                  as.double(qxtrmu),
                  integer(n),
                  integer(1),	
                  as.integer(n),
                  as.double(alpha),
                  PACKAGE="diffractometry"
                  )
#	
###spectral density###
#
  fn<- tmp[[2]]*medx
  nind<-tmp[[8]]
  ind<-tmp[[7]][1:nind]
  soj<-diff(ind)
  vol<-fn[ind[1:(nind-1)]]
#	
###plot scale###
#
  if(fig){
    plot(abs(x),ylab="REGRESSION FUNCTION",col="grey")
    lines(fn,t="s")
  }
###OUTPUTS###
#
# fn	     = regression function
#
# pks	     = number of peaks
#
####
  list(fn=fn,ind=ind,nind=nind,soj=soj,vol=vol)

}






"extrempmreg.fn" <- function(ts)
  {
### INPUT ###
#
# ts: Ergebnis von fnpreg(m)
#
###
    tmp <- .C("extremfnpreg",
              as.double(ts$fn),
              as.integer(ts$pkloc),
              as.integer(ts$pks),
              as.integer(ts$knts),
              integer(ts$pks),
              integer(ts$pks),
              integer(ts$pks),
#              integer(le=ts$pks),
#              integer(le=ts$pks),
#              integer(le=ts$pks),
              PACKAGE="diffractometry"
              )
    list(exvalindl=tmp[[5]][tmp[[7]]>0], exvalindr=tmp[[6]][tmp[[7]]>0], minmax=tmp[[7]][tmp[[7]]>0],pkloc=ts$pkloc[tmp[[7]]>0])
  }




"wsspoiss" <- function(x,y,tau=2.3, weights,sigma, thresh, mr=1,glob=0,q=2.0,nit,fn, sqfn)
  {
    ### INPUT ###
    #
    # x: x-Koordinaten der Daten
    # y: y-Koordinaten der Daten
    #
    # fn bzw. sqfn: TS Anpassung zur Skalierung
    #
    ####
    print(c("q:", q))
    n <- length(y)
    if(length(x)!=n)
      stop("x and y length differ")
    b <- double((n-1))
    for (i in 1:(n-1)){
      b[i]<- (y[i]-y[i+1])/sqrt(2)
    }
    if(missing(sigma))
      sigma <- mad(b,center=0)
    print(c("sigma:", sigma))
    if(missing(thresh))
      thresh <- sigma*sqrt(tau*log(n))
    print(c("thresh:", thresh))
    if(missing(weights)){
      w <- 1/((max(x)-min(x))^2*n*100)
      weights <- rep(w,le=n)
    }
    else{
      if(length(weights)!=n)
        stop("y and weights lengths differ")
    }
    if((mr==2)&&(missing(nit))){
      print("Number of iterations not specified!")
      break
    }
    else
      if(missing(nit))
        nit=1
    if(missing(sqfn)){
      sqfn <- sqrt(fn)
      sqfn[sqfn<sigma] <- sigma
    }
    weights <- weights/(sqfn^2)
    if((max(y)/mad(y,center=0))>5) shrtint <- 1
    else
      shrtint <- 0
    XXX <- .C("wsspoisschngd",
              as.double(x),
              as.double(y),
              as.integer(n),
              as.double(weights),
              as.double(thresh/sigma),
              reg=double(n),
              as.integer(glob),
              as.integer(mr),
              as.integer(nit),
              as.double(q),
              as.double(sqfn),
              as.integer(shrtint),
              as.integer(0),								## Memory
              PACKAGE="diffractometry")
    list(weights=XXX[[4]], reg=XXX[[6]], thresh=thresh, sigma=sigma,sqfn=sqfn,mem=XXX[[13]])
  }



"exber.maxwdth" <- function(x,y,datx,gamma=1,maxwdth=5,xval,xex)
  {
### INPUT ###
#
# x: output von extrempmreg.fn fuer TS
# y: output von wsspoiss
# maxwdth: maximale Peakbreite, Standardwert= 5°
# xval: x-Koordinaten der Daten
# xex: x-Koordinaten der in extrempmreg.fn bestimmten Peakintervall
#
###
    nex <- length(x$exvalindl)
    n <- length(y$reg)
    if(missing(datx))
      datx <- 1:n
    diffwss <- diff(y$reg)/diff(datx)
    dx <- (datx[-1]+datx[-n])/2
    diff2wss <- diff(diffwss)/diff(dx)
    indl <- c(1,x$exvalindl,n-1)
    indr <- c(1,x$exvalindr,n-1)
    minmax <- c(0,x$minmax,0)
    print(c("Median", median(abs(diffwss))))
    print(c("MAD", mad(abs(diffwss))))
    gr <- gamma*median(abs(diffwss))
    exent <- c(0,rep(1,le=nex),0)
    tmp <- .C("exber_maxwdth",
              as.double(y$reg),
              as.double(diffwss),
              as.double(gr),
              as.integer(n),
              as.integer(minmax),
              as.integer(nex),
              as.integer(indl),
              as.integer(indr),
              as.integer(exent),
              as.double(maxwdth),
              as.double(xval),
              as.double(xex),
              PACKAGE="diffractometry")
    tmpex <- tmp
    exent <- tmp[[9]]
    indl <- tmp[[7]][exent>0]
    indr <- tmp[[8]][exent>0]
    nex <- length(indl)
    if(length(indl)>1){
      tmp <- indl[-1]-indr[-nex]
      if(min(tmp)<15){
        indlmrgd <- indl[-(2:nex)[tmp<15]]
        indrmrgd <- indr[-(1:(nex-1))[tmp<15]]
        if(min(tmp)<2){
          indl <- indl[-(2:nex)[tmp<2]]
          indr <- indr[-(1:(nex-1))[tmp<2]]
        }
      }
      else{
      indlmrgd <- indl
      indrmrgd <- indr
    }
    }
    else{
      indlmrgd <- indl
      indrmrgd <- indr
    }
    list(indlsep=indl, indrsep=indr,indl=indlmrgd,indr=indrmrgd)#,midmerg=midmerg,imerg=imerg)
  }



"indextremw" <- function(index)
  {
    nex <- length(index$indl)
    ind <- 0
    for(i in 1:nex)
      ind <- c(ind,index$indl[i]:index$indr[i])
    ind[-1]
  }



"indextremw.c" <- function(index)
  {
    nex <- length(index$indl)
    ind <- rep(0,le=index$indl[nex])
    tmp <- .C("indextremw",
              as.integer(nex),
              as.integer(index$indl),
              as.integer(index$indr),
              as.integer(ind),
              PACKAGE="diffractometry")
    nmaxind <- tmp[[1]]+1
    tmp[[4]][1:nmaxind]
  }



"basiserg" <- function(data.x, data.y,anp,pindl,pindr,ind)
  {
    ### INPUT
    #
    # data.x: x-Koordinaten der Daten
    # data.y: y-Koordinaten der Daten
    # anp: Basislinienanpassung
    # pindl: linke Indizes der Extremwertintervalle
    # pindr: rechte Indizes der Extremwertintervalle
    # ind: Indizes der Extremwerte
    #
    ###
    np <- length(pindl)    # Anzahl Peaks
    n <- length(data.x)      # Datenlänge
    basisl <- double(n)
    basisl[-ind] <- anp
    peaks <- double(n)
    peaks[-ind] <- 0
    indl <- integer(np)
    indr <- integer(np)
    tmp <- .C("basiserg",
              as.integer(np),
              as.integer(n),
              as.double(data.x),
              as.double(data.y),
              as.double(basisl),
              as.double(peaks),
              as.integer(indl),
              as.integer(indr),
              as.integer(pindl),
              as.integer(pindr),
              PACKAGE="diffractometry"
              )
    indl <- tmp[[7]][tmp[[7]]>0]
    indr <- tmp[[8]][tmp[[8]]>0]
    list(basisl=tmp[[5]], peaks=tmp[[6]], indl=indl, indr=indr)
  }


"pkcnt" <- function(indl,indr,ts)
  {
    ### INPUT ###
    #
    # indl: linke Indizes der Peakintervalle
    # indr: rechte Indizes der Peakintervalle
    # ts: output von fnpreg(m)
    #
    ###
    nind <- length(indl)
    tmp <- .C("pkcnt",
              as.integer(indl),
              as.integer(indr),
              as.integer(nind),
              as.integer(ts$pkloc),
              as.integer(ts$pks),
              as.double(c(ts$fn[1],ts$fn[ts$pkloc[1]])),
#              integer(le=nind),
              as.integer(rep(0,nind)),
              PACKAGE="diffractometry")
    tmp[[7]]
  }
