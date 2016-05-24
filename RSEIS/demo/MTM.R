
###############  reproduce the figures, more or less from Lees and Park, 1995


##############  Slepian Tapers:
 nwin = 5
     npi = 3
     npoints = 900
     sleps = get.slepians(npoints, nwin, npi)

     matplot(sleps, type='l', xlab="Index", ylab="Taper Amplitude")
     legend('topleft', legend=1:nwin, lty=1:nwin, col=1:nwin)
title("Slepian Tapers: Figure 1, Lees and Park, 1995")
readline("To Continue Hit Enter Key\n")
##############
##############  EXAMPLES:



f1 =20
f2 = 30
dt = 0.01

t = seq(from=0, to=6, by=dt)


noise =  runif(length(t), 0, 2)

y = 2*sin(2*pi*f1*t) + 3*sin(2*pi*f2*t) + noise 

par(mfrow=c(1,1))
plot(t,y, type='l')

title("Figure 2(synthetic example 1): Lees and Park, 1995")


readline("To Continue Hit Enter Key\n")
###############
##########################  make a function to get the single taper spectrum

singleTaper<-function(y, dt, tappercent=0.1 )
  {
    if(missing(tappercent)) tappercent=0.1
    N= length(y)	
    fn = 1/(2*dt)
    tapy = rsspec.taper(y, p=tappercent)
    ##tapy = tapy-mean(tapy)
    
    Y = fft(tapy)
    Pyy = (Mod(Y)^2)/(N*N)
    ##  Pyy = Y * Conj(Y)
    n = floor(length(Pyy)/2)
    Syy = Pyy[1:n]
    f = (0:(length(Syy)-1))*fn/length(Syy)
#####    plot(f, Syy, type='l', xlab="frequency", ylab="Power Density", log='')
    
    invisible(list(f=f, syy=Syy))
    
  }
#######################  make a function to create the plots
DOLEESPARK<-function(y, dt, tappercent=0.1)
  {

    y = y-mean(y)

####

    nn=next2(length(y))

    Mspec =   mtapspec(y, dt, klen=nn,  MTP=list(kind=2,nwin=5, npi=3,inorm=1)  )
    f = Mspec$freq

    amp1 = Mspec$spec[1:length(f)]
    amp = 10*log10(amp1 )

#### 
    singy = singleTaper(y, dt, tappercent=tappercent)
    
    stap1 = 10*log10(singy$syy )

#### 
    squig = list(y=y, dt=dt)
    ZIM = autoreg(squig , numf=length(Mspec$freq) , pord = 80, PLOT=FALSE)

    AutoR = 10*log10(ZIM$amp )

#### 
    frange = range(c( f ))
    amprange  = range(c(amp, stap1, AutoR))

####
   ## dev.new()

    
    opar <- par(no.readonly = TRUE)

    par(mfrow=c(3,1))
    
    plot(frange, amprange, type='n',ylab="Spectrum Amplitude", xlab="Frequency")


     kdof = 2*Mspec$mtm$nwin-2
    ppoints  =  c(50.0, 90.0, 95.0, 99.0, 99.5, 99.9)
    myf = qf(ppoints/100, 2, kdof)

    hivals1 = Mspec$freq[which(Mspec$Fv[1:length(Mspec$freq)]>myf[3] & Mspec$dof[1:length(Mspec$freq)]>9.0  )]
    hivals2 = Mspec$freq[which(Mspec$Fv[1:length(Mspec$freq)]>myf[4] & Mspec$dof[1:length(Mspec$freq)]>9.0  )]

    lightgreen = rgb(.8,1,.8)
    darkgreen  = rgb(.1, .8, .1)

    
    abline(v=hivals1, col=lightgreen , lty=2)
    abline(v=hivals2, col=darkgreen, lty=2)

    
    lines(f, amp)

    lines(ZIM$freq,AutoR, col='red')
    lines(singy$f,stap1, col='blue')
#### 
    title("Lees and Park, 1995")
    legend("topright", legend=c("multitaper",  "Autoregressive", "Single-taper"), lty=c(1,1,1), col=c("black", "red", "blue"))

###   show diagnostics from MTM estimate:



    plot(Mspec$freq,Mspec$dof[1:length(Mspec$freq)],type='l',xlab="Frequency",ylab="Effective Degrees of Freedom")
      abline(v=hivals1, col=lightgreen, lty=2)
     abline(v=hivals2, col=darkgreen, lty=2)

    plot(Mspec$freq, Mspec$Fv[1:length(Mspec$freq)],type='l',xlab="Frequency",ylab="F-test")
    abline(v=hivals1, col=lightgreen, lty=2)
     abline(v=hivals2, col=darkgreen, lty=2)

    u=par("usr")
#####  see Percival and Walden p. 499-500 for degrees of freedom
   
    abline(h=myf,lty=2 )
    text(u[2],myf,ppoints, pos=4, xpd=TRUE)

    par(opar)
  }

###############  show various spectrum estimates


DOLEESPARK(y, dt, tappercent=0.05)

readline("To Continue Hit Enter Key\n")

###############   next example shows 1000
###   traces aded together each with the same amplitude spectrum, shifted by a linear phase
###  equivalent to one sample.

dt = 1
cut =  0.125
NF = 1024*8
ff = seq(from=0, to=1, length=NF)

mf = rep(1, NF)
mf[ff>0.125&ff<=0.875] = 0

K = (NF/2)

j = K/2
phz1 = 2*pi*j*ff
FIL = runif(1)* mf * as.complex(cos(phz1), sin(phz1) )

G = FIL
ag = fft(G, inverse = TRUE) / (NF/2)
y = Re(ag)[1:(NF/2)]

####  plot(seq(from=0, length=(NF/2), by=1), y, type='l', xlab='time', ylab="amplitude")

####  readline("To Continue Hit Enter Key\n")
v=NF/4

N = 1000
MY = rep(0, N)
for(i in (v-1000):(v) )
{
  i1 = i
  i2 =  i1+N-1
 #### print(paste(i1,i2, length(i1:i2) ))
  MY =   MY+runif(1)*y[i1:i2]
}

MY = MY/N
MY = detrend(MY)

plot(seq(from=0, length=length(MY), by=1), MY, type='l', xlab='time', ylab="amplitude")
title("Input Signal 2")



readline("To Continue Hit Enter Key\n")




dt2=1
DOLEESPARK(MY, dt2, tappercent=0.05)
readline("To Continue Hit Enter Key\n")
#######################
########  apply same method to the Delta-O18 data
data(OH)

plot( seq(from=0, length=length(OH$JSTR[[1]]), by=3000), OH$JSTR[[1]], type='l', xlab="time, yrs", ylab="Delta O18")
title("Delta O18")

readline("To Continue Hit Enter Key\n")


DOLEESPARK(OH$JSTR[[1]], 3, tappercent=0.05)

###title("Delta-O18")


readline("DONE: To Continue Hit Enter Key\n")
