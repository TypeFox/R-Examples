################################################################################
##  seewave by Jerome Sueur, Caroline Simonis & Thierry Aubin
##  http://rug.mnhn.fr/seewave/
################################################################################

################################################################################
##                                ACI
################################################################################

ACI <- function (
    wave,
    f,
    wl = 512,
    ovlp = 0,
    wn = "hamming",
    flim = NULL,
    nbwindows = 1
    )
    
{
    ## INPUT
    input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; rm(input)

    ## SPECTROGRAM MATRIX
    z <- sspectro(wave, f, wn=wn, ovlp=ovlp)

    ## FREQUENCY SELECTION 
    if(!is.null(flim))
        {
            flim<-flim*1000*wl/f
            z <- z[flim[1]:flim[2], ]
        }

    acis <- array(0, nbwindows)

    for(j in 1:nbwindows)
        {
            ## WINDOWING
            l <- dim(z)[2]
            z1 <- z[,floor(l/nbwindows*(j-1) + 1):floor(l/nbwindows*j)]

            ## ACI COMPUTATION  
            t <- array(0, dim(z1)-c(0,1))
            for (i in 1:(dim(t)[1])) {t[i,] <- abs(diff(z1[i,]))/sum(z1[i,])}
            acis[j] <- sum(t)
        }

    return(sum(acis))
}


################################################################################
##                                ACOUSTAT
################################################################################

acoustat <- function (
            wave,
            f,
            wl = 512,
            ovlp = 0,
            wn = "hanning",
            tlim = NULL, 
            flim = NULL, 
            aggregate = sum,
            fraction = 90,
            plot = TRUE,
            type = "l",
            ...) 

{
    ## INPUT
    input <- inputw(wave = wave, f = f)
    wave <- input$w
    f <- input$f
    if (!is.null(tlim)) {
        wave <- cutw(wave, f = f, from = tlim[1], to = tlim[2])
    }
    n <- length(wave)
    rm(input)

    ## STFT
    m <- sspectro(wave, f = f, wl = wl, ovlp = ovlp, wn = wn)

    ## FREQUENCY SELECTION AND FREQUENCY AXIS
    if (is.null(flim)) {freq <- seq(0, (f/2) - (f/wl), length.out = nrow(m))/1000}
    else {
        fl1 <- flim[1] * nrow(m) * 2000 / f
        fl2 <- flim[2] * nrow(m) * 2000 / f
        m <- m[(fl1:fl2) + 1, ]
        freq <- seq(flim[1], flim[2], length.out = nrow(m))
    }

    ## TIME AXIS
    time <- seq(0, n/f, length.out = ncol(m))

    ## CONTOURS
    t.cont <- apply(m, MARGIN = 2, FUN = aggregate)
    f.cont <- apply(m, MARGIN = 1, FUN = aggregate)
    t.cont <- t.cont/sum(t.cont)
    f.cont <- f.cont/sum(f.cont)
    t.cont.cum <- cumsum(t.cont)
    f.cont.cum <- cumsum(f.cont)

    ## STATISTICS
    P <- fraction/100
    proportions <- as.matrix(c((1 - P)/2, 0.5, P + (1 - P)/2))
    t.quantiles <- apply(proportions, MARGIN = 1, function(x) time[length(t.cont.cum[t.cont.cum <= 
        x]) + 1])
    f.quantiles <- apply(proportions, MARGIN = 1, function(x) freq[length(f.cont.cum[f.cont.cum <= 
        x]) + 1])
    time.P1 <- t.quantiles[1]
    time.M <- t.quantiles[2]
    time.P2 <- t.quantiles[3]
    time.IPR <- time.P2 - time.P1
    freq.P1 <- f.quantiles[1]
    freq.M <- f.quantiles[2]
    freq.P2 <- f.quantiles[3]
    freq.IPR <- freq.P2 - freq.P1
    results <- list(time.contour = cbind(time = time, contour = t.cont), 
        freq.contour = cbind(frequency = freq, contour = f.cont), 
        time.P1 = time.P1, time.M = time.M, time.P2 = time.P2, 
        time.IPR = time.IPR, freq.P1 = freq.P1, freq.M = freq.M, 
        freq.P2 = freq.P2, freq.IPR = freq.IPR)

    ## PLOT AND RETURN
    if (plot) {
        def.par <- par(no.readonly = TRUE)
        on.exit(par(def.par))
        par(mfrow = c(2, 1))
        plot(x = time, y = t.cont, type = type, xlab = "Time (s)", 
            ylab = "Probability", main = "Time envelope", ...)
        segments(x0 = t.quantiles[1], y0 = 0, y1 = t.cont[time == 
            time.P1], col = 2)
        segments(x0 = t.quantiles[2], y0 = 0, y1 = t.cont[time == 
            time.M], col = 2, lwd = 2)
        segments(x0 = t.quantiles[3], y0 = 0, y1 = t.cont[time == 
            time.P2], col = 2)
        plot(x = freq, y = f.cont, type = type, xlab = "Frequency (kHz)", 
            ylab = "Probability", main = "Frequency spectrum", 
            ...)
        segments(x0 = f.quantiles[1], y0 = 0, y1 = f.cont[freq == 
            freq.P1], col = 2)
        segments(x0 = f.quantiles[2], y0 = 0, y1 = f.cont[freq == 
            freq.M], col = 2, lwd = 2)
        segments(x0 = f.quantiles[3], y0 = 0, y1 = f.cont[freq == 
            freq.P2], col = 2)
        invisible(results)
    }
    else return(results)
}



################################################################################
##                                ADDSILW
################################################################################

addsilw<-function(
                  wave,
                  f,
                  at = "end",
                  choose = FALSE,
                  d = NULL,
                  plot = FALSE,
                  output = "matrix",
                  ...)

{
  input<-inputw(wave=wave,f=f); wave<-input$w; f<-input$f; bit <- input$bit; rm(input)

  switch(at,
         start = {at<-0},
         middle = {at<-nrow(wave)/(2*f)},
         end = {at<-nrow(wave)/f}
         )
  
  if(is.null(d)) stop("silence duration has to be set with the argument 'd'")

  if(choose)
    { 
      cat("choose position on the wave\n")
      if(.Platform$OS.type == "windows") flush.console()
      oscillo(wave,f=f)
      coord<-locator(n=1)
      at<-coord$x[1]; abline(v=at,col=2,lty=2)
    }

  pos<-round(at*f)
  wave1<-wave[1:pos,1]
  sil<-rep(0,d*f)
  wave2<-wave[-(1:pos),1]

  wave3 <- c(wave1,sil,wave2)

  wave3 <- outputw(wave=wave3, f=f, format=output)
  
  if(plot)
    {
      oscillo(wave=wave3,f=f,...)
      invisible(wave3)
    }
  else return(wave3)
}

################################################################################
##                                AFILTER                                       
################################################################################

afilter<-function(
                  wave,
                  f,
                  threshold = 5,
                  plot = TRUE,
                  listen = FALSE,                  
                  output = "matrix",
                  ...
                  )

{
  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; bit <- input$bit ; rm(input)

  t1<-max(abs(wave))*(threshold/100)
  wave1 <- ifelse(abs(wave)<=t1,yes=0,no=1)
  wave2 <- wave[,1]*wave1[,1]
  wave2 <- outputw(wave=wave2, f=f, format=output)

  if(plot)
    {
      oscillo(wave=wave2,f=f,...)
      invisible(wave2)
    }
  else
    {
      if(listen) {listen(wave2,f=f)}
      return(wave2)
    }
}


################################################################################
##                                AKAMATSU
################################################################################


akamatsu <- function(Lx,         # L in cm
                     Ly,
                     Lz,
                     mode = c(1,1,1),
                     c = 148000,    # in cm/s
                     plot = FALSE,
                     xlab = "Frequency (kHz)",
                     ylab = "Attenuation distance (cm)",
                     ...)

{          
## if loop to be sure that the largest length is the right one
    if (Lx > Ly){
        a  <- Lx
        Lx <- Ly
        Ly <- a
    }

    Fres <- (c*sqrt((mode[1]/Lx)^2+(mode[2]/Ly)^2+(mode[3]/Lz)^2))/2 # resonance frequency of the aquarium 
    Fcut <- (c*sqrt((1/Ly^2)+(1/Lz^2)))/2 # frequency of the high pass filter 
    F <- list(res=Fres, cut=Fcut)

## Plot showing the attenuation distance in function of the frequency 
    if(plot == TRUE){
        f <- seq(0, Fcut, Fcut/100)
        f <- f[-length(f)]
        D <- (2*log(10)*c)/(4*pi*Fcut*sqrt(1-(f/Fcut)^2))
        plot(f/1000, D, xlab = xlab, ylab = ylab,...)
        invisible(F)
    }

    return(F)
}


################################################################################
##                                AMA
################################################################################

ama<-function(
              wave,
              f,
              envt = "hil",
              wl = 512,
              plot = TRUE,
              type = "l",
              ...)

{
  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; rm(input)
  enve<-env(wave=wave, f=f, envt = envt, plot = FALSE)
  meanspec(wave=enve, f=f, wl=wl, plot=plot, type=type,...)
}


################################################################################
##                                AR                                        
################################################################################

AR <- function(...,
               datatype = "objects",    # if "file" -> only .wav and .mp3 read by tuneR functions redWave and readMP3
               envt = "hil",
               msmooth = NULL,
               ksmooth = NULL,
               ssmooth = NULL,
               pattern = "[wav]$|[mp3]$")
    
{
    ## data as R objects
    if(datatype=="objects"){
        if(is.character(list(...)[[1]])) stop("The first argument should not be a character string with 'datatype=object'")
        ## number of objects = number of arguments - missing arguments
        n <- nargs() - !missing(datatype) - !missing(envt) - !missing(pattern) 
        ## number of objects = number of arguments - NULL arguments
        if(is.null(msmooth)) m <- 0 else m <- 1
        if(is.null(ksmooth)) k <- 0 else k <- 1
        if(is.null(ssmooth)) s <- 0 else s <- 1
        n <- n - m - k- s
        if(n <= 1) warning("Computing the AR index on a single object does not make a lot of sense as AR is ranked index")
        name <- M <- Ht <- AR <- numeric(n)
        for(i in 1:n){ # does not use the function M() to avoid to compute the envelope twice
            wave <- list(...)[[i]]
            if(class(wave)=="Wave"|class(wave)=="WaveMC") {bit <- wave@bit}
            else if(class(wave)=="audioSample") {bit <- wave$bits}
            if(!is.null(names(list(...))[i])) {name[i] <- names(list(...))[i]}
            env <- env(wave, envt=envt, msmooth=msmooth, ksmooth=ksmooth, ssmooth=ssmooth, plot=FALSE)
            env.norm <- env/sum(env)
            M[i] <- median(env)*2^(1-bit)
            Ht[i] <- th(env.norm)
        }
    }
    ## data as .wav or .mp3 files
    ## data should be a path to a directory containing either .wav or .mp3 files
    else if(datatype=="files"){
        if(!is.character(...)) stop("The first argument should be character string when 'datatype='files''")
        files <- dir(..., pattern=pattern)
        n <- length(files)
        name <- M <- Ht <- AR <- numeric(n)
        for(i in 1:n){ # does not use the function M() to avoid to compute the envelope twice
            filename <- basename(files[i])
            nch <- nchar(filename)
            extension <- substr(filename, start = nch-3, stop = nch)
            name[i] <- substr(filename, start=1, stop=nch-4)
            if(extension == ".wav") {wave <- readWave(files[i])}
            else if (extension == ".mp3") {wave <- readMP3(files[i])}
            env <- env(wave, envt=envt, msmooth=msmooth, ksmooth=ksmooth, ssmooth=ssmooth, plot=FALSE)
            env.norm <- env/sum(env)
            M[i] <- median(env)*2^(1-wave@bit)
            Ht[i] <- th(env.norm)
        }
    }
    
    ## OUTPUT
    AR <- (rank(M)*rank(Ht))/(n^2)
    res <- data.frame(M, Ht, AR)
    if(!is.numeric(name)) rownames(res) <- name
    return(res)
}


################################################################################
##                                ATTENUATION                                        
################################################################################

attenuation <- function(lref,
                        dref = 1,
                        dstop,
                        n,
                        plot = TRUE,
                        xlab = "Distance (m)",
                        ylab = "dB",
                        type="l",
                        ...
                        )

{
    d <- seq(from=dref, to=dstop, length.out=n)
    dB <- lref-(20*log10(d/dref))

    if(plot)
        {
            plot(x=d, y=dB, xlab=xlab, ylab=ylab, type=type,...)
            invisible(data)
        }
    else {return(d)}
}


################################################################################
##                                AUTOC                                         
################################################################################

autoc <- function(
                wave,
                f,
                wl = 512,
                fmin = 1,
                fmax = f/2,
                threshold = NULL,
                plot = TRUE,
                xlab = "Time (s)",
                ylab = "Frequency (kHz)",
                ylim = c(0,f/2000),
                pb = FALSE,
                ...)

{
  input <- inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; rm(input)

  if(!is.null(threshold)) wave <- afilter(wave=wave,f=f,threshold=threshold,plot=FALSE)

  lag.min <- round(wl*(fmin/(f/2)))
  lag.max <- round(wl*(fmax/(f/2)))   # in that case we consider 2*fmax to be sure to get a result 
  if(lag.max==wl) {lag.max <- lag.max-1}
  n <- nrow(wave)
  step <- seq(1,n-2*wl,wl)
  N <- length(step) 
  if(pb) {pbar <- txtProgressBar(min=0, max=N, style = 3)}
  R <- matrix(data=numeric((lag.max+1)*N), lag.max+1, N)
  R[is.nan(R)] <- 0  # the use of 'threshold' can produce NaN
  for (i in step)
    {
      R[,which(step==i)] <- as.vector(acf(wave[i:(wl+i-1)], lag.max=lag.max, plot=FALSE)$acf)
    }
  R[is.nan(R)] <- 0 # need to do it twice as acf() can produce NaN
  tfond<-numeric(N)
  excl <- 1:lag.min

  for (i in 1:N)
    {
    tfond[i] <- fpeaks(R[-excl,i], f=NA, nmax=1)[1]
    if(pb) {setTxtProgressBar(pbar, i)}
    }
#  for (k in 1:N) {tfond[k]<-which.max(R[-excl, k])}
#  tfond <- ifelse(tfond==1 | tfond==nrow(R), yes=NA, no=tfond)
  tfond <- tfond + length(excl)
  y<-f/tfond/1000

  x<-seq(0,n/f,length.out=N)  

  if(plot)
    {
      plot(x=x, y=y,
           xlab = xlab,
           ylab = ylab, ylim = ylim,
           las = 1,
           ...)
      invisible(cbind(x,y))
    }

  else return(cbind(x,y))
  if(pb) close(pbar)
}

################################################################################
##                                CCOH                                        
################################################################################

ccoh<-function(
               wave1,
               wave2,
               f,
               wl = 512,
               ovlp = 0,
               plot = TRUE,
               grid = TRUE,
               scale = TRUE,
               cont = FALSE,
               collevels = seq(0,1,0.01),
               palette = reverse.heat.colors,
               contlevels = seq (0,1,0.01),
               colcont = "black",
               colbg = "white",
               colgrid = "black",
               colaxis = "black",
               collab = "black",
               xlab = "Time (s)",
               ylab = "Frequency (kHz)",
               scalelab = "Coherence",
               main = NULL, 
               scalefontlab = 1,
               scalecexlab =0.75,
               axisX = TRUE,
               axisY = TRUE,
               flim = NULL,
               flimd = NULL,
               ...)

{
  input1<-inputw(wave=wave1,f=f) ; wave1<-input1$w ; f<-input1$f ; rm(input1)
  wave2<-inputw(wave=wave2,f=f)$w

  n1<-nrow(wave1)
  n2<-nrow(wave2)
  if(n1 != n2) stop("'wave 1' and 'wave 2' must have the same length")
  n<-n1


                                        # dynamic vertical zoom (modifications of analysis parameters)
  if(!is.null(flimd))
    {
                                        # zoom magnification
      mag<-round((f/2000)/(flimd[2]-flimd[1]))
                                        # new parameters
      wl<-wl*mag
      if(ovlp==0) ovlp<-100
      ovlp<-100-round(ovlp/mag)
                                        # use of normal flim for following axis modifications
      flim<-flimd
    }

  step<-seq(1,n-wl,wl-(ovlp*wl/100))	# coherence windows

  z1<-matrix(data=numeric((wl)*length(step)),wl,length(step))

  for(i in step)
    {
      z1[,which(step==i)]<-spec.pgram(cbind(wave1[i:(wl+i-1),],
                  wave2[i:(wl+i-1),]), spans = c(3,3), fast=FALSE, taper=FALSE, plot=FALSE)$coh
    }

  z<-z1[1:(wl/2),]							  

                                        # X axis settings
  X<-seq(0,n/f,length.out=length(step))

                                        # vertical zoom
  if(is.null(flim))
    {
      Y<-seq(0,f/2000,length.out=nrow(z))
    }
  else
    {
      fl1<-flim[1]*nrow(z)*2000/f
      fl2<-flim[2]*nrow(z)*2000/f
      z<-z[fl1:fl2,]
      Y<-seq(flim[1],flim[2],length.out=nrow(z))
    }
  
  Z<-t(z)
  
  if(plot)
    {
      Zlim<-range(Z, finite = TRUE) 
      
      if(scale)
        {
          def.par <- par(no.readonly = TRUE)
          on.exit(par(def.par))
          layout(matrix(c(1, 2), ncol = 2, byrow=TRUE), widths = c(6, 1))
          par(mar=c(5,4.1,1,0),las=1,cex=1,bg=colbg,col=colaxis,col.axis=colaxis,col.lab=collab)
          filled.contour.modif2(x=X ,y=Y, z=Z, levels=collevels, nlevels=20,
                                plot.title=title(main=main,xlab=xlab,ylab=ylab), color.palette=palette,axisX=axisX, axisY=axisY)
          if(grid) grid(nx=NA, ny=NULL, col=colgrid)
          if(colaxis != colgrid) abline(h=0,col=colaxis) else abline(h=0,col=colgrid)
          par(mar=c(5,1,4.5,3),las=0)
          dBscale(collevels=collevels,palette=palette,fontlab=scalefontlab,
                  cexlab=scalecexlab,collab=collab,textlab=scalelab,colaxis=colaxis)
        }
      
      if(scale==FALSE)
        {
          par(las=1, col=colaxis, col.axis=colaxis, col.lab=collab,,bg=colbg,...)
          filled.contour.modif2(x=X ,y=Y, z=Z, levels=collevels, nlevels=20,
                                plot.title=title(main=main,xlab=xlab,ylab=ylab), color.palette=palette, axisX=axisX, axisY=axisY,
                                col.lab=collab,colaxis=colaxis)		
          if(grid) grid(nx=NA, ny=NULL, col=colgrid)
          if(colaxis != colgrid) abline(h=0,col=colaxis) else abline(h=0,col=colgrid)
        }

      if(cont) 
        {
          contour(X,Y,Z,add=TRUE,
                  levels=contlevels,nlevels=5,col=colcont,...)
        }  
      invisible(list(time=X, freq=Y, coh=Z))
    }
  else return(list(time=X, freq=Y, coh=Z))
}




################################################################################
##                                CEPS                                         
################################################################################

ceps<-function(
               wave,
               f,
               wl = 512,
               at = NULL,
               from = NULL,
               to = NULL,
               tidentify = FALSE,   # identify in seconds
               fidentify = FALSE,   # identify in Hz
               col = "black",
               cex = 1,
               plot = TRUE,
               qlab = "Quefrency (bottom: s, up: Hz)",
               alab = "Amplitude",
               qlim = NULL,
               alim = NULL,
               type= "l",
               ...)

{
  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; rm(input)

  if(!is.null(from)|!is.null(to))
    {
      if(is.null(from) && !is.null(to)) {a<-1; b<-round(to*f)}
      if(!is.null(from) && is.null(to)) {a<-round(from*f); b<-length(wave)}
      if(!is.null(from) && !is.null(to))
        {
          if(from>to) stop("'from' cannot be superior to 'to'")
          if(from==0) {a<-1} else a<-round(from*f)
          b<-round(to*f)
        }
      wave<-as.matrix(wave[a:b,])
    }

  if(!is.null(at))
    {
      if(wl==FALSE) stop("Argument 'wl' has to be set up, for instance wl=512")
      c<-round(at*f)
      wl2<-wl%/%2
      wave<-as.matrix(wave[(c-wl2):(c+wl2),])
    }

  n<-nrow(wave)
  N<-round(n/2)

  z1<-Re(fft(log(abs(fft(wave[,1]))),inverse=TRUE))
  z<-z1[1:N]

  x<-seq(0,N/f,length.out=N)

  if(plot)
    {
      plot(x=x,y=z,
           xlab=qlab,xaxt="n",xaxs="i",
           ylab=alab,yaxt="n",yaxs="i",
           type=type,col=col,cex=cex,xlim=qlim,ylim=alim,...)
      if(!is.null(qlim)) E<-qlim[2] else E<-N/f
      X<-seq(0,E,length.out=7)
      axis(side=1,at=X, labels=round(X,3))
      axis(side=3,at=X, labels=round(1/X,3))

      if(tidentify)
        {
          cat("time identification: choose points on the cepstrum\n")
          if(.Platform$OS.type == "windows") flush.console()
          id<-identify(x=x,y=z,labels=round(x,5),tolerance=0.15,col="red")
          return(round(x[id],5))
        }
      
      if(fidentify)
        {
          cat("frequency identification: choose points on the cepstrum\n")
          if(.Platform$OS.type == "windows") flush.console()
          id<-identify(x=x,y=z,labels=round(1/round(x,5),1),tolerance=0.15,col="red")
          return(round(1/round(x[id],5),1))
        }
      results <- cbind(x,z)
      invisible(results)
    }

  else
    {
      results <- cbind(x,z)
      return(results)
    }
}


################################################################################
##                                   CEPSTRO                                    
################################################################################

cepstro<-function(
                  wave,
                  f,
                  wl = 512,
                  ovlp = 0,
                  plot = TRUE,
                  grid = TRUE,
                  scale = TRUE,
                  cont = FALSE,
                  collevels = seq(0,1,0.01),
                  palette = reverse.heat.colors,
                  contlevels = seq (0,1,0.01),
                  colcont = "black",
                  colbg = "white",
                  colgrid = "black",
                  colaxis = "black",
                  collab = "black",
                  xlab = "Time (s)",
                  ylab = "Quefrency (ms)",
                  scalelab = "Amplitude",
                  main = NULL,
                  scalefontlab = 1,
                  scalecexlab =0.75,
                  axisX = TRUE,
                  axisY = TRUE,
                  ...)

{
  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; rm(input)
  wave<-ifelse(wave==0,yes=1e-6,no=wave)

  n<-nrow(wave)
  p<-round(n/2)
  step<-seq(1,n-wl,wl-(ovlp*wl/100))

  N<-length(step)
  WL<-wl%/%2
  z1<-matrix(data=numeric(wl*N),wl,N)
  for(i in step)
    {z1[,which(step==i)]<-Re(fft(log(abs(fft(wave[i:(wl+i-1),]))),inverse=TRUE))}
  z2<-z1[1:WL,]
  z<-ifelse(z2=="NaN"|z2=="-Inf"|z2<=0,yes=0,no=z2)
  Z<-t(z/max(z))

  if(plot)
    {
      X<-seq(0,n/f,length.out=length(step))
      Y<-seq(0,WL/f,length.out=nrow(z))*1000
      Zlim<-range(Z, finite = TRUE)
      if(scale)
        {
          def.par <- par(no.readonly = TRUE)
          on.exit(par(def.par))
          layout(matrix(c(1, 2), ncol = 2, byrow=TRUE), widths = c(6, 1))
          par(mar=c(5,4.1,1,0),las=1,cex=1,col=colaxis,col.axis=colaxis,col.lab=collab,bg=colbg)
          filled.contour.modif2(x=X ,y=Y, z=Z, levels=collevels, nlevels=20,
                                plot.title=title(main=main,xlab=xlab,ylab=ylab),
                                color.palette=palette,axisX=axisX, axisY=axisY)
          if(grid) grid(nx=NA, ny=NULL, col=colgrid)
          if(colaxis != colgrid) abline(h=0,col=colaxis) else abline(h=0,col=colgrid)
          par(mar=c(5,1,4.5,3),las=0)
          dBscale(collevels=collevels,palette=palette,fontlab=scalefontlab,
                  cexlab=scalecexlab,collab=collab,textlab=scalelab,colaxis=colaxis)
        }

      if(scale==FALSE)
        {
          par(las=1, col=colaxis, col.axis=colaxis, col.lab=collab,bg=colbg,...)
          filled.contour.modif2(x=X ,y=Y, z=Z, levels=collevels, nlevels=20,
                                plot.title=title(main=main,xlab=xlab,ylab=ylab),
                                color.palette=palette, axisX=axisX, axisY=axisY,
                                col.lab=collab,colaxis=colaxis)
          if(grid) grid(nx=NA, ny=NULL, col=colgrid)
          if(colaxis != colgrid) abline(h=0,col=colaxis) else abline(h=0,col=colgrid)
        }

      if(cont)
        {
          contour(X,Y,Z,add=TRUE,
                  levels=contlevels,nlevels=5,col=colcont,...)
        }
      invisible(list(time=X, quef=Y, amp=Z))
    }

  else{return(list(time=X, quef=Y, amp=Z))}
}



################################################################################
##                                COH                                         
###############################################################################

coh<-function(
              wave1,
              wave2,
              f,
              plot =TRUE,
              xlab = "Frequency (kHz)",
              ylab = "Coherence",
              xlim = c(0,f/2000),
              type = "l",
              ...
              )

{
  input1<-inputw(wave=wave1,f=f) ; wave1<-input1$w ; f<-input1$f ; rm(input1)
  wave2<-inputw(wave=wave2,f=f)$w

  n1<-nrow(wave1)
  n2<-nrow(wave2)
  if(n1 != n2) stop("'wave 1' and 'wave 2' must have the same length")

  Y<-spec.pgram(cbind(wave1, wave2), fast=FALSE, taper=FALSE,
                spans = c(3,3),plot=FALSE)$coh
  X<-seq(0,f/2000,length.out=nrow(Y))

  if(plot)
    {
      plot(x=X,y=Y,xlab=xlab,ylab=ylab,xlim=xlim,type=type,...)
      invisible(cbind(X,Y))
    }
  else return(cbind(X,Y))
}


################################################################################
##                                COMBFILTER                                         
###############################################################################

combfilter <- function(wave,
                 f,
                 alpha,
                 K,
                 units = c("samples", "seconds"),
                 plot = FALSE,
                 output="matrix",
                 ...)
{
    ## INPUT
    input <- inputw(wave=wave,f=f) ; wave <- input$w ; f <- input$f ; rm(input)
    wave <- wave/max(wave)
    n <- length(wave)
    K <- switch(match.arg(units), samples = K, seconds = round(K*f))

    ## FILTER
    wave1 <- c(rep(0,K), wave[1:(n-K)]) + alpha*wave
    wave2 <- outputw(wave=wave1, f=f, format=output)

    ## FREQUENCY RESPONSE
    w <- seq(0, pi, length.out=n)
    H <- sqrt(1+alpha^2 + 2*alpha*cos(w*K))
    
    ## PLOT AND RETURN
    if(plot){
       def.par <- par(no.readonly = TRUE)
       on.exit(par(def.par))
       layout(matrix(c(1, 2), ncol = 2, byrow=TRUE), widths = c(4, 1))
        spectro(wave2, f=f, scale=FALSE, osc=FALSE, ...)
        par(mar=c(5.1,0,4.1,1))
        freq <- seq(0, f/2000, length.out=n)
        plot(x=H, y=freq, xlim=c(0, max(H)), type="l", xaxs="i", yaxs="i", xlab="Linear amplitude", ylab="", yaxt="n")
        invisible(wave2)
    }
    else return(wave2)
}


################################################################################
##                                CONVSPL                                         
###############################################################################

convSPL<-function(
                  x,
                  d = 1,
                  Iref = 10^-12,
                  pref = 2*10^-5 
                  )

{
  P<-4*pi*(d^2)*Iref*(10^(x/10))
  I<-Iref*(10^(x/10))
  p<-pref*(10^(x/20))
  conv<-list(P = P, I = I, p = p)
  return(conv)
}


################################################################################
##                                CORENV                                         
################################################################################

corenv<-function(
                 wave1,
                 wave2,
                 f,
                 envt="hil",
                 msmooth = NULL,
                 ksmooth = NULL,
                 plot = TRUE,
                 plotval = TRUE,
                 method = "spearman",
                 col = "black",
                 colval = "red",
                 cexval = 1,
                 fontval = 1,
                 xlab = "Time (s)",
                 ylab = "Coefficient of correlation (r)",
                 type= "l",
                 pb = FALSE,
                 ...)

{
  input1<-inputw(wave=wave1,f=f) ; wave1<-input1$w ; f<-input1$f ; rm(input1)
  wave2<-inputw(wave=wave2,f=f)$w

  n<-nrow(wave1)

  cat("please wait...\n")
  if(.Platform$OS.type == "windows") flush.console()
  
  x<-env(wave=wave1,f=f,envt=envt,msmooth=msmooth,ksmooth=ksmooth,plot=FALSE)
  y<-env(wave=wave2,f=f,envt=envt,msmooth=msmooth,ksmooth=ksmooth,plot=FALSE)

  nx<-nrow(x)
  ny<-nrow(y)

  if(nx != ny) stop("'wave 1' and 'wave 2' must have the same length")

  meanx<-mean(x)
  meany<-mean(y)
  diffx<-x-meanx
  r1<-numeric(nx)
  r2<-numeric(nx)

  if(pb) {pbar <- txtProgressBar(min=0, max=nx-1, style = 3)}
  for (i in 0:(nx-1))
    {
    r1[i+1]<-cor(x=x,y=c(y[(i+1):ny],rep(0,i)),method = method)
    r2[i+1]<-cor(x=x,y=c(rep(0,i),y[1:(ny-i)]),method = method)
    if(pb) {setTxtProgressBar(pbar, i)}
    }

  r2<-r2[-1]
  r<-c(rev(r1),r2)
  rmax<-max(r,na.rm=TRUE)
  offset<-which.max(r)
  if(offset<=(length(r)/2)) {offsetp<-which.max(r)} else {offsetp<-which.max(r)-1}
  if(!is.null(msmooth)|!is.null(ksmooth)) {F<-f*nx/n; offsett<-(offsetp-nx)/F}
  else{offsett<-(offsetp-n)/f}

  if(offsetp < nx){p<-cor.test(x=x, y=c(y[(nx-offsetp+1):ny],rep(0,nx-offsetp)),method = method)}
  else {p<-cor.test(x=x, y=c(rep(0,offsetp-nx),y[1:(ny-(offsetp-nx))]), method = method)}
  p<-p$p.value

  X<-seq(-n/f,n/f,length.out=2*nx-1)
  corr<-list(r = cbind(X,r), rmax = rmax, p = p, t = offsett)
  
  if(plot)
    {
      plot(x = X, y = r, xlab = xlab, ylab = ylab, col = col, type = type,...)
      if(plotval)
        {
          mtext(paste(
                      "rmax = ", as.character(round(rmax,2)),
                      ", offset = ", as.character(round(offsett,3)), "s", sep=" "),
                side=3, line=-2, col=colval, cex=cexval, font=fontval)
          segments(x0=round(offsett,3),y0=min(r), x1=round(offsett,3), y1=rmax, col=colval, lty=2)
          segments(x0=-n/f,y0=rmax, x1=round(offsett,3), y1=rmax, col=colval, lty=2)
          points(x=round(offsett,3),y=rmax,pch=19,cex=1, col=colval)
        }
      invisible(corr)
    }

  else
    {
      return(corr)
    }
 if(pb) close(pbar)
}



################################################################################
##                                CORSPEC                                         
################################################################################

corspec <- function(
                    spec1,
                    spec2,
                    f = NULL,
                    mel = FALSE,
                    plot = TRUE,
                    plotval = TRUE,
                    method = "spearman",
                    col = "black",
                    colval = "red",
                    cexval = 1,
                    fontval = 1,
                    xlab = NULL,
                    ylab = "Coefficient of correlation (r)",
                    type ="l",
                    ...
                    )

{
    options(warn=-1) ## to remove warning messages from cor.test()
    if(is.matrix(spec1) && ncol(spec1)==2) s1<-spec1[,2] else s1 <- spec1
    if(is.matrix(spec2) && ncol(spec2)==2) s2<-spec2[,2] else s2 <- spec2

    n1<-length(s1)
    n2<-length(s2)

    if(n1 != n2) stop("'spec1' and 'spec2' must have the same length")
    if(any(s1 < 0) | any(s2 < 0)) stop("data does not have to be in dB")
    if(sum(s1-s2)==0) stop("'spec1' and 'spec2' are the same object")
    
    mean1<-mean(s1)
    mean2<-mean(s2)
    diffx<-s1-mean1
    r1<-numeric(n1)
    r2<-numeric(n2)

    for (i in 0:(n1-1))
        {
            r1[i]<-cor(x=s1,y=c(s2[(i+1):n2],rep(0,i)),method = method)
        }

    for (i in 0:(n1-1))
        {
            r2[i+1]<-cor(x=s1,y=c(rep(0,i),s2[1:(n2-i)]),method = method)
        }

    r2<-r2[-1]
    r<-c(rev(r1),r2)

    rmax<-max(r,na.rm=TRUE)
    offset<-which.max(r)
    if(offset<=(length(r)/2)) {offsetp<-which.max(r)} else {offsetp<-which.max(r)-1}

    if(offsetp < n1)
        {
            p<-cor.test(x=s1, y=c(s2[(n1-offsetp+1):n2],rep(0,n1-offsetp)),
                        method = method)
        }
    else
        {
            p<-cor.test(x=s1, y=c(rep(0,offsetp-n1),s2[1:(n2-(offsetp-n1))]),
                        method = method)
        }
    p<-p$p.value
    

    if (!is.null(f) & mel) {f <- 2*mel(f/2)}
    if(is.null(f))
        {
            if(is.vector(spec1) & is.vector(spec2))
                {
                    if(plot) warning("'f' is missing, frequency of maximum correlation cannot be computed - No plot produced", call.=FALSE)
                    else warning("'f' is missing, frequency of maximum correlation cannot be computed", call.=FALSE)
                    corr <- list(r = r, rmax = rmax, p = p, f = NA)
                    return(corr)
                }
            else
                {
                    if(is.matrix(spec1)) f<-spec1[nrow(spec1),1]*2000*nrow(spec1)/(nrow(spec1)-1)
                    else if(is.matrix(spec2)) f<-spec2[nrow(spec2),1]*2000*nrow(spec2)/(nrow(spec2)-1)
                }
        }

    range<-c(0,f/2000)
    offsetf<-((range[2]-range[1])*(offsetp-n1))/n1
    X<-seq(-range[2],range[2],length.out=2*n1-1)

    corr<-list(r = cbind(X,r), rmax = rmax, p = p, f = offsetf)

    if(plot)
        {if (mel) {w <- " kmel"}
         else {w <- " kHz"}
         if(is.null(xlab)) {if (mel) xlab <- "Frequency (kmel)" else xlab <- "Frequency (kHz)"}
         plot(x = X, y = r, xlab = xlab, ylab = ylab, col = col, type = type,...)
         if(plotval)
             {
                 mtext(paste(
                     "rmax = ", as.character(round(rmax,2)),
                     ", offset = ", as.character(round(offsetf,2)), w, sep=" "),
                       side=3, line=-2, col=colval, cex=cexval, font=fontval)
                 segments(x0=offsetf,y0=min(r), x1=offsetf, y1=rmax, col=colval, lty=2)
                 segments(x0=-range[2],y0=rmax, x1=offsetf, y1=rmax, col=colval, lty=2)
                 points(x=offsetf,y=rmax,pch=19,cex=1, col=colval)
             }
         invisible(corr)
     }
    
    else  
        {
            return(corr)
        }
    options(warn=0) ## to reset warning messages
}


################################################################################
##                                COVSPECTRO                                         
################################################################################


covspectro<-function(
                     wave1,
                     wave2,
                     f,
                     wl = 512,
                     wn = "hanning",
                     n,
                     plot = TRUE,
                     plotval = TRUE,
                     method = "spearman",
                     col = "black",
                     colval = "red",
                     cexval = 1,
                     fontval = 1,
                     xlab = "Time (s)",
                     ylab = "Normalised covariance (cov)",
                     type ="l",
                     pb = FALSE,
                     ...
                     )

{
  input1<-inputw(wave=wave1,f=f) ; wave1<-input1$w ; f<-input1$f ; rm(input1)
  wave2<-inputw(wave=wave2,f=f)$w

  if(n>21)
    {
      cat("please wait...")
      if(.Platform$OS.type == "windows") flush.console()
    }
  
  if(n %% 2 != 1) stop("'n' must be odd")

  n<-(n%/%2)+1

  n1<-nrow(wave1)
  n2<-nrow(wave2)

  if(n1 != n2) stop("'wave 1' and 'wave 2' must have the same length")

  step1<-seq(1,n1-wl,wl); lstep1<-length(step1)
  step2<-round(seq(1,n2,length.out=n))

                                        # wave not time shifted
  spectro1<-sspectro(wave=wave1,f=f,wl=wl,wn=wn)
  WL<-wl%/%2
  spectro2a<-array(numeric(WL*lstep1*n),dim=c(WL,lstep1,n))
  spectro2b<-array(numeric(WL*lstep1*n),dim=c(WL,lstep1,n))
  cov1<-numeric(n)
  cov2<-numeric(n)

                                        # wave time shifted
                                        # successive spectrograms
                                        # covariance of spectrogram1/spectra1 with spectrogram2/spectra1,
                                        # spectrogram1/spectra2 with spectrogram2/spectra2 and so on
                                        # diagonal of the cov matrix and mean of this diagonal
                                        # one mean cov for comparaison between 2 spectrograms for(i in step2)
  if(pb) {pbar <- txtProgressBar(min=0, max=n2, style = 3)}
  for (i in step2)
    {
      spectro2a[,,which(step2==i)]<-sspectro(wave=as.matrix(c(wave2[i:n2],rep(0,i-1))),f=f,wl=wl,wn=wn)
      spectro2a<-ifelse(spectro2a=="NaN",yes=0,no=spectro2a)
      cov1[which(step2==i)]<-mean(diag(cov(x=spectro1,y=spectro2a[,,which(step2==i)],method = method)))
      spectro2b[,,which(step2==i)]<-sspectro(wave=as.matrix(c(rep(0,i),wave2[1:(n2-i)])),f=f,wl=wl,wn=wn)
      spectro2b<-ifelse(spectro2b=="NaN",yes=0,no=spectro2b)
      cov2[which(step2==i)]<-mean(diag(cov(x=spectro1,y=spectro2b[,,which(step2==i)],method = method)))
      if(pb) {setTxtProgressBar(pbar, i)}
    }

                                        # to normalise the covariance we need covmax that is the autocovariance of spectro1
  covmax<-mean(diag(cov(x=spectro1,y=spectro1,method = method)))

                                        # discard the first value of cov2 that is already computed in cov1
  cov2<-cov2[-1]
  cov3<-c(rev(cov1),cov2)
  cov4<-cov3/covmax
  cov4max<-max(cov4)
  offset<-which.max(cov4)
  offsetp<-which.max(cov4)
  offsett<-(((offsetp*n1)/n)-n1)/f
  covar<-list(cov = cov4, covmax = cov4max, t = offsett)

  if(plot)
    {
      x<-seq(-n1/f,n1/f,length.out=(2*n)-1)
      plot(x = x, y = cov4, xlab = xlab, ylab = ylab, col = col, type = type,...)
      if(plotval)
        {
          mtext(paste(
                      "covmax = ", as.character(round(cov4max,2)),
                      ", offset = ", as.character(round(offsett,3)), "s", sep=" "),
                side=3, line=-2, col=colval, cex=cexval, font=fontval)
          segments(x0=round(offsett,3),y0=min(cov4), x1=round(offsett,3), y1=cov4max, col=colval, lty=2)
          segments(x0=-n1/f,y0=cov4max, x1=round(offsett,3), y1=cov4max, col=colval, lty=2)
          points(x=round(offsett,3),y=cov4max,pch=19,cex=1, col=colval)
        }
      invisible(covar)
    }

  else
    {
      return(covar)
    }
  if(pb) close(pbar)

}


################################################################################
##                               CREST
################################################################################

crest <- function(
                  wave,
                  f,
                  plot = FALSE,
                  col = 2,
                  cex = 3,
                  symbol = "*",
                  ...
                  )

{  
  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; rm(input)
  max <- max(abs(wave))
  loc.max <- which(wave == max(wave))/f
  c <- max/rms(wave)

  if(plot)
    {
      if(which.max(abs(range(wave)))==1) {max2 <- -max} else {max2 <- max}# when the max is actually a negative value
      oscillo(wave=wave, f=f,...)
      text(x=loc.max, y=max2, labels=symbol, cex=cex, col=col)
    }
  return(list(C=c, val=max, loc=loc.max))
}


################################################################################
##                                CSH
################################################################################


csh<-function(
              wave,
              f,
              wl = 512,
              wn = "hanning",
              ovlp = 0,
              fftw = FALSE,
              threshold = NULL,
              plot = TRUE,
              xlab = "Times (s)",
              ylab = "Spectral Entropy",
              ylim = c(0,1.1),
              type = "l",
              ...)

{
  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; rm(input)

                                        # threshold
  if(!is.null(threshold)) wave<-afilter(wave=wave,f=f,threshold=threshold,plot=FALSE)

                                        # STFT (see function spectro())
  n<-nrow(wave)
  step<-seq(1,n-wl,wl-(ovlp*wl/100))
  z <- stft(wave=wave, f=f, wl=wl, zp=0, step=step, wn=wn, fftw=fftw)
                                        # sh applied to the Fourier matrix
  h<-apply(z, MARGIN=2, FUN=sh)

  t<-seq(0,n/f,length.out=length(step))

                                        # graphic
  if(plot)
    {
      plot(x=t, y=h,
           xaxs = "i", xlab = xlab,
           yaxs = "i", ylab = ylab, ylim = ylim,
           type = type,
           ...)
      invisible(cbind(t,h))
    }
  else
    return(cbind(t,h))
}


################################################################################
##                                CUTSPEC                                       
################################################################################

cutspec <- function(spec,
                     f = NULL,
                     flim,
                     mel = FALSE,
                     norm = FALSE,
                     PMF = FALSE
                     )

{
    if(norm & PMF) stop("'norm' and 'PMF' should not be both set to TRUE")

    if(is.vector(spec))
        {
            if(is.null(f)) stop("'f' is missing and is necessary when 'spec' is a vector")
            if(mel) {
                f <- 2*mel(f/2)
		flim <- mel(flim*1000)/1000
            } 
            wl <- length(spec)*2
            if(flim==f/2000) {i <- 1; j=0} else {i <- j <- 1}           ## to avoid to search over Nyquist - 1
            specut <- spec[(flim[1]*1000*wl/f+i):(flim[2]*1000*wl/f+j)] ## JS added +1 because first value at '0 Hz'
            if(norm) {specut <- specut/max(specut)}
            if(PMF)  {specut <- specut/sum(specut)}     
        }

    else if(is.matrix(spec))
        {
            if(ncol(spec)>2) {stop("'spec' should not have more than two columns")}
            if(is.null(f)) {
                f <- spec[nrow(spec),1]*2000*nrow(spec)/(nrow(spec)-1)
                if (mel) {flim <- mel(flim*1000)}
            }
            wl <- nrow(spec)*2
            if(flim[2]==f/2000) {i <- 1; j=0} else {i <- j <- 1}           ## to avoid to search over Nyquist - 1
            specut <- spec[(flim[1]*1000*wl/f+i):(flim[2]*1000*wl/f+j), , drop=FALSE] ## JS added +1 because first value at '0 Hz'
            if(norm) {specut[,2] <- specut[,2]/max(specut[,2])}
            if(PMF)  {specut[,2] <- specut[,2]/sum(specut[,2])}
        }

    return(specut)
}


################################################################################
##                                CUTW                                         
################################################################################

cutw<-function(
               wave,
               f,
               from = NULL,
               to = NULL,
               choose = FALSE,
               plot = FALSE,
               marks = TRUE,               
               output = "matrix",
               ...)

{
  input <- inputw(wave=wave,f=f) ; wave <- input$w ; f <- input$f ; bit <- input$bit ; rm(input)

  if(choose)
    { 
      cat("choose start and end positions on the wave\n")
      if(.Platform$OS.type == "windows") flush.console()
      oscillo(wave,f=f)
      coord<-locator(n=2)
      from<-coord$x[1]; a<-round(from*f) ; abline(v=from,col=2,lty=2)
      to<-coord$x[2]; b<-round(to*f); abline(v=to,col=2,lty=2)
    }
  else if(!is.null(from)|!is.null(to))
    {
      if(is.null(from) && !is.null(to)) {a<-1; b<-round(to*f)}
      if(!is.null(from) && is.null(to)) {a<-round(from*f); b<-length(wave)}
      if(!is.null(from) && !is.null(to))
        {
          if(from>to) stop("'from' cannot be superior to 'to'")
          if(from==0) {a<-1} else a<-round(from*f)
          b<-round(to*f)
        }
    }

  wavecut <- as.matrix(wave[a:b,])
#  wavecut<-wavecut1/max(abs(wavecut1))

  wavecut <- outputw(wave=wavecut, f=f, bit=bit, format=output)

  if(plot)
    {
      def.par <- par(no.readonly = TRUE)
      on.exit(par(def.par))
      par(mfrow=c(2,1),oma=c(0,0.1,0,0))
      oscillo(wave,f=f,...)
      title(main="original")
      if(marks)
        {
          abline(v=from, col="red", lty=2)
          abline(v=to, col="red", lty=2)
        }
      oscillo(wavecut,f=f,...)
      title(main="cut")
      invisible(wavecut)
    }
  else return(wavecut)
}


################################################################################
##                                DBSCALE                                        
################################################################################

dBscale<-function
(
 collevels,
 palette = spectro.colors,
 side = 4,
 textlab = "Amplitude\n(dB)",
 cexlab = 0.75,
 fontlab = 1,
 collab = "black",
 colaxis = "black",
 ...
 )

{
  plot.new()
  levels<-collevels
  col <- palette(length(collevels) - 1)
  par(las=1)
  
  if(side == 2 | side == 4)
    {    
      plot.window(xlim = c(0, 1), ylim = range(collevels), xaxs = "i",
                  yaxs = "i")
      mtext(textlab, side=3, outer=FALSE, line=1.5, adj=0, font=fontlab, cex=cexlab, col=collab)
      rect(xleft=0, ybottom=levels[-length(levels)], xright=0.95, ytop=levels[-1],
           col = col, lty=0, border = TRUE)
      segments(x0=0,y0=max(collevels),x1=0.95,y1=max(collevels),col=colaxis)
      segments(x0=0,y0=min(collevels),x1=0.95,y1=min(collevels),col=colaxis)          
      abline(v=c(0,0.95),col=colaxis)
      if(side == 2) axis(2,col=colaxis,col.axis=colaxis,...)
      if(side == 4) axis(4,pos=0.95,col=colaxis,col.axis=colaxis,...)
    }

  if(side == 1  | side == 3)
    {    
      plot.window(xlim = range(collevels), ylim = c(0, 1), xaxs = "i",
                  yaxs = "i")
      mtext(textlab, side=3, outer=FALSE, line=1.5, adj=0, font=fontlab, cex=cexlab, col=collab)
      rect(xleft=levels[-length(levels)], ybottom=0, xright=levels[-1], ytop=0.95, col = col, lty=0)
      segments(x0=min(collevels),y0=0,x1=min(collevels),y1=0.95,col=colaxis)
      segments(x0=max(collevels),y0=0,x1=max(collevels),y1=0.95,col=colaxis)       
      abline(h=c(0,0.95),col=colaxis)
      if(side == 1) axis(1,col=colaxis,col.axis=colaxis,...)
      if(side == 3) axis(3,pos=0.95,col=colaxis,col.axis=colaxis,...)
    }    
}    


################################################################################
##                                DBWEIGHT
################################################################################

dBweight <- function(
                    f,
                    dBref=NULL
                    )
  
{
                                        # dBA
  num <- (12200^2*f^4)
  den <- (f^2 + 20.6^2)* sqrt((f^2 + 107.7^2)*(f^2 + 737.9^2))*(f^2 + 12200^2)
  A <- 2 + 20*log10(num/den)
  A <- ifelse(is.infinite(A), yes=NA, no=A)

                                        # dBB
  num <- (12200^2*f^3);
  den <- (f^2 + 20.6^2)*sqrt((f^2 + 158.5^2))*(f^2 + 12200^2)
  B <- 0.17 + 20*log10(num/den)
  B <- ifelse(is.infinite(B), yes=NA, no=B)

                                        #dBC
  num <- (12200^2*f^2);
  den <- (f^2 + 20.6^2)*(f^2 + 12200^2)
  C <- 0.06 + 20*log10(num/den)
  C <- ifelse(is.infinite(C), yes=NA, no=C)

                                        #dBD
  a <- f/(6.8966888496476*10^-5);
  h <- ((1037918.48 - f^2)^2 + (1080768.16*f^2))/((9837328 - f^2)^2 + (11723776*f^2))
  b <- sqrt(h/((f^2 + 79919.29)*(f^2 + 1345600)))
  D <- 20*log10(a*b)
  D <- ifelse(is.infinite(D), yes=NA, no=D)


                                        # result
  result <- list(A=A, B=B, C=C, D=D)
  if(!is.null(dBref)) result <- list(A=dBref+A, B=dBref+B, C=dBref+C, D=dBref+D)
  return(result)
}


################################################################################
##                                DELETEW
################################################################################

deletew<-function(
                  wave,
                  f,
                  from = NULL,
                  to = NULL,
                  choose = FALSE,
                  plot = FALSE,
                  marks = TRUE,
                  output = "matrix",
                  ...)

{
  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; bit <- input$bit ; rm(input)

  if(choose)
    { 
      cat("choose start and end positions on the wave\n")
      if(.Platform$OS.type == "windows") flush.console()
      oscillo(wave,f=f)
      coord<-locator(n=2)
      from<-coord$x[1]; a<-round(from*f) ; abline(v=from,col=2,lty=2)
      to<-coord$x[2]; b<-round(to*f); abline(v=to,col=2,lty=2)
    }
  else if(!is.null(from)|!is.null(to))
    {
      if(is.null(from) && !is.null(to)) {a<-1; b<-round(to*f)}
      if(!is.null(from) && is.null(to)) {a<-round(from*f); b<-length(wave)}
      if(!is.null(from) && !is.null(to))
        {
          if(from>to) stop("'from' cannot be superior to 'to'")
          if(from==0) {a<-1} else a<-round(from*f)
          b<-round(to*f)
        }
    }
  wavecut <- wave[-(a:b),]

  wavecut <- outputw(wave=wavecut, f=f, format=output)

  if(plot)
    {
      def.par <- par(no.readonly = TRUE)
      on.exit(par(def.par))
      par(mfrow=c(2,1))
      oscillo(wave,f=f,k=1,j=1,...)
      title(main="original")
      if(marks)
        {
          abline(v=from, col="red", lty=2)
          abline(v=to, col="red", lty=2)
        }
      oscillo(wavecut,f=f,k=1,j=1,...)
      title(main="after deletion")
      invisible(wavecut)
    }
  else return(wavecut)
}


################################################################################
##                                DFREQ                                         
################################################################################

dfreq <- function(
                wave,
                f,
                wl = 512,
                wn = "hanning",
                ovlp = 0,
                fftw = FALSE,
                at = NULL,
                tlim = NULL,
                threshold = NULL,
                bandpass = NULL,
                clip = NULL,  
                plot = TRUE,
                xlab = "Times (s)",
                ylab = "Frequency (kHz)",
                ylim = c(0,f/2000),
                ...)

{
  # error messages
  if(!is.null(at) && ovlp != 0) stop("The 'ovlp' argument cannot bue used in conjunction with the arguement 'at'.")
  if(!is.null(clip)) {if(clip <=0 | clip >= 1) stop("'clip' value has to be superior to 0 and inferior to 1")}

  # input
  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; rm(input)

  # time limits
  if(!is.null(tlim)) wave<-cutw(wave,f=f,from=tlim[1],to=tlim[2])                                      
  
  # amplitude threshold
  if(!is.null(threshold)) {wave<-afilter(wave=wave,f=f,threshold=threshold,plot=FALSE)}

  # Position(s)
  n<-nrow(wave)
  if(!is.null(at))
       {
       step <- at*f
       N <- length(step)
       if(step[1]==0) {step[1] <- 1}
       if(step[N] == n | step[N] == n+1 | step[N] == n-1) {step[N] <- n-wl}
       x <- c(0, at, n/f)
       }
  else {
       step <- seq(1,n-wl,wl-(ovlp*wl/100))
       N <- length(step)
       x <- seq(0, n/f, length.out=N)
       }

  # Fourier
  step <- round(step)
  y1 <- stft(wave=wave, f=f, wl=wl, zp=0, step=step, wn=wn, fftw=fftw) 
  
  # bandpass filter, values outside the bandpass limits are replaced by 0 values
  if(!is.null(bandpass))
    {
      if(length(bandpass)!=2) stop("'The argument 'bandpass' should be a numeric vector of length 2'")
      if(bandpass[1] > bandpass[2]) stop("The first element of 'bandpass' has to be inferior to the second element, i.e. bandpass[1] < bandpass[2]")
      if(bandpass[1] == bandpass[2]) stop("The limits of the bandpass have to be different")
      lowlimit <-round((wl*bandpass[1])/f)
      upperlimit <-round((wl*bandpass[2])/f)
      y1[-(lowlimit:upperlimit),] <- 0
   }

  # Maximum search
  maxi <- apply(y1, MARGIN=2, FUN=max)
  y2 <- apply(y1, MARGIN=2, FUN=which.max)
  y2[which(maxi==0)] <- NA
 
  # discards peaks with an amplitude lower than the clip value
  if(!is.null(clip))
    {
      maxi <- apply(y1, MARGIN=2, FUN=max)
      y2[which(maxi < clip)] <- NA
    }
  # converts into frequency
  y <- (f*y2)/(1000*wl)

  if(!is.null(at)) {y <- c(NA, y, NA)}

  if(plot)
    {
      plot(x=x, y=y,
           xaxs="i", xlab = xlab,
           yaxs="i", ylab = ylab, ylim = ylim,
           ...)
      invisible(cbind(x,y))      
    }
  else
    return(cbind(x,y))
}


################################################################################
##                                DIFFCUMSPEC                                        
################################################################################

diffcumspec <- function(spec1,
                        spec2,
                        f = NULL,
                        mel = FALSE,
                        plot = FALSE,
                        type = "l",
                        lty = c(1,2),
                        col = c(2,4,8),
                        flab = NULL,
                        alab = "Cumulated amplitude",
                        flim = NULL, 
                        alim = NULL,
                        title = TRUE,
                        legend = TRUE,
                        ...
                        )

    {
        leg<-c(as.character(deparse(substitute(spec1))),as.character(deparse(substitute(spec2))))

        if(is.matrix(spec1) && ncol(spec1)==2) s1<-spec1[,2] else s1 <- spec1
        if(is.matrix(spec2) && ncol(spec2)==2) s2<-spec2[,2] else s2 <- spec2

        n1<-length(s1)
        n2<-length(s2)

        if(n1 != n2) stop("spec1 and spec2 must have the same length")
        
        if(any(is.na(s1)) | any(is.na(s2)))
            {
                D <- NA
                warning("The data set contains 'NA' values. The returned values have been set to NA.", call.=FALSE)
                return(D)
            } 
        else
            {
                if(any(s1 < 0) | any(s2 < 0)) stop("spectra (spec 1 and/or spec 2) do not have to be in dB")
                if(sum(s1)==0) {warning(paste("Caution!, spec1 is a null spectrum"), call.=FALSE)}
                if(sum(s2)==0) {warning(paste("Caution!, spec2 is a null spectrum"), call.=FALSE)}

		s1<-s1/sum(s1)
  		s2<-s2/sum(s2)
                cum.s1 <- cumsum(s1)
                cum.s2 <- cumsum(s2)     
                D <- sum(abs(cum.s1 - cum.s2))/n1

                if(plot)
                    {
                        if (!is.null(f) & mel) {f <- 2*mel(f/2)}
                        if(is.null(f))
                            {
                                if(is.vector(spec1) & is.vector(spec2)) stop("'f' is missing")  
                                else
                                    {
                                        if(is.matrix(spec1)) f<-spec1[nrow(spec1),1]*2000*nrow(spec1)/(nrow(spec1)-1)
                                        else if(is.matrix(spec2)) f<-spec2[nrow(spec2),1]*2000*nrow(spec2)/(nrow(spec2)-1)
                                    }
                            }
			x  <- seq(0, f/2*(n1-1)/n1, length.out=n1)/1000
                        st <- 0
                        en <- (f/2*(n1-1)/n1)/1000
      			if(is.null(alim)) alim<-c(0,1)
      			if(is.null(flim)) flim<-c(0,f/2000)
			if(is.null(flab)) {if (mel) flab <- "Frequency (kmel)" else flab <- "Frequency (kHz)"}
                        plot(x=x, y=cum.s1,
                             col=col[1], lty=lty[1], type = "n",
                             xlim = flim, xaxs = "i", xlab = flab, 
                             ylim = alim, yaxs = "i", ylab = alab, ...)
                        polygon(x=c(seq(st,en,length.out=n1),seq(en,st,length.out=n1)),
                                y=c(cum.s1,rev(cum.s2)),
                                col=col[3], border=NA)
                        lines(x=x, y=cum.s1, type=type, lty=lty[1], col=col[1])
                        lines(x=x, y=cum.s2, type=type, lty=lty[2], col=col[2])
                        if(legend) legend("topleft", col=c(col[1],col[2]) , lty=c(lty[1], lty[2]), legend=leg, bty="n")
                        if(title) title(main=paste("D = ", round(D, 3)))
                        invisible(D)
                    }
                else return(D)
            }
    }


################################################################################
##                                DIFFENV                                        
################################################################################

diffenv <- function(
                  wave1,
                  wave2,
                  f,
                  envt = "hil",
                  msmooth = NULL,
                  ksmooth = NULL,
                  plot = FALSE,
                  lty1 = 1,
                  lty2 = 2,
                  col1 = 2,
                  col2 = 4,
                  cold = 8,
                  xlab = "Time (s)",
                  ylab = "Amplitude",
                  ylim = NULL,
                  legend = TRUE,
                  ...
                  )

{
  leg<-c(as.character(deparse(substitute(wave1))),as.character(deparse(substitute(wave2))))

  input1<-inputw(wave=wave1,f=f) ; wave1<-input1$w ; f<-input1$f ; rm(input1)
  wave2<-inputw(wave=wave2,f=f)$w

  env1<-env(wave=wave1,f=f,envt=envt,msmooth=msmooth,ksmooth=ksmooth,plot=FALSE)
  env2<-env(wave=wave2,f=f,envt=envt,msmooth=msmooth,ksmooth=ksmooth,plot=FALSE)

  n1<-length(env1)
  n2<-length(env2)
  if(n1 != n2) stop("wave1 and wave2 should have the same length")

  if(!is.null(msmooth)|!is.null(ksmooth)) {f<-f*n1/nrow(wave1)}

  envPMF1<-env1/sum(env1)
  envPMF2<-env2/sum(env2)

  denv<-sum(abs(envPMF1-envPMF2))/2

  if(plot)
    {
      x<-seq(0,n1/f,length.out=n1)
      if(is.null(ylim)) ylim<-c(0,max(envPMF1,envPMF2))
      plot(x=x, y=envPMF1, type="n",
           xaxs="i", xlab = xlab, 
           ylim=ylim, yaxs="i", ylab=ylab, ...)
      polygon(x=c(x,rev(x)),
              y=c(envPMF1,rev(envPMF2)),
              col=cold,
              border=NA)
      lines(x, envPMF1, lty=lty1, col=col1)
      lines(x, envPMF2, lty=lty2, col=col2)
      box()
      if(legend) {legend("topleft", col=c(col1,col2),lty=c(lty1,lty2),legend=leg)}
      invisible(denv)
    }
  return(denv)
}


################################################################################
##                                DIFFSPEC                                         
################################################################################

diffspec <- function(spec1,
                     spec2,
                     f = NULL,
                     mel = FALSE,
                     plot = FALSE,
                     type = "l",
                     lty = c(1,2),
                     col = c(2,4,8),
                     flab = NULL,
                     alab = "Amplitude",
                     flim = NULL,
                     alim = NULL,
                     title = TRUE,
                     legend = TRUE,
                     ...
                     )

{
    leg<-c(as.character(deparse(substitute(spec1))),as.character(deparse(substitute(spec2))))

    if(is.matrix(spec1) && ncol(spec1)==2) s1<-spec1[,2] else s1 <- spec1
    if(is.matrix(spec2) && ncol(spec2)==2) s2<-spec2[,2] else s2 <- spec2

    n1<-length(s1)
    n2<-length(s2)

    if(n1 != n2) stop("spec1 and spec2 must have the same length")

    if(any(is.na(s1)) | any(is.na(s2)))
        {
            dspec <- NA
            warning("The data set contains 'NA' values. The returned values have been set to NA.", call.=FALSE)
            return(dspec)
        } 
    else
	{
            if(any(s1 < 0) | any(s2 < 0)) stop("spectra (spec 1 and/or spec 2) do not have to be in dB")
            if(sum(s1)==0) {warning(paste("Caution!, spec1 is a null spectrum"), call.=FALSE)}
            if(sum(s2)==0) {warning(paste("Caution!, spec2 is a null spectrum"), call.=FALSE)}

            s1<-s1/sum(s1)
            s2<-s2/sum(s2)

            dspec<-sum(abs(s1-s2))/2

            if(plot)
                {
                    if (!is.null(f) & mel) {f <- 2*mel(f/2)}
                    if(is.null(f))
                        {
                            if(is.vector(spec1) & is.vector(spec2)) stop("'f' is missing")  
                            else
                                {
                                    if(is.matrix(spec1)) f<-spec1[nrow(spec1),1]*2000*nrow(spec1)/(nrow(spec1)-1)
                                    else if(is.matrix(spec2)) f<-spec2[nrow(spec2),1]*2000*nrow(spec2)/(nrow(spec2)-1)
                                }
                        }
                    x<-seq(0, f/2*(n1-1)/n1, length.out=n1)/1000
                    st <- 0
                    en <- (f/2*(n1-1)/n1)/1000
                    if(is.null(alim)) alim<-c(0,max(s1,s2))
                    if(is.null(flim)) flim<-c(0,f/2000)
                    if(is.null(flab)) {if (mel) flab <- "Frequency (kmel)" else flab <- "Frequency (kHz)"}
                    plot(x=x, y=s1, type="n",
                         xlim=flim, xlab=flab, xaxs="i",
                         ylim=alim, ylab=alab, yaxs="i",...)
                    polygon(x=c(seq(st,en,length.out=n1),seq(en,st,length.out=n1)),
                            y=c(s1,rev(s2)),
                            col=col[3],
                            border=NA)
                    lines(x=x, y=s1, type=type, lty=lty[1], col=col[1])
                    lines(x=x, y=s2, type=type, lty=lty[2], col=col[2])
                    box()
                    if(legend) legend("topleft", col=c(col[1],col[2]),lty=c(lty[1],lty[2]),legend=leg, bty="n")
                    if(title) title(main=paste("dspec = ", round(dspec, 3)))
                    invisible(dspec)
                }
            else return(dspec)
        }
}



################################################################################
##                               DIFFWAVE                                        
################################################################################

diffwave<-function(
                   wave1,
                   wave2,
                   f,
                   wl = 512,
                   envt= "hil",
                   msmooth = NULL,
                   ksmooth = NULL
                   )

{
  input1<-inputw(wave=wave1,f=f) ; wave1<-input1$w ; f<-input1$f ; rm(input1)
  wave2<-inputw(wave=wave2,f=f)$w

                                        # spectral difference
  spec1<-meanspec(wave=wave1,f=f,wl=wl,PMF=TRUE,plot=FALSE)
  spec2<-meanspec(wave=wave2,f=f,wl=wl,PMF=TRUE,plot=FALSE)
  DF<-diffspec(spec1=spec1,spec2=spec2,f=f,plot=FALSE)

                                        # temporal difference 
  DE<-diffenv(wave1=wave1,wave2=wave2,f=f,msmooth=msmooth,ksmooth=ksmooth,plot=FALSE)

  z<-DF*DE
  return(z)
}


################################################################################
##                               DISCRETS                                        
################################################################################

discrets <- function (x,
                      symb = 5,
                      collapse = TRUE,
                      plateau = 1
                      ) 

{
    ## DATA INPUT
    if (symb != 3 && symb != 5) 
        stop("'symb' should be set to 3 or 5")
    x <- inputw(wave = x, f = NULL)$w
    n <- length(x)
    ## 5 SYMBOLS
    if(symb==5)
        {
            ## from the second point to the n-1 point
            s<-character(n-2)
            if(plateau == 1){     ## Cazelles et al definition, a plateau is encoded as PFP      
                for (i in 1:(n-2)){
                    if(x[i]<=x[i+1]   & x[i+1]<x[i+2])  s[i+1]<-"I"  # increase
                    if(x[i]<=x[i+2]   & x[i+2]<=x[i+1]) s[i+1]<-"P"  # peak
                    if(x[i+1]<x[i]    & x[i]<=x[i+2])   s[i+1]<-"T"  # trough
                    if(x[i+1]<x[i+2]  & x[i+2]<=x[i])   s[i+1]<-"T"  # trough
                    if(x[i+2]<x[i]    & x[i]<=x[i+1])   s[i+1]<-"P"  # peak
                    if(x[i+2]<=x[i+1] & x[i+1]<x[i])    s[i+1]<-"D"  # decrease
                    if(x[i]==x[i+1]   & x[i+1]==x[i+2]) s[i+1]<-"F"  # flat
                }
            }
            else if (plateau == 2) { ## Lellouch definition, a plateau is encoded as IFD
                for (i in 1:(n - 2)) {
                    if (x[i] < x[i + 1] & x[i + 1] > x[i + 2])   s[i + 1] <- "P" # peak
                    if (x[i] > x[i + 1] & x[i + 1] < x[i + 2])   s[i + 1] <- "T" # trough
                    if (x[i] <= x[i + 1] & x[i + 1] <= x[i + 2]) s[i + 1] <- "I" # increase
                    if (x[i] >= x[i + 1] & x[i + 1] >= x[i + 2]) s[i + 1] <- "D" # decrease
                    if (x[i] == x[i + 1] & x[i + 1] == x[i + 2]) s[i + 1] <- "F" # flat
                }
            }
        }
    else if (symb == 3) {
        s <- character(n - 1)
        for (i in 1:(n - 1)) {
            if (x[i] == x[i + 1]) 
                s[i + 1] <- "F"
            if (x[i] < x[i + 1]) 
                s[i + 1] <- "I"
            if (x[i] > x[i + 1]) 
                s[i + 1] <- "D"
        }
    }
    s <- s[-1]
    if (collapse) {s <- paste(s, collapse = "")}
    return(s)
}


################################################################################
##                                DRAWENV
################################################################################

drawenv<-function(
                  wave,
                  f,
                  n=20,
                  plot=FALSE,
                  listen = FALSE,                  
                  output = "matrix"
                  )

{
                                        # input
  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; bit <- input$bit ; rm(input)

  wave<-wave/max(abs(wave))
  wave<-rmoffset(wave, f=f)

                                        # interactive graph
  oscillo(wave=wave,f=f)
  cat("choose points on the positive amplitude side of the wave\nto change the amplitude profile (amplitude envelope)\n")
  if(.Platform$OS.type == "windows") flush.console()
  coord<-locator(n=n,type="p",col=2)

                                        # coordinates ; ordered following x if positions are not localised in order along the x-time axis
  X<-coord$x ; x<-round(X[order(X)]*f)
  Y<-coord$y ; y<-Y[order(X)]
  if(any(X<0)) stop("point localization cannot be on the negative part of the time axis")
  if(any(X>(nrow(wave)/f))) stop("point localization cannot be outside the positive part of the time axis")
  if(any(Y<0)) stop("point localization cannot be on the negative part of the amplitude axis")

                                        # profile generation
  profile<-numeric(nrow(wave))
  profile[1:x[1]]<-seq(0,y[1],length.out=x[1])
  for(i in 1:(length(x)-1))
    {
      profile[x[i]:x[i+1]]<-seq(y[i],y[i+1],length.out=x[i+1]-x[i]+1)
    }
  profile[x[length(x)]:length(profile)]<-seq(y[length(x)],0,length.out=length(profile)-x[length(x)]+1)

                                        # new wave generation
  wave2<-rmam(wave,f=f)
  wave2<-wave2/max(abs(wave2))
  wave3<-wave2[,1]*profile

  wave3 <- outputw(wave=wave3, f=f, format=output)

  if(plot)
    {
      dev.new()
      oscillo(wave3,f=f)
      if(listen) {listen(wave3,f=f)}
      invisible(wave3)
    }
  else
    {
      if(listen) {listen(wave3,f=f)}
      return(wave3)
    }
}


################################################################################
##                               DRAWFILTER
################################################################################

drawfilter <- function(f,   
                       n = 256,
                       continuous = TRUE,
                       discrete = TRUE
                       )
    
{
    if (!(continuous || discrete)) {stop("continuous or discrete must be TRUE")}
    f <- f/2000 
    plot.new()
    par(xaxs="i", yaxs="i", las=1)  
    plot(x=NA, xlim=c(0,f), ylim=c(0,1), type="n", xlab="toto", ylab="tata") 
    grid() 
    axis(1)
    axis(2)
    d <- seq(0,f,length.out=n+1)[1:n]
    g <- array(0, n)
    if (continuous)
        {cat("Please draw a continuous shape\n") 
         if (.Platform$OS.type == "windows") 
             flush.console()
         a <- locator(n=n, type = "l")
         if (length(a$x) < 2) {warning ("Please choose at least two points for the continuous filter\n")}
         else
             {if (any (diff(a$x) <= 0)) {stop("Values should increase along the x (frequency) axis")} 
              e <- array(0, n)
              for (i in 1:n)
                  {if (length(which(d[i] >= a$x)) == 0)
                       {e[i] <- 0}
                   else e[i] <- max(which(d[i] > a$x))}
              for (i in 1:n)
                  {if (e[i] == 0 || e[i] == length(a$x)) {g[i] <- 0}
                   else {h <- a$y[e[i]] + (a$y[e[i]+1] - a$y[e[i]])*(d[i] - a$x[e[i]])/(a$x[e[i]+1] - a$x[e[i]])
                         g[i] <- (h+abs(h))/2}
               }
          }
     }
    if (discrete)
        {cat("Please choose discrete frequencies\n")
         if (.Platform$OS.type == "windows") 
             flush.console()
         b <- locator(n=n, type = "p")
         if (length(b$x) < 1) {warning ("Locate at least one point for the discrete filter\n")}
         j <- round(b$x*n/f)
         j[j>=-5 & j < 0] <- 0
         j[j==n] <- n-1
         k <- b$y
         p <- j[j>=0 & j<=n-1 & k>=0 & k<=1.03]
         q <- k[j>=0 & j<=n-1 & k>=0 & k<=1.03]
         ext <- length(b$x) - length(p)
         red <- 0
         m <- array(0,n)
         if (length(p) >=1)
             {for(i in 1:length(p))
                  {if (m[p[i]+1] != 0) {red <- red + 1}
                   m[p[i]+1] <- m[p[i]+1] + q[i]
               }}
         if (ext > 0) {warning (cat(as.character(ext), "Discrete points out of the frame were not taken in account\n"))}
         if (red > 0) {warning (cat(as.character(red), "Duplications in discrete frequencies were found. For such frequencies, the sum of input amplitudes was used\n"))}
         g[m>0] <- m[m>0]}
    return(cbind(freq=d,amp=g))  ## JS changed the name of the columns to make it more clear
}


################################################################################
##                               DURATION
################################################################################

duration <- function(wave, f){
    input <- inputw(wave=wave,f=f)
    return(length(input$w) / input$f)
}

################################################################################
##                               DYNOSCILLO
################################################################################

dynoscillo <- function(
                  wave,
                  f,
                  wd = NULL,
                  wl = NULL,
                  wnb = NULL,
                  title = TRUE,
                  ...)                  

{
# STOP MESSAGES
  if(is.null(wl) && is.null(wnb) && is.null(wd)) stop("Either 'wd', 'wl' or 'wnb' has to be set.")
  if(!is.null(wl) && !is.null(wnb) && !is.null(wd)) stop("'wd', 'wl', 'wnb' cannot be used in the same time. Please choose one.")
  if(!is.null(wl) && !is.null(wnb)) stop("'wl' and 'wnb' cannot be used in the same time. Please choose one.")
  if(!is.null(wl) && !is.null(wd)) stop("'wl' and 'wd' cannot be used in the same time. Please choose one.")
  if(!is.null(wnb) && !is.null(wd)) stop("'wnb' and 'wd' cannot be used in the same time. Please choose one.")
  
# INPUT
  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; rm(input)

# SECTIONS
  n <- nrow(wave)
  if(!is.null(wd)) {wl <- wd*f}
  if(!is.null(wnb)) {wl <- round(n/wnb); wd <- wl/f}
  if(!is.null(wf)) {wd <- wl/f}
  step <- seq(from = 0, to = n-wl, by = wl)/f
  lstep <- length(step)
  pos <- 1:lstep
  
# PLOT  
      plot.osc<-function(panel)
        {
          with(panel,
                  {
                     max <- max(abs(wave))
                     soscillo(wave = wave, f = f, from = step[pos], to = step[pos]+wd, ylim=c(-max, max), tickup=max,...)
                     if(title) title(main=paste(pos,"/",lstep, "\n [wd = ", round(wd,3), " s, ", "wl = ", wl, "]", sep=""))
                  }
               )
          panel
        }

      osc.panel <- rpanel::rp.control("Position")
      rpanel::rp.slider(osc.panel, pos, from = 1, to = lstep, resolution = 1,
                title = "Position along the signal", action = plot.osc)
}


################################################################################
##                               DYNSPEC
################################################################################

dynspec<-function(
                  wave,
                  f,
                  wl = 512,
                  wn = "hanning",
                  zp = 0,
                  ovlp = 0,
                  fftw = FALSE,
                  norm = FALSE,
                  dB = NULL,
                  dBref = NULL,
                  plot = TRUE,
                  title = TRUE,
                  osc = FALSE,
                  flab = "Frequency (kHz)",
                  alab = "Amplitude",
                  alim = NULL,
                  flim = c(0,f/2000),
                  type ="l",
                  from = NULL,
                  to = NULL,
                  envt = NULL,
                  msmooth = NULL,
                  ksmooth = NULL,
                  colspec = "black",
                  coltitle = "black",
                  colbg = "white",
                  colline = "black",
                  colaxis = "black",
                  collab = "black",
                  cexlab = 1,
                  fontlab = 1,
                  colwave = "black",
                  coly0 = "lightgrey",
                  colcursor = "red",
                  bty = "l"
                  )

{
                                        # STOP MESSAGES
  if(is.logical(dB)) stop("'dB' is no more a logical. Please see the documentation: help(spec).")
  if(!is.null(dB) && all(dB!=c("max0","A","B","C","D")))
    stop("'dB' has to be one of the following character strings: 'max0', 'A', 'B', 'C' or 'D'")

                                        # INPUT
  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; rm(input)

  if(!is.null(from)|!is.null(to))
    {
      if(is.null(from) && !is.null(to)) {a<-1; b<-round(to*f)}
      if(!is.null(from) && is.null(to)) {a<-round(from*f); b<-length(wave)}
      if(!is.null(from) && !is.null(to))
        {
          if(from > to) stop("'from' cannot be superior to 'to'")
          if(from==0) {a<-1} else a<-round(from*f)
          b<-round(to*f)
        }
      wave<-as.matrix(wave[a:b,])
    }

  n <- nrow(wave)
  step <- seq(1,n-wl,wl-(ovlp*wl/100))
  lstep <- length(step)
                                        # STFT
  z <- stft(wave=wave, f=f, wl=wl, zp=0, step=step, wn=wn, norm=norm, fftw=fftw)
 
                                        # DB
  if(!is.null(dB))
    {
      if(is.null(dBref)) z<-20*log10(z) else z<-20*log10(z/dBref)
      if(dB == "max0") z <- z
      if(dB == "A") z <- dBweight(x*1000, dBref = z)$A 
      if(dB == "B") z <- dBweight(x*1000, dBref = z)$B 
      if(dB == "C") z <- dBweight(x*1000, dBref = z)$C 
      if(dB == "D") z <- dBweight(x*1000, dBref = z)$D
    }
  
                                        # FREQUENCY DATA 
  x<-seq(0, (f/2)-(f/wl),length.out=nrow(z))/1000

  if(plot)
    {
      if(is.null(alim))
        {
          alim<-c(0,max(z)+0.05)
          if(norm && is.null(dB)) alim<-c(0,1.1)
          if(!is.null(dB)) alim<-c(min(z),10)
        }     
      pos<-1:lstep
      poslabel<-numeric(lstep)
      for (i in 1:lstep){poslabel[i]<-round((step[i]+((step[2]-step[1])/2))/f,3)}

      plot.spec<-function(panel)
        {
          with(panel,
               {
                 par(bg=colbg, col=colline)
                 if(osc | !is.null(envt)) {layout(c(1,2),heights=c(2.5,1)); par(mar=c(4.5,4,3,2))}
                 plot(x=x,y=z[,pos],
                      xaxs = "i", xlab = flab, xlim = flim,
                      yaxs = "i", yaxt = "s", ylab = alab, ylim = alim,
                      col = colspec, col.axis = colaxis,
                      col.lab = collab, cex.lab = cexlab, font.lab = fontlab,
                      type = type, las = 1)

                 if(title)
                   {
                     nc <- ncol(z)
                     title(main=paste(pos, "/", nc, " - Position along the signal = ", poslabel[pos], "s", sep=""),
                           col.main=coltitle)
                   }
                 if(title==FALSE) title(main="")
                 if(is.character(title)) title(main = paste(title), col.main = coltitle)

                 if(osc)
                   {
                     par(mar=c(4.5,4,0.5,2))
                     max <- max(abs(wave))
                     soscillo(wave = wave, f = f,
                              colwave = colwave, collab = collab,
                              cexlab = cexlab, fontlab = fontlab, colline = colline,
                              colaxis = colaxis, coly0 = coly0, bty = bty,
                              tickup=max(abs(wave),na.rm=TRUE),
                              ylim=c(-max,max))
                     abline(v=poslabel[pos], col=colcursor)
                   }
                 else if(!is.null(envt))
                   {
                     par(mar=c(4.5,4,0.5,2))
                     env(wave = wave, f = f, k=1, j=1,
                         envt = envt, msmooth = msmooth, ksmooth = ksmooth,
                         colwave = colwave, collab = collab,
                         cexlab = cexlab, fontlab = fontlab, colline = colline,
                         colaxis = colaxis, coly0 = coly0, bty = bty)
                     abline(v=poslabel[pos], col=colcursor)
                   }
               }
               )
          panel
        }

      spec.panel <- rpanel::rp.control("Position")
      rpanel::rp.slider(spec.panel, pos, from=1, to=lstep, resolution=1,
                title = "Position along the signal", action=plot.spec)
    }
  else return(list(time=seq(0,n/f,length.out=length(step)), freq=x, amp=z))
}


################################################################################
##                               ECHO
################################################################################

echo <- function(
               wave,
               f,
               amp,
               delay,
               plot= FALSE,
               listen = FALSE,               
               output = "matrix",
               ...
               )

{
  cat("Please wait...\n")
  if(.Platform$OS.type == "windows") flush.console()

  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; bit <- input$bit ; rm(input)

  wave<-wave/max(abs(wave))
  n<-nrow(wave)
  namp<-length(amp)
  ndelay<-length(delay)
  delayp<-delay*f
  if(namp!=ndelay) stop("'namp' and 'ndelay' arguments should have the same length")
  if(any(namp)>1)  stop("'namp' argument cannot be > 1")
  if(any(namp)<0)  stop("'namp' argument cannot be negative")

  pulse<-c(1,numeric(delayp[namp]+n-1))
  for(i in 1:namp) pulse[delayp[i]]<-amp[i]
  pulse<-rev(pulse)
  wave2<-convolve(wave[,1],pulse,type="open")
                                        # delete the zero padded at the end of the wave by convolve()
  wave3<-wave2[1:(delayp[namp]+n)]

  wave3 <- outputw(wave=wave3, f=f, format=output)

  if(plot)
    {
      oscillo(wave=wave3,f=f,...)
      if(listen) {listen(wave3,f=f)}
      invisible(wave3)
    }
  else
    {
      if(listen) {listen(wave3,f=f)}
      return(wave3)
    }
}


################################################################################
##                               ENV                                        
################################################################################

env <- function(wave,
                 f,
                 envt = "hil",
                 msmooth = NULL,   
                 ksmooth = NULL,
                 ssmooth = NULL,
                 fftw = FALSE,
                 norm = FALSE,	
                 plot = TRUE,
                 k=1,
                 j=1,
                 ...)

{
    ## STOP MESSAGES
    if(!is.null(msmooth) & !is.null(ksmooth) & !is.null(ssmooth)) stop("Please use one of the smoothing arguments, not all of them together.") 
    if(!is.null(msmooth) & !is.null(ksmooth)) stop("Please use one of the smoothing arguments, not two of them.") 
    if(!is.null(msmooth) & !is.null(ssmooth)) stop("Please use one of the smoothing arguments, not two of them.") 
    if(!is.null(ksmooth) & !is.null(ssmooth)) stop("Please use one of the smoothing arguments, not two of them.") 

    ## INPUT
    input <- inputw(wave=wave,f=f) ; wave <- input$w ; f <- input$f ; rm(input)
    n <- nrow(wave)

    if(envt=="hil"){
        wave <- Mod(hilbert(wave,f=f,fftw=fftw))
    } 
    else if(envt=="abs"){
        if(fftw==TRUE && is.null(ksmooth)) warning("Setting 'fftw=TRUE' and 'ksmooth=NULL' has no meaning with 'envt=abs'.")
        wave <- abs(wave)
    }

    if(!is.null(msmooth))
        {
            if(msmooth[1] == 0)   stop("'smooth' window length cannot be equal to 0")
            if(msmooth[1] == 1)   stop("'smooth' window length cannot be equal to 1")
            if(msmooth[2] == 100) stop("'smooth' window overlap cannot be equal to 100")
            step<-seq(1,n-msmooth[1],msmooth[1]-(msmooth[2]*msmooth[1]/100))
            wave<-apply(as.matrix(step),1,function(x) mean(wave[x:(x + msmooth[1])]))
            wave<-as.matrix(wave)
            f<-f*(nrow(wave)/n)
        }

    if(!is.null(ksmooth))
        {
            if(fftw == TRUE){                                             
                kernapply.fftw <- function (x, k) {
                    if (!is.vector(x)) stop("'x' is not a vector") 
                    if (!is.tskernel(k)) stop("'k' is not a kernel")
                    m <- k$m
                    if (length(x) <= 2L * m)
                        stop("'x' is shorter than kernel 'k'")
                    if (m == 0L) return(x)
                    else {
                        n <- length(x)
                        w <- c(k[0L:m], rep_len(0, n - 2L * m - 1L), k[-m:-1L])
                        y <- fftw::IFFT(fftw::FFT(as.double(x)) * fftw::FFT(as.double(w)))/n
                        if (is.numeric(x)) y <- Re(y)
                        return(y[(1L + m):(n - m)])
                    }
                }
            wave <- kernapply.fftw(as.vector(wave), ksmooth)
            }
            else{wave <- kernapply(as.matrix(wave), ksmooth)}
            wave <- as.matrix(wave)
            f <- f*(nrow(wave)/n)
        }

    if(!is.null(ssmooth))
        {
            wave <- as.matrix(sumsmooth(wave, wl=ssmooth))
        } 

    if(plot)
        {
            oscillo(wave=wave,f=f,k=k,j=j,...)
            if(norm) wave<-wave/max(abs(wave))
            invisible(wave)
        }
    else
        {
            if(norm) wave<-wave/max(abs(wave))
            return(wave)
        }
}


################################################################################
##                               EXPORT
################################################################################

export<-function(
                 wave,
                 f = NULL,
                 filename = NULL,
                 header = TRUE, 
                 ...)

{
  if(is.null(filename)) {filename <- paste(as.character(deparse(substitute(wave))),".txt",sep="")}

  input<-inputw(wave=wave,f=f) ; wave<-input$w ;
  if(is.null(f)) {f<-input$f}
  rm(input)

  wave<-wave/(max(abs(wave))*1.5) # this avoids overclipping problems
  n<-nrow(wave)
  if(header==TRUE) {header<-paste("[ASCII ",f,"Hz, Channels: 1, Samples: ",n,", Flags: 0]", sep="")}
  else 
    {
      if(header==FALSE) header<-FALSE
      if(is.character(header)) header<-header
    }
  write.table(x=wave, file=filename, row.names=FALSE, col.names=header, quote=FALSE, ...)
}


################################################################################
##                                FADEW
################################################################################

fadew<-function(
                wave,
                f,
                din = 0,
                dout = 0,
                shape = "linear",
                plot = FALSE,
                listen = FALSE,
                output = "matrix",
                ...
                )

{
  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; bit <- input$bit ; rm(input)

  wave<-wave/max(abs(wave))
  n<-nrow(wave)
  ndin<-din*f
  ndout<-dout*f
  nfade<-ndin+ndout

  if(din==0 && dout==0) stop("Please specify at least a fade in or a fade out duration.")
  if(nfade>n) stop("The sum of fade in and fade out durations cannot be longer than wave length.")
  if(ndin>n) stop("Fade in duration cannot be longer than wave length.")
  if(ndout>n) stop("Fade in duration cannot be longer than wave length.")

  IN<-seq(0,1,length=ndin)
  OUT<-seq(0,1,length=ndout)

  if(shape=="exp")
    {
      IN<-exp(IN) ; IN<-IN-1 ; IN<-IN/max(IN)   # pb si din ou dout sont ==0
      OUT<-exp(OUT) ; OUT<-OUT-1 ; OUT<-OUT/max(OUT)
    }

  if(shape=="cos")
    {
      if(din == 0) IN<-integer(0)
      else {IN<-cos(rev(IN)) ; IN<-IN-min(IN) ; IN<-IN/max(IN)}   # pb si din ou dout sont ==0
      if(dout == 0) OUT<-integer(0)
      else {OUT<-cos(rev(OUT)) ; OUT<-OUT-min(OUT) ; OUT<-OUT/max(OUT)}
    }

  MID<-rep(1,nrow(wave)-(length(IN)+length(OUT)))

  FADE<-c(IN,MID,rev(OUT))

  wave2<-wave*FADE
  wave2<-wave2/max(abs(wave2))

  wave2 <- outputw(wave=wave2, f=f, format=output)

  if(plot)
    {
      oscillo(wave=wave2,f=f,...)
      if(listen) {listen(wave2,f=f)}
      invisible(wave2)
    }
  else
    {
      if(listen) {listen(wave2,f=f)}
      return(wave2)
    }
}


################################################################################
##                                FBANDS
################################################################################

fbands <- function(spec,
                   f = NULL,
                   bands = 10,
                   width = FALSE,
                   mel = FALSE,
                   plot=TRUE,
                   xlab= NULL,
                   ylab = "Relative amplitude",
                   ...
                   )
    
    {
        ## stop messages and data
        if(is.matrix(spec) && ncol(spec)!=2){stop("If 'spec' is a numeric matrix it should be a two-column matrix with the first colum describing the frequency x-axis and the second column describing the amplitude y-axis")}
        
        if(is.vector(spec))
            {
                ##                if (!is.null(f) & mel) {f <- 2*mel(f/2)}  ## supression JS
                if(is.null(f)){stop("If 'spec' is a numeric vector describing the amplitude only, the sampling frequency 'f' of the original signal should be provided (for instance (f = 44100)")}
                ##                if (mel) {f <- 2*mel(f/2)} ## supression JS
                N <- length(spec)
                spec <- cbind(seq(0, f*(N-1)/(2*N), length=N)/1000, spec)
            }

        if(is.null(bands)) stop("The argument 'bands' cannot be NULL.")
        if(any(bands<0)) stop("The argument 'bands' cannot include any negative value")
        
        n <- nrow(spec)
        
        if(length(bands)==1)     # a number of windows with a similar length
            {
                if(n/bands <= 2) {stop("Decrease the number of frequency bands")}
                bands <- seq(spec[1,1], spec[nrow(spec),1], length.out=bands+1)
            }
        
        k <- length(bands)
        res <- names <- freq <- wid <- numeric(k-1)
        for (i in 1:(k-1)){
            tmp <- cutspec(spec, flim=c(bands[i],bands[i+1]))
            ## remove the last line to avoid overlap between successive bands
            ## except for the last band :   [0,1[  [1,2[ ... [n, f/2]
            if(i < k-1) {tmp <- tmp[-nrow(tmp),]} 
            res[i] <- sum(tmp[,2])
            if(i==k-1) {names[i] <- paste("[", round(bands[i],1),"-", round(bands[i+1],1),"]",sep="")}
            else {names[i] <- paste("[", round(bands[i],1),"-", round(bands[i+1],1),"[",sep="")}
            freq[i] <- round(mean(c(bands[i],bands[i+1])),1)
            if(width) wid[i] <- bands[i+1]-bands[i] else wid <- 1
        }
        res <- cbind(freq, res/sum(res))
        colnames(res) <- c("freq","height")
        if(plot){if (is.null(xlab)) {if (mel) xlab <- "Frequency (kmel)" else xlab <- "Frequency (kHz)"}
                 barplot(height=res[,2], names.arg=names, width=wid , xlab=xlab, ylab=ylab,...)
                 
                 invisible(res)
             }
        else return(res)
    }



################################################################################
##                                FDOPPLER
################################################################################

fdoppler<-function(
                   f,
                   c = 340,
                   vs,
                   vo = 0,
                   movs = "toward",
                   movo = "toward"
                   )

{
  F<-f*((c-vo)/(c-vs))
  if(movs == "toward" && movo == "away") F<-f*((c+vo)/(c-vs))
  if(movs == "away" && movo == "toward") F<-f*((c-vo)/(c+vs))
  if(movs == "away" && movo == "away")   F<-f*((c+vo)/(c+vs))
  return(F)
}

################################################################################
##                                FFILTER
################################################################################

ffilter <- function(
                  wave,
                  f,
                  from = NULL,
                  to = NULL,
                  bandpass = TRUE,
                  custom = NULL,
                  wl = 1024,
                  ovlp = 75, 
                  wn="hanning",
                  fftw = FALSE,
                  output = "matrix"
                  )

{
  input <- inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; bit <- input$bit ; rm(input)

  n <- nrow(wave)
  step <- seq(1,n-wl,wl-(ovlp*wl/100))

  ## STFT
  z <- stft(wave=wave, f=f, wl=wl, zp=0, step=step, wn=wn, fftw=fftw, complex=TRUE, scale=FALSE)

  if(!is.null(custom))
    {
      if(is.matrix(custom)) custom <- custom[,2]
      if((length(custom)) != wl/2) stop("custom filter length has to be equal to 'wl'/2")
      z <-z*(custom/max(custom))
    }
  else
    {
      if(is.null(from) && is.null(to)) stop("At least one of the 'from' and 'to' arguments has to be set")
      if(is.null(from)) from<-0
      if(is.null(to)) to<-f/2
      if(from >= to) stop("'from' should be strictly inferior to 'to'")
      F<-ceiling(wl*(from/f))
      T<-floor(wl*(to/f))
      if(bandpass) {z[-c(F:T),]<-0} else {z[F:T,]<-0}
    }

  ## ISTFT
  res <- istft(z, wl=wl, ovlp = ovlp, wn=wn, output=output, f=f)
  
  return(res)
}


################################################################################
##                                FIELD
################################################################################

field<-function(f,d)
{
  c<-wasp(f=f)$c
  k<-f/c
  kd<-k*d
  if(length(d)==1)
    {
      if(kd<0.1) decision<-as.character("You are probably in the near-field, see documentation")
      if(kd>=0.1 & kd<1) decision<-as.character("You are probably at the limit between near-field and far-field, see documentation")
      if(kd>=1) decision<-as.character("You are probably in the far-field, see documentation")
      results<-list(kd=kd,d=decision)
    }
  else results<-list(kd=kd)
  return(results)
}


################################################################################
##                                FIR
################################################################################

fir<-function(
              wave,
              f,
              from = NULL,
              to = NULL,
              bandpass = TRUE,
              custom = NULL,
              wl = 512,
              wn = "hanning",
              listen = FALSE,
              output = "matrix"
              )

{
                                        # input
  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; bit <- input$bit ; rm(input)

                                        # frequency limits of the filter
  if(is.null(from)) from <- 0
  if(is.null(to)) to <- f/2
  from <- round((from * wl)/f)
  to <- round((to * wl)/f)
  n <- nrow(wave)

                                        # frequency response of the filter
  if(!is.null(custom))
    {
      if(is.matrix(custom)) custom<-custom[,2]
      if((length(custom)) != wl/2) stop("custom filter length has to be equal to 'wl'/2")
      if(bandpass)  {filtspec1 <- c(custom, rev(custom))}
      else                   {filtspec1 <- 1 - c(custom, rev(custom))}
      filtspec1 <- filtspec1/max(filtspec1)
    }
  else
    {
      filtspec1 <- rep(1, wl/2)
      if(bandpass)  {filtspec1[-(from:to)] <- 0}
      else                   {filtspec1[from:to] <- 0}
      filtspec1 <- c(filtspec1, rev(filtspec1))
    }

                                        # generation of filter pulse
  pulse1 <- Re(fft(filtspec1, inverse = TRUE)/length(filtspec1))
  pulse1 <- pulse1/max(pulse1)
  pulse2 <- c(pulse1[((wl/2)+1):wl], pulse1[-((wl/2+1):wl)])

                                        # window shape
  W <- ftwindow(wl = wl, wn = wn)

                                        # filter by convolution between the signal and pulse
  wave2<-convolve(wave[,1],pulse2*W,type="filter")

                                        # adds 0s before and after the signal to compensate for the reduction of wave length
  wave2<-c(rep(0,wl%/%2),wave2,rep(0,wl%/%2-1))

                                        # delete any potential offset
  wave2<-wave2-mean(wave2)

  wave2 <- outputw(wave=wave2, f=f, format=output)

  if(listen) {listen(wave2,f=f)}
  return(wave2)
}


################################################################################
##                                FMA
################################################################################

fma<-function(
              wave,
              f,
              threshold = NULL,
              plot = TRUE,
              ...)

{
  ifreq<-ifreq(wave, f=f, threshold=threshold, plot = FALSE)$f
  ifreq<-na.omit(ifreq)
  spec(ifreq[,2],f=f,plot=plot,...)
}


################################################################################
##                                FPEAKS
################################################################################

fpeaks <- function(spec,
                   f = NULL,
                   nmax = NULL,
                   amp = NULL,
                   freq = NULL,
                   threshold = NULL,
                   mel = FALSE,			
                   plot = TRUE,
                   title = TRUE,
                   xlab= NULL,
                   ylab = "Amplitude",
                   labels = TRUE,
                   legend = TRUE,
                   collab = "red",
                   ...
                   )
    
{
    ## STOP MESSAGES AND DATA
    if(is.matrix(spec))
        {
            if(ncol(spec)!=2) stop("If 'spec' is a numeric matrix it should be a two column matrix with the first colum describing the frequency x-axis and the second column describing the amplitude y-axis")
            N <- nrow(spec)
        }

    if(is.vector(spec))
        {
            N <- length(spec)
            if(is.null(f))
                {
                    stop("If 'spec' is a numeric vector describing the amplitude only, the sampling frequency 'f' of the original signal should be provided (for instance (f = 44100)")
                }
            if(!is.null(f) && !is.na(f))
                {
                    if (mel) {f <- 2*mel(f/2)}
                    spec <- cbind(seq(0, f/2-f/(N*2), length=N)/1000, spec)
                }
            if(!is.null(f) && is.na(f))
                {
                    if (mel) {f <- 2*mel(f/2)}
                    spec <- cbind(1:N,spec)
                    plot <- FALSE                               
                }
        }
    
    ## remove any flatness in the spec that would generate errors
    ## this is done by comparing successive points and if they have the same value a minor addition is operated
    flat <- round(N/20)     # the maximum number of sucessive points with the same amplitude (length of the flat part)
    spec.tmp <- c(spec[,2], rep(NA,flat))
    for(i in 1:(N-(flat+1)))
        {
            ref<-spec.tmp[i]
            for(j in 1:flat)
                {
                    if(spec.tmp[i+j]==ref)
                        {
                            spec.tmp[i+j]<-spec.tmp[i+j]+0.00001*spec.tmp[i+j]
                        }
                }
        }
    
    spec <- cbind(spec[,1],spec.tmp[1:N])

    ## DISCRETISATION
    sym <- discrets(spec[,2], symb=5, plateau=2, collapse=FALSE)
    if(sym[1]=="I") sym[1] <- "T"    # the sequence should start with a valley, otherwhise everything is shifted
    if(sym[1]=="P") sym[1] <-"D"     # the sequence should not start with a peak
    sym <- c(NA,sym,NA)  
    peaks <- which(sym=="P")
    n <- length(peaks)

    ## PEAKS SEARCH
    if(n==0) # no peaks
        {
            res <- NA
            plot <- FALSE
        }  
    else
        {
            if(!is.null(amp) | !is.null(nmax))    # left and right slope of the peak
                {
                    p <- numeric(n)
                    v <- numeric(n+1)
                    pic <- c(0, peaks, N+1)
                    for (i in 1:n)   {p[i] <- spec[peaks[i],2]}                      # peak i
                    for (i in 1:n+1) {v[i] <- min(spec[(pic[i]+1):(pic[i+1]-1), 2])} # valley i         
                    diffvp <- p - v[1:n]     # difference peaki - valleyi = left slope of the peak
                    diffpv <- p - v[2:(n+1)] # difference peaki - valleyi+1 = right slope of the peak
                }

            if(!is.null(nmax) && n!=0)
                {
                    if(!is.null(amp) | !is.null(freq) | !is.null(threshold)) {cat("Caution! The argument 'nmax' overrides the arguments 'amp', 'freq', and 'threshold'")}
                    if(n < nmax) {cat(paste("There are", n, "peaks only (< nmax ="), nmax,")")}
                    if(nmax==1)
                        {
                            tmp <- spec[peaks,, drop=FALSE]
                            res <- tmp[which.max(tmp[,2]), , drop=FALSE] 
                        }
                    else
                        {
                            alt <- cbind(peaks, diffvp, diffpv)
                            leftorder <- alt[order(-alt[,2]), , drop=FALSE]  # data ordered following the left peak slope
                            rightorder <- alt[order(-alt[,3]), , drop=FALSE] # data ordered following the right peak slope
                            left <- leftorder[,1]                            # left peak slopes ordered
                            right <- rightorder[,1]                          # right peak slopes ordered                                       
                            l <- 0             
                            i <- 1
                            while(l[i] < nmax)
                                {                                            # search matching between left and right peak slopes
                                    comp <- left[1:i] %in% right[1:i]        # returns the lowest number of left and right slopes for which slopes are max
                                    l <- c(l,length(comp[comp==TRUE]))
                                    i <- i+1
                                }
                            peaks0 <- left[1:(i-1)]
                            if(l[i] > nmax)
                                {
                                    error <- l[i] - nmax
                                    peaks0 <- peaks0[1:(length(peaks0)-error)]
                                }
                            peaks <- peaks0[comp]
                            res <- matrix(na.omit(spec[peaks,]), ncol = 2)   # remove NA, coerce into a two-column matrix
                            colnames(res) <- c("freq","amp")
                            res <-  res[order(res[,1]),]                     # reorder according to frequency values
                        }
                }
            ## AMPLITUDE SELECTION
            else
                {
                    if(!is.null(amp))  
                        {
                            if(length(amp)!=2) stop("The length of 'amp' should equal to 2.")
                            for(i in 1:n)
                                {
                                    if(!is.na(diffvp[i]) && !is.na(diffpv[i])
                                       && diffvp[i] > 0 && diffpv[i] > 0
                                       && diffvp[i] >= amp[1] && diffpv[i] >= amp[2])
                                        peaks[i] <- peaks[i]
                                    else
                                        peaks[i] <- NA
                                }
                        }
                    
                    ## FREQUENCY SELECTION
                    if (!is.null(freq)) {
                        freq <- freq/1000
                        diffpeak <- numeric(n - 1)
                        for (i in 1:(n - 1)) {
                            peak1 <- spec[peaks[i], 1]
                            peak2 <- spec[peaks[i + 1], 1]
                            diffpeak[i] <- peak2 - peak1
                            if (!is.na(diffpeak[i]) && diffpeak[i] <= freq) 
                                if (spec[peaks[i + 1], 2] > spec[peaks[i], 
                                                                 2]) 
                                    peaks[i] <- NA
                                else peaks[i + 1] <- NA
                        }
                    }
                    
                    if(!is.null(f) && is.na(f)) {res <-  peaks}
                    else
                        {
                            res <- matrix(na.omit(spec[peaks,]),ncol = 2) # remove NA, coerce into a 2 column matrix
                            colnames(res) <- c("freq","amp")
                        }
                    if(!is.null(threshold))
                        {
                            res <- res[res[,2]>threshold, ,drop=FALSE]      
                        }
                }
        }

    ## PLOT
    if(plot)
        {if (is.null(xlab)) {if (mel) xlab <- "Frequency (kmel)" else xlab <- "Frequency (kHz)"}


         plot(spec, type="l", xlab = xlab, ylab = ylab, xaxs="i", yaxt="n", ...)
         if(title) {
             if(nrow(res) == 1){text.title <- "peak detected"} else {text.title <- "peaks detected"}
             title(main=paste(nrow(res), text.title))
         }
         points(res, col = collab)
         if(labels & nrow(res)!=0) text(res, labels=round(res[,1],2), pos=3, col=collab)
         if(!is.null(threshold))
             {
                 abline(h=threshold, col=collab, lty=2)
                 mtext(paste(threshold), side=2, line=0.5, at=threshold, las=1, col=collab)
             }
         if(legend)
             {
                 if(!is.null(nmax)){text.legend <- paste("nmax=",nmax,sep="")}
                 else {
                     if(is.null(amp)){amp[1] <- amp[2] <- "-"} else amp <- round(amp,2)
                     if(is.null(freq)){freq <- "-"}
                     if(is.null(threshold)){threshold <- "-"}
                     text.legend <- c(paste("amp=", amp[1], "/", amp[2], sep=""),
                                      paste("freq=", freq, sep=""),
                                      paste("threshold=", threshold, sep="")
                                      )
                 }
                 legend("topright", pch=NA, legend = text.legend, bty="n", text.col="darkgrey")
             }
         invisible(res)
     }
    else return(res)
}


################################################################################
##                                FTWINDOW
################################################################################

ftwindow<-function(
                   wl,
                   wn = "hamming"
                   )

{
  if(wn=="bartlett")  w<-bartlett.w(wl)
  if(wn=="blackman")  w<-blackman.w(wl)
  if(wn=="flattop")   w<-flattop.w(wl)
  if(wn=="hamming")   w<-hamming.w(wl)
  if(wn=="hanning")   w<-hanning.w(wl)
  if(wn=="rectangle") w<-rectangle.w(wl)
  return(w)
}

################################################################################
##                                   FUND
################################################################################

fund <- function (
    wave,
    f,
    wl = 512,
    ovlp = 0,
    fmax,
    threshold = NULL, 
    at = NULL,
    from = NULL,
    to = NULL,
    plot = TRUE,
    xlab = "Time (s)",
    ylab = "Frequency (kHz)", 
    ylim = c(0, f/2000),
    pb = FALSE,
    ...)
    
{
    ## ERROR MESSAGES
    if(!is.null(at)){
        if(!is.null(threshold)) stop("The 'threshold' argument cannot be used with the argument 'at'.")
        if(ovlp!=0) stop("The 'overlap' argument should equal to 0 when using the argument 'at'.")
        if(!is.null(from) | !is.null(to)) stop("The 'from' and/or 'to' arguments cannot be used when using the argument 'at'.")
        if(pb) stop("No progress bar can be displayed when using the argument 'at'.")
        if(plot) {plot <- FALSE; warning("When the argument 'at' is used, the argument 'plot' is automatically turned to 'FALSE'.")}
    }
        
        ## INPUT
        input <- inputw(wave = wave, f = f)
        wave <- input$w
        f <- input$f
        rm(input)
        WL <- wl%/%2

        ## FROM-TO SELECTION
        if (!is.null(from) | !is.null(to)) {
            if (is.null(from) && !is.null(to)) {
                a <- 1
                b <- round(to * f)
            }
            if (!is.null(from) && is.null(to)) {
                a <- round(from * f)
                b <- length(wave)
            }
            if (!is.null(from) && !is.null(to)) {
                if (from > to) 
                    stop("'from' cannot be superior to 'to'")
                if (from == 0) {
                    a <- 1
                }
                else {
                    a <- round(from * f)
                }
                b <- round(to * f)
            }
            wave <- as.matrix(wave[a:b, ])
        }

        ## AT SELECTION
        if (!is.null(at)) {
            c <- round(at * f)
            wave <- as.matrix(wave[(c - WL):(c + WL), ]) 
        }
        
        ## THRESHOLD
        if (!is.null(threshold)) {
            wave <- afilter(wave = wave, f = f, threshold = threshold, 
                            plot = FALSE)
        }

        ## SLIDING WINDOW
        wave <- ifelse(wave == 0, yes = 1e-06, no = wave)
        n <- nrow(wave)
        step <- seq(1, n - wl, wl - (ovlp * wl/100))
        N <- length(step)
        z1 <- matrix(data = numeric(wl * N), wl, N)

        ## CEPSTRUM AND OPTIONAL PROGRESS BAR
        if (pb) {
            pbar <- txtProgressBar(min = 0, max = n, style = 3)
        }
        for (i in step) {
            z1[, which(step == i)] <- Re(fft(log(abs(fft(wave[i:(wl + 
                           i - 1), ]))), inverse = TRUE))
            if (pb) {
                setTxtProgressBar(pbar, i)
            }
        }

        ## FUNDAMENTAL FREQUENCY TRACKING
        z2 <- z1[1:WL, ]
        z <- ifelse(z2 == "NaN" | z2 == "-Inf", yes = 0, no = z2)
        z <- as.matrix(z)
        fmaxi <- f%/%fmax
        tfund <- numeric(N)
        for (k in 1:N) {
            tfund[k] <- which.max(z[-c(1:fmaxi), k])
        }
        tfund <- as.numeric(ifelse(tfund == 1, yes = NA, no = tfund))  ## "NA" passé en NA
        ffund <- f/(tfund + fmaxi - 1)
        if(!is.null(at)) {x <- at} else {x <- seq(0, n/f, length.out = N)}
        y <- ffund/1000
        res <- cbind(x, y)

        ## PLOT AND RETURN
        if (plot) {
            plot(x = x, y = y,
                 xaxs = "i", xlab = xlab,
                 yaxs = "i", ylab = ylab, ylim = ylim,
                 las = 1, ...)
            invisible(res)
        }       
        else {return(res)}
        if (pb) close(pbar)
    }


################################################################################
##                                   GGSPECTRO
################################################################################

ggspectro <- function(wave, f, tlab = "Time (s)", flab = "Frequency (kHz)", alab = "Amplitude\n(dB)\n", ...)
{
    ## data
    spectrogram <- spectro(wave, f=f, plot=FALSE, ...)
    frequency <- rep(spectrogram$freq, times=ncol(spectrogram$amp))
    time <- rep(spectrogram$time, each=nrow(spectrogram$amp))
    amplitude <- as.vector(spectrogram$amp)
    df <- data.frame(time, frequency, amplitude)
    ## ggplot object
    ## note: aes_string to pass R CMD check if not get a NOTE 'no visible global function definition for aes'
    v <- ggplot2::ggplot(df, ggplot2::aes_string(x="time", y="frequency", z = "amplitude")) + ggplot2::xlab(tlab) + ggplot2::ylab(flab) + ggplot2::scale_fill_continuous(alab)
    rm(spectrogram)
    return(v)
}


################################################################################
##                                   H
################################################################################

H<-function(
            wave,
            f,
            wl = 512,
            envt = "hil",
            msmooth = NULL,
            ksmooth = NULL
            )

{
  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; rm(input)
                                        # spectral entropy
  spec<-meanspec(wave=wave,f=f,wl=wl,plot=FALSE)
  SH<-sh(spec)

                                        # temporal entropy
  enve<-env(wave=wave,f=f,envt=envt,msmooth=msmooth,ksmooth=ksmooth,plot=FALSE)
  TH<-th(enve)

  z<-SH*TH
  return(z)
}



################################################################################
##                                HILBERT
################################################################################

hilbert <- function(wave,
                     f,
                     fftw = FALSE  ## ++
                     )    
{
    wave <- inputw(wave = wave,f = f)$w
    n <- nrow(wave)
    h <- rep(0,n)
    if(n > 0 & 2*floor(n/2) == n){
        h[c(1, n/2 + 1)] <- 1
        h[2:(n/2)] <- 2
    } else {
        h[1] <- 1;
        h[2:((n+1)/2)] <- 2
    }

    if(fftw == FALSE){ht <- fft(fft(wave)*h, inverse = TRUE)/n}            ## ++
    else {ht <- fftw::IFFT(fftw::FFT(as.double(wave))*h, scale = FALSE)/n} ## ++

    return(as.matrix(ht))
}


################################################################################
##                                IFREQ
################################################################################

ifreq<-function(
                wave,
                f,
                phase = FALSE,
                threshold = NULL,
                plot = TRUE,
                xlab = "Time (s)",
                ylab = NULL,
                ylim = NULL,
                type = "l",
                ...
                )

{
  if(class(wave)=="Wave") {f <- wave@samp.rate}
  wave<-hilbert(wave, f=f)
  n<-nrow(wave)
                                        # instantaneous phase
  phi<-Arg(wave)
                                        # instantaneous unwrapped phase
  phi2<-unwrap(phi[,1])
                                        # instantaneous frequency
  ifreq<-numeric(length(phi2)-1)
  for(i in 1:(length(phi2)-1)){ifreq[i]<-(f/1000)*(abs(phi2[i+1]-phi2[i]))/(2*pi)}
                                        # because the lenghth of ifreq is n-1 (the last point cannot be computed)
                                        # we build the n point as equals to the n-1 point 
  ifreq<-c(ifreq,ifreq[length(ifreq)-1])

  if(!is.null(threshold))
    {
      wavet<-afilter(wave=wave,f=f,threshold=threshold,plot=FALSE)
      phi[which(wavet[,1]==0)]<-"NA"
      ifreq[which(wavet[-n,1]==0)]<-"NA"
    }

  xf<-seq(0,n/f,length.out=n)
  xp<-seq(0,n/f,length.out=n)

  results <- list(f=cbind(xf,ifreq), p=cbind(xp,phi))
  
  if(plot)
    {
      if(phase == FALSE)
        {
          if(is.null(ylab)) {ylab<-"Frequency (kHz)"}
          if(is.null(ylim)) {ylim<-c(0,f/2000)}
          
          plot(x=xf, y=ifreq,
               xaxs="i", xlab=xlab,
               yaxs="i", ylab=ylab, ylim=ylim, 
               type=type, ...)
        }
      else
        {
          if(is.null(ylab)) {ylab<-"Phase (rad)"}
          if(is.null(ylim)) {ylim<-c(-pi,pi)}
          plot(x=xp, y=phi, 
               xaxs="i", xlab=xlab,
               yaxs="i", ylab=ylab, ylim=ylim, 
               type=type, ...)    
        }  
      invisible(results)
    }
  else return(results)
}


################################################################################
##                                ISTFT
################################################################################

istft <- function(
                  stft,
                  f,
                  wl,
                  ovlp = 75,
                  wn = "hanning",
                  output = "matrix"
                  )
    
    {
        ## stop messages
        if(!is.complex(stft)) stop("The object stft should be of complex mode (ie Re + Im.")

        h <- wl*(100-ovlp)/100  # get h (hop) parameter from ovlp (overlap)      
        coln <- ncol(stft)         
        xlen <- wl + (coln-1)*h   
        x <- numeric(xlen)          
        win <- ftwindow(wl=wl, wn=wn)

        ## perform IFFT and weighted-OLA
        for(b in seq(0, h*(coln-1), by=h)) 
            {
                ## Extract FFT points and build the mirror image of the FFT.
                ## The Nyquist point is built with the real part of the Nyquist-1 point and a O imaginary part
                X <- stft[, 1+b/h]   
                mirror <- rev(X[-1])
                mirror <- complex(real=Re(mirror), imaginary=-Im(mirror))
                X <- c(X, complex(real=Re(X[length(X)]), imaginary=0), mirror)
                ## IFFT
                xprim <- Re(fft(X, inverse = TRUE)/length(X))
                ## weighted-OLA  
                x[(b+1):(b+wl)] <-  x[(b+1):(b+wl)] + xprim*win
            }
        ## find W0
        W0 <- sum(win^2)
        ## scale the weighted-OLA
        x  <- x*h/W0                     
        ## return
        x <- outputw(wave = x, f = f, format = output)
    }


################################################################################
##                                ITAKURA.DIST
################################################################################

itakura.dist <- function(spec1,
                         spec2,
                         scale = FALSE
                          )

    {
        if(is.matrix(spec1) && ncol(spec1)==2) spec1<-spec1[,2]
        if(is.matrix(spec2) && ncol(spec2)==2) spec2<-spec2[,2]

        n1<-length(spec1)
        n2<-length(spec2)
        if(n1 != n2) stop("spec1 and spec2 must have the same length")

        if(any(is.na(spec1)) | any(is.na(spec2))) {
            D1 <- D2 <- D <- NA
            warning("The data set contains some 'NA' values. The returned values have been set to 'NA'.", call.=FALSE)
        } 
        else
            {
                if(any(spec1<0, na.rm=TRUE) | any(spec2<0, na.rm=TRUE)) stop("Data do not have to be in dB")
                if(sum(spec1)==0) {warning("Caution!, spec1 is a null spectrum", call.=FALSE)}
                if(sum(spec2)==0) {warning("Caution!, spec2 is a null spectrum", call.=FALSE)}
                spec1[spec1==0]<-1e-7
                spec2[spec2==0]<-1e-7
                spec1 <- spec1/sum(spec1)  # PMF
                spec2 <- spec2/sum(spec2)  # PMF
                r <- spec1/spec2
                D1 <- sum(r - log(r) - 1)/n1
                D2 <- sum(1/r - log(1/r) - 1)/n2
                if (scale == TRUE) {
                    D1 <- D1/n1
                    D2 <- D2/n2
                }
                D <- (D1+D2)/2
            }
        return(list(D1=D1, D2=D2, D=D))
    }



################################################################################
##                                Q
################################################################################

Q <- function (spec, f = NULL, level = -3, plot = TRUE, colval = "red",
    cexval = 1, fontval = 1, flab = "Frequency (kHz)", alab = "Relative amplitude (dB)",
    type = "l", ...)
{
    if (is.null(f)) {
        if (is.vector(spec))
            stop("'f' is missing")
        else if (is.matrix(spec))
            f <- spec[nrow(spec), 1] * 2000
    }
    if (is.matrix(spec))
        spec <- spec[, 2]
    range <- c(f/2000/length(spec), f/2000)
    if (max(spec) == 1)
        stop("data must be in dB")
    if (which.max(spec) == 1)
        stop("maximal peak cannot be the first value of the spectrum")
    n0 <- length(spec)
    spec1 <- approx(spec, n = 102400)
    spec1 <- as.matrix(spec1$y)
    n1 <- nrow(spec1)
    level2 <- round(max(spec1[, 1]), 1) + level
    f0 <- which.max(spec1[, 1])
    specA <- as.matrix(spec1[1:f0, 1])
    nA <- nrow(specA)
    specB <- as.matrix(spec1[f0:nrow(spec1), 1])
    f1 <- which(round(specA, 1) == level2)
    f1khz <- ((f1[length(f1)]/n1) * (range[2] - range[1])) +
        range[1]
    f2 <- which(round(specB, 1) == level2) + (nA - 1)
    f2khz <- ((f2[1]/n1) * (range[2] - range[1])) + range[1]
    Q <- f0/(f2[1] - f1[length(f1)])
    results <- list(Q=Q, fmin = f1khz, fmax = f2khz, bdw = f2khz-f1khz)
    if (plot) {
        x <- seq(range[1], range[2], length.out = n0)
        plot(x = x, y = spec, xlab = flab, ylab = alab, type = type,
            ...)
        arrows(f1khz, level2, f2khz, level2, length = 0.1, col = colval,
            code = 3, angle = 15)
        text(paste("Q =", as.character(round(Q, 2))), x = f2khz,
            y = level2, pos = 4, col = colval, cex = cexval,
            font = fontval)
        invisible(results)
    }
    return(results)
}



################################################################################
##                                KL.DIST
################################################################################

kl.dist <- function(
                    spec1,
                    spec2,
                    base = 2
                    )

  {
    if(is.matrix(spec1) && ncol(spec1)==2) spec1<-spec1[,2]
    if(is.matrix(spec2) && ncol(spec2)==2) spec2<-spec2[,2]

    n1<-length(spec1)
    n2<-length(spec2)
    if(n1 != n2) stop("spec1 and spec2 must have the same length")

    if(any(is.na(spec1)) | any(is.na(spec2))) {
      D1 <- D2 <- D <- NA
      warning("The data set contains some 'NA' values. The returned values have been set to 'NA'.", call.=FALSE)
    } 
    else
      {
        if(any(spec1<0, na.rm=TRUE) | any(spec2<0, na.rm=TRUE)) stop("Data do not have to be in dB")
        if(sum(spec1)==0) {warning("Caution!, spec1 is a null spectrum", call.=FALSE)}
        if(sum(spec2)==0) {warning("Caution!, spec2 is a null spectrum", call.=FALSE)}
        spec1[spec1==0]<-1e-7
        spec2[spec2==0]<-1e-7
        spec1 <- spec1/sum(spec1)  # PMF
        spec2 <- spec2/sum(spec2)  # PMF
        D1 <- sum(spec1 * log(spec1/spec2, base=base))
        D2 <- sum(spec2 * log(spec2/spec1, base=base))             
      }
    D  <- (D1+D2)/2    
    return(list(D1=D1, D2=D2, D=D))
  }


################################################################################
##                                KS.DIST
################################################################################

ks.dist <- function(spec1,
                    spec2,
                    f = NULL,
		    mel = FALSE,
                    plot = FALSE,
                    type = "l",
                    lty = c(1,2),
		    col = c(2,4),
                    flab = NULL,
                    alab = "Cumulated amplitude",
                    flim = NULL,
                    alim = NULL,
                    title = TRUE,
                    legend = TRUE,
                    ...
                    )
    
    {
        leg<-c(as.character(deparse(substitute(spec1))),as.character(deparse(substitute(spec2))))

        if(is.matrix(spec1) && ncol(spec1)==2) s1<-spec1[,2] else s1 <- spec1
        if(is.matrix(spec2) && ncol(spec2)==2) s2<-spec2[,2] else s2 <- spec2

        n1<-length(s1)
        n2<-length(s2)

        if(n1 != n2) stop("spec1 and spec2 must have the same length")

        if(any(is.na(s1)) | any(is.na(s2)))
            {
                D <- F <- NA
                warning("The data set contains 'NA' values. The returned values have been set to NA.", call.=FALSE)
                res <- list(D=D, F=F)
                return(res)
            } 
        else
            {
                if(any(s1 < 0) | any(s2 < 0)) stop("spectra (spec 1 and/or spec 2) do not have to be in dB")
                if(sum(s1)==0) {warning(paste("Caution!, spec1 is a null spectrum"), call.=FALSE)}
                if(sum(s2)==0) {warning(paste("Caution!, spec2 is a null spectrum"), call.=FALSE)}

                s1<-s1/sum(s1)
                s2<-s2/sum(s2)
                cum.s1 <- cumsum(s1)
                cum.s2 <- cumsum(s2)
                
                diff <- abs(cum.s1 - cum.s2)
                D <- max(diff)
                if (!is.null(f) & mel) {f <- 2*mel(f/2)}
                if(is.null(f))
                    {
                        if(is.vector(spec1) & is.vector(spec2))
                            {
                                if(plot) warning("'f' is missing, F cannot be computed - No plot produced", call.=FALSE) else warning("'f' is missing, F cannot be computed", call.=FALSE)
                                res <- list(D=D, F=NA)
                                return(res)
                            }
                        else
                            {
                                if(is.matrix(spec1)) f<-spec1[nrow(spec1),1]*2000*nrow(spec1)/(nrow(spec1)-1)
                                else if(is.matrix(spec2)) f<-spec2[nrow(spec2),1]*2000*nrow(spec2)/(nrow(spec2)-1)
                            }
                    }

                x <- seq(0,f/2*(n1-1)/n1, length.out=n1)/1000
                pos <- which.max(diff)
                F <- x[pos]
                res <- list(D=D, F=F)

                if(plot)
                    {if (mel) {w <- " kmel"}
                     else {w <- " kHz"}
                     if(is.null(alim)) alim<-c(0,1)
                     if(is.null(flim)) flim<-c(0,f/2000)
                     if(is.null(flab)) {if (mel) flab <- "Frequency (kmel)" else flab <- "Frequency (kHz)"}
                     x <- seq(0,(f/2)*(n1-1)/n1, length.out=n1)/1000
                     plot(x=x, y=cum.s1,
                          col=col[1], lty=lty[1], type = type,
                          xaxs="i", xlab=flab, xlim=flim,
                          yaxs="i", ylim=alim, ylab=alab,...)
                     lines(x=x, y=cum.s2, col=col[2], lty=lty[2], type = type)
                     segments(x0=F, y0= cum.s1[pos], x1=F, y1=cum.s2[pos], lwd=2)
                     if(title) title(main=paste("D = ", round(D, 3), "\n F = ", round(F, 3), w))
                     if(legend) legend("topleft", col=col, lty=lty, legend=leg, bty="n")
                     invisible(res)
                 }
                else return(res)
            }
    }


################################################################################
##                                LISTEN
################################################################################

listen <- function(
                 wave,
                 f = NULL,
                 from = NULL,
                 to = NULL,
                 choose = FALSE
                 )

{
  input <- inputw(wave=wave,f=f) ; wave <- input$w ;
  if(is.null(f)) {f <- input$f}
  rm(input)

  if(choose)
    { 
      cat("choose start and end positions on the wave\n")
      if(.Platform$OS.type == "windows") flush.console()
      oscillo(wave,f=f)
      coord<-locator(n=2)
      from<-coord$x[1]; a<-round(from*f) ; abline(v=from,col=2,lty=2)
      to<-coord$x[2]; b<-round(to*f); abline(v=to,col=2,lty=2)
      wave<-wave[a:b]
    }
  else if(!is.null(from)|!is.null(to))
    {
      if(is.null(from) && !is.null(to)) {a<-1; b<-round(to*f)}
      if(!is.null(from) && is.null(to)) {a<-round(from*f); b<-length(wave)}
      if(!is.null(from) && !is.null(to))
        {
          if(from>to) stop("'from' cannot be superior to 'to'")
          if(from==0) {a<-1} else a<-round(from*f)
          b<-round(to*f)
        }
      wave<-as.matrix(wave[a:b])
    }

  wave <- Wave(left=wave, samp.rate=f, bit=16)
  wave <- normalize(wave, unit="16")
  tuneR::play(wave)
}


################################################################################
##                                LFS
################################################################################

lfs <- function(
    wave,
    f,
    shift,
    wl = 1024,
    ovlp = 75,
    wn="hanning",
    fftw = FALSE,
    output = "matrix"
    )

{
    ## STOP MESSAGES
    input <- inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; bit <- input$bit ; rm(input)
    n <- nrow(wave)

    if(shift == 0) stop("'shift' value cannot be equal to 0")
    if(shift>f/2) stop("Positive 'shift' value cannot exceed half of the sampling frequency")
    if(shift<(-f/2)) stop("Negative 'shift' value cannot be less than half of the sampling frequency")
    if(abs(shift)<f/wl) stop("'shift' value cannot be less than the frequency resolution (f/wl)")
    if(wl>n*2) stop("'wl' value is too high, wl should less than length(wave)*2")

    step <- seq(1,n-wl,wl-(ovlp*wl/100))
    FSH <- abs(shift)

    ## STFT
    z1 <- stft(wave=wave, f=f, wl=wl, zp=0, step=step, wn=wn, fftw=fftw, complex=TRUE, scale=FALSE)

    S <- round(wl*(FSH/f))

    ## generate a 0 matrix corresponding to the frequency shift to apply
    z2a <- matrix(data=0, nrow=S, ncol=ncol(z1))

    ## first case: the frequency shift is positive
    if(shift > 0)
        {
            z2b <- z1[c(1:(wl/2-S)),]
            z2c <- rbind(z2a,z2b)
        }

    ## second case: the frequency shift is negative
    if(shift < 0)
        {
            z2b <- z1[-c(1:S),]
            z2c <- rbind(z2b,z2a)
        }

    ## ISTFT
    res <- istft(z2c, wl=wl, ovlp = ovlp, wn=wn, output=output, f=f)
    return(res)
}



################################################################################
##                                LOCALPEAKS
################################################################################

localpeaks <- function(spec,
                       f = NULL,
                       bands = 10,
		       mel = FALSE,
                       plot = TRUE,
                       xlab= NULL,
                       ylab = "Amplitude",
                       labels = TRUE, ...)
    
    {
        ## stop messages
        if(is.matrix(spec) && ncol(spec)!=2){stop("If 'spec' is a numeric matrix it should be a two-column matrix with the first colum describing the frequency x-axis and the second column describing the amplitude y-axis")}
        
        if(is.vector(spec))
            {
                if(is.null(f)){stop("If 'spec' is a numeric vector describing the amplitude only, the sampling frequency 'f' of the original signal should be provided (for instance (f = 44100)")}
                N <- length(spec)
                spec <- cbind(seq(0, f/2-f/(N*2), length=N)/1000, spec)
            }

        if(is.null(bands)) stop("The argument 'bands' cannot be NULL.")
        if(any(bands<0)) stop("The argument 'bands' cannot include any negative value")
        
        n <- nrow(spec)

        if(length(bands)==1)     # a number of windows with a similar length
            {
                if(n/bands <= 2) {stop("Decrease the number of frequency bands")}
                bands <- seq(spec[1,1], spec[nrow(spec),1], length.out = bands+1)
            }
        
        k <- length(bands)
        res <- matrix(numeric((k-1)*2), nrow = k-1, ncol = 2)
        for(i in 1:(k-1))
            {
                tmp <- cutspec(spec, flim=c(bands[i],bands[i+1]))
                ## remove the last line to avoid overlap between successive bands
                ## except for the last band :   [0,1[  [1,2[ ... [n, f/2]
                if(i < k-1) {tmp <-tmp[-nrow(tmp),]} 
                res[i,] <- fpeaks(tmp, nmax=1, plot=FALSE)
            }
        
        colnames(res) <- c("freq","amp")
        
        if(plot)
            {
                if (is.null(xlab)) {
                    if (mel) xlab <- "Frequency (kmel)"
                    else xlab <- "Frequency (kHz)"
                }
                plot(spec, type="l", xlab = xlab, ylab = ylab, xaxs="i", yaxt="n")
                points(res, col = "red")
                abline(v=bands, col="grey")
                if(labels) text(res, labels=round(res[,1],2), pos=3, col="red")
                box()
                invisible(res)
            }
        else {return(res)}
    }


################################################################################
##                                LOGSPEC.DIST
################################################################################

logspec.dist <- function(spec1,
                         spec2,
                         scale = FALSE
                         )

    {
        
        if(is.matrix(spec1) && ncol(spec1)==2) spec1<-spec1[,2]
        if(is.matrix(spec2) && ncol(spec2)==2) spec2<-spec2[,2]

        n1<-length(spec1)
        n2<-length(spec2)
        if(n1 != n2) stop("spec1 and spec2 must have the same length")

        if(any(is.na(spec1)) | any(is.na(spec2))) {
            D <- NA
            warning("The data set contains some 'NA' values. The returned values have been set to 'NA'.", call.=FALSE)
        } 
        else
            {
                if(any(spec1<0, na.rm=TRUE) | any(spec2<0, na.rm=TRUE)) stop("Data do not have to be in dB")
                if(sum(spec1)==0) {warning("Caution!, spec1 is a null spectrum", call.=FALSE)}
                if(sum(spec2)==0) {warning("Caution!, spec2 is a null spectrum", call.=FALSE)}
                spec1[spec1==0]<-1e-7
                spec2[spec2==0]<-1e-7
                spec1 <- spec1/sum(spec1)  # PMF
                spec2 <- spec2/sum(spec2)  # PMF
                D <- sqrt(sum((10*log10(spec1/spec2))^2))
                if(scale==TRUE) {D <- D/sqrt(n1)}
            }
        return(D)
    }


################################################################################
##                                M
################################################################################

M <- function(wave,
              f,
              envt = "hil",
              plot = FALSE,
              ...)
    
{
    ## INPUT
    input <- inputw(wave=wave,f=f) ; wave <- input$w ; f <- input$f ; bit <- input$bit ; rm(input)
    ## ENVELOPE
    env <- env(wave, f=f, envt=envt, plot=plot, ...)
    ## INDEX
    M <- median(env)*(2^(1-bit))    
    return(M)
}


################################################################################
##                                MEANDB
################################################################################

meandB<-function(x, level="IL")
{
  if(level=="IL") {a <- 10} else {a <- 20}
  return(a*log10(mean(10^(x/a))))
}


################################################################################
##                               MEANSPEC
################################################################################

meanspec <-function(
                    wave,
                    f,
                    wl = 512,
                    wn = "hanning",
                    ovlp = 0,
                    norm = TRUE, 
                    fftw = FALSE,
                    PSD =  FALSE,
                    PMF = FALSE,
                    dB = NULL,
                    dBref = NULL,
                    from = NULL,
                    to = NULL,
                    identify = FALSE,
                    col = "black",
                    cex = 1,
                    plot = 1,
                    flab = "Frequency (kHz)",
                    alab = "Amplitude",
                    flim = NULL,
                    alim = NULL,
                    type ="l",
                    ...)

{
                                        # STOP MESSAGES  
  if(!isTRUE(norm) & PMF) stop("'PMF' can be computed only if 'norm' is TRUE")
  if(!isTRUE(norm) & !is.null(dB)) stop("dB are computed on normalised spectra only, 'norm' should be turned to TRUE")
  if(!is.null(dB) & PMF) stop("'PMF' cannot be in 'dB'")
  if(!is.null(dB) & PSD) stop("PSD cannot be in dB")
  if(is.null(dB) & !is.null(dBref)) stop("'dB' cannot be NULL  when 'dBref' is not NULL")
  if(is.logical(dB)) stop("'dB' is no more a logical. Please see the documentation: help(spec).")
  if(!is.null(dB) && all(dB!=c("max0","A","B","C","D")))
    stop("'dB' has to be one of the following character strings: 'max0', 'A', 'B', 'C' or 'D'")
  if(!is.null(wl) & wl%%2 == 1) stop("'wl' has to be an even number.")
  
                                        # INPUT
  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; rm(input)

                                        # FROM-TO SELECTION
  if(!is.null(from)|!is.null(to))
    {
      if(is.null(from) && !is.null(to)) {a<-1; b<-round(to*f)}
      if(!is.null(from) && is.null(to)) {a<-round(from*f); b<-length(wave)}
      if(!is.null(from) && !is.null(to))
        {
          if(from>to) stop("'from' cannot be superior to 'to'")
          if(from==0) {a<-1} else {a<-round(from*f)}
          b<-round(to*f)
        }
      wave<-as.matrix(wave[a:b,])
    }

                                        # STFT
  n <- nrow(wave)
  step <- seq(1, n-wl, wl-(ovlp*wl/100))
  y <- stft(wave=wave, f=f, wl=wl, zp=0, step=step, wn=wn, fftw=fftw, scale=norm)
  y <- rowMeans(y) # mean computation by rows faster than  y<-apply(y,MARGIN=1,mean)     
  if(norm) {y <- y/max(y)}                               
  y <- ifelse(y==0,yes=1e-6,no=y) # replaces 0 values in spectra that cannot be processed by log10())

                                        # FREQUENCY DATA 
  x<-seq(0, (f/2)-(f/wl), length.out= wl%/%2) / 1000

                                        # PSD and PMF OPTIONS
  if(PSD) y<-y^2
  if(PMF) y<-y/sum(y)

                                        # DB
    if(!is.null(dB))
    {
      if(is.null(dBref)) {y <- 20*log10(y)} else {y <- 20*log10(y/dBref)}
      if(dB!="max0")
          {
            if(dB == "A") y <- dBweight(x*1000, dBref = y)$A 
            if(dB == "B") y <- dBweight(x*1000, dBref = y)$B 
            if(dB == "C") y <- dBweight(x*1000, dBref = y)$C 
            if(dB == "D") y <- dBweight(x*1000, dBref = y)$D
          }
    }
  
                                        # LIMITS OF THE AMPLITUDE AXIS
  if(is.null(alim))
    {
      if(is.null(dB)) {alim<-c(0,1.1)} else {alim <- c(min(y, na.rm=TRUE), max(y, na.rm=TRUE) + 20)}
      if(PMF | !isTRUE(norm)) alim<-c(0,max(y, na.rm=TRUE))
    }


                                          # HORIZONTAL PLOT
  if(plot == 1)
    {
      if(!is.null(dB))
        {
          plot(x=x, y=y,
               xaxs = "i", xlab = flab, xlim = flim,
               yaxs = "i", yaxt = "s", ylab = alab, ylim = alim,
               col = col, cex = cex,
               type = type, las = 1,
               ...)
        }
      else
        {         
          if(isTRUE(norm))
            {
              yaxt<-"n"
              ylab<-alab
              if(isTRUE(PMF)) {yaxt="s"}
            }
          else
            {
              yaxt<-"s"
              ylab<-" "
            }

          plot(x=x,y=y,
               xaxs="i", xlab=flab, xlim = flim,
               yaxs="i", yaxt=yaxt, ylab = ylab, ylim=alim,
               col = col, cex = cex,
               type = type, las = 1,
               ...)
        }

      if(identify)
        {
          cat("Choose points on the spectrum\n")
          if(.Platform$OS.type == "windows") flush.console()
          id<-identify(x=x,y=y,labels=round(x,2),tolerance=0.15,col="red")
          id.freq<-x[id]
          id.amp<-y[id]
          coord<-list(freq = id.freq ,amp = id.amp)
          return(coord)
        }
    }

                                        # VERTICAL PLOT
  if(plot == 2)
    {
      if(!is.null(dB))
        {
          plot(x=y,y=x,
               xaxs = "i", xlab = alab, xlim = alim,
               yaxs = "i", yaxt = "s", ylab = flab, ylim = flim,
               col = col, cex = cex,
               type = type, las = 1,
               ...)
        }
      else
        {
          if(isTRUE(norm))
            {
              xaxt <- "n"
              xlab <- alab
              if(isTRUE(PMF)) {xaxt = "s"}
            }
          else
            {
              xaxt<-"s"
              xlab<-" "
            }
          
          plot(x=y,y=x,
               xaxs = "i", xaxt = xaxt, xlab = xlab, xlim = alim,
               yaxs = "i", ylab = flab, ylim = flim,
               col = col, cex = cex,
               type = type, las = 1,
               ...)
        }

      if(identify)
        {
          cat("choose points on the spectrum\n")
          if(.Platform$OS.type == "windows") flush.console()
          id<-identify(x=y,y=x,labels=round(x,2),tolerance=0.15,col="red")
          id.freq<-x[id]
          id.amp<-y[id]
          coord<-list(freq = id.freq ,amp = id.amp)
          return(coord)
        }
    }


                                        # INVISIBLE RETURN DATA  
  if(plot == 1 | plot == 2)
    {
      spec<-cbind(x,y)	
      invisible(spec)
    }
                                        # DATA RETURN WHEN NO PLOT  
  else if(plot == FALSE) 
    {
      spec<-cbind(x,y)
      return(spec)
    }
}


################################################################################
##                                MEL
################################################################################

mel<-function(
              x,
              inverse = FALSE
              )

{
  y<-1127.01048*log(1+(x/700))
  if(inverse) y<-700*(exp(x/1127.01048)-1)
  return(y)
}


################################################################################
##                                MICSENS
################################################################################

micsens<-function(x,sref=1,inverse=FALSE){
  if(inverse==FALSE)
    {
      s<-x/1000
      S<-20*log10(s/sref)
    }
  else {S<-1000*sref*10^(x/20)}
  return(S)
}


################################################################################
##                                MOREDB
################################################################################

moredB<-function(x, level="IL")
{
  if(level=="IL") {a <- 10} else {a <- 20}
  return(a*log10(sum(10^(x/a))))
}


################################################################################
##                                MUTEW
################################################################################

mutew<-function(
                wave,
                f,
                from = NULL,
                to = NULL,
                choose = FALSE,
                plot = TRUE,
                output = "matrix",
                ...
                )

{
  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; bit <- input$bit ; rm(input)

  n<-nrow(wave)
  
  if(choose)
    { 
      cat("choose start and end positions on the wave\n")
      if(.Platform$OS.type == "windows") flush.console()
      oscillo(wave,f=f)
      coord<-locator(n=2)
      from<-coord$x[1]; a<-round(from*f) ; abline(v=from,col=2,lty=2)
      to<-coord$x[2]; b<-round(to*f); abline(v=to,col=2,lty=2)
      wave.muted<-as.matrix(c(wave[1:(a-1),],rep(0,length(a:b)),wave[(b+1):n,]))
    }

  else if(!is.null(from)|!is.null(to))
    {
      if(is.null(from) && !is.null(to))
        {
          b<-round(to*f)
          wave.muted<-as.matrix(c(rep(0,b),wave[(b+1):n,]))
        }

      if(!is.null(from) && is.null(to)) 
        {
          a<-round(from*f)
          wave.muted<-as.matrix(c(wave[1:(a-1),],rep(0,length(a:n))))
        }

      if(!is.null(from) && !is.null(to))
        {
          if(from > to) stop("'from' cannot be superior to 'to'")
          if(from == 0) {a<-1; b<-round(to*f)}
          else {
            a<-round(from*f)
            b<-round(to*f)}
          wave.muted<-as.matrix(c(wave[1:(a-1),],rep(0,length(a:b)),wave[(b+1):n,]))
        }
    }

  wave.muted <- outputw(wave=wave.muted, f=f, format=output)

  if(plot)
    {
      oscillo(wave.muted,f=f,...)
      invisible(wave.muted)
    }
  else 
    {
      return(wave.muted)
    }

}


################################################################################
##                                NDSI
################################################################################

NDSI <- function(
                x,
                anthrophony = 1,
                biophony = 2:8
                )
    
{
    ## INPUT
    if(!is.matrix(x)) stop("'x' should be a two-column matrix obtained with the function sounscapespec().")
    else if(ncol(x) <1 | ncol(x) > 2) stop("'x' should be a two-column matrix obtained with the function sounscapespec().")
    ## INDEX
    alpha <- x[anthrophony, 2]
    beta <- sum(x[biophony, 2])
    res <- (beta - alpha) / (beta + alpha)
    names(res) <- NULL # to remove the colnames 'amplitude' of the soundscapespec result
    return(res)
}


################################################################################
##                                NOISEW
################################################################################

noisew<-function(
                 f,
                 d,
                 type = "unif",
                 listen = FALSE,
                 output = "matrix"
                 )

{
  if(type == "unif") wave <- as.matrix(runif(d*f,min=-1,max=1))
  if(type == "gaussian") wave <- as.matrix(rnorm(d*f))
  wave <- outputw(wave=wave, f=f, format=output)
  if(listen) {listen(wave,f=f)}
  return(wave)
}


################################################################################
##                                NOTEFREQ
################################################################################


notefreq <- function(note,
                     ref=440,
                     octave=3
                     )
  {
    if(is.character(note))
      {
        n <- nchar(note)
        if(n > 2) stop("'note' cannot be a character vector with more than 2 characters")
        if(any(note == c("E#","Fb", "B#", "Cb"))) stop("This note does not exist")
        if(n==2)
          {
            notesplit <- unlist(strsplit(note, split=NULL))
            if(notesplit[2]=="b") names <- c("C","Db","D","Eb","E","F","Gb","G","Ab","A","Bb","B")
            if(notesplit[2]=="#") names <- c("C","C#","D","D#","E","F","F#","G","G#","A","A#","B") 
          }
        else names <- c("C","C#","D","D#","E","F","F#","G","G#","A","A#","B") 
        note <- which(names==note)
      }
    f <- ref*2^((octave-3) + ((note-10)/12))
    return(f)
  }


################################################################################
##                                OCTAVES
################################################################################

octaves <- function(x, below=3, above=3)
  {
    y <- numeric(below)
    z <- numeric(above)

    for(i in 1:below)  {y[i] <- x/(2^i)}
    for(i in 1:above)  {z[i] <- x*2^i}
    res <- c(rev(y), x, z)
    return(res)
  }


################################################################################
##                                OSCILLO
################################################################################

oscillo <- function(
    wave,
    f,
    from = NULL,
    to = NULL,
    scroll = NULL,
    zoom = FALSE,
    k=1,
    j=1,
    cex = NULL,
    labels = TRUE,
    tlab = "Time (s)",
    alab = "Amplitude",            
    byrow = TRUE,
    identify = FALSE,
    nidentify = NULL, 
    plot = TRUE,
    colwave = "black",
    coltitle = "black",
    cextitle = 1.2,
    fonttitle = 2,
    collab = "black",
    cexlab = 1,
    fontlab = 1,
    colline = "black",
    colaxis = "black",
    cexaxis = 1,
    fontaxis = 1,
    coly0 = "lightgrey",
    tcl = 0.5,
    title = FALSE,
    xaxt= "s",
    yaxt= "n",
    type = "l",
    bty = "l"
)

{
    input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; rm(input)

    p<-k*j

    if(!is.null(from) && is.na(from)) {from <- NULL}
    if(!is.null(to) && is.na(to))     {to <- NULL}
    if(is.null(from) && is.null(to))  {a<-0; b<-length(wave); from<-0; to<-length(wave)/f}
    if(is.null(from) && !is.null(to)) {a<-1; b<-round(to*f); from<-0}
    if(!is.null(from) && is.null(to)) {a<-round(from*f); b<-length(wave); to<-length(wave)/f}
    if(!is.null(from) && !is.null(to))
        {
            if(from>to) stop("'from' cannot be superior to 'to'")
            if(from==0) {a<-1} else {a<-round(from*f)}
            b<-round(to*f)
        }
    wave<-as.matrix(wave[a:b,])
    n<-nrow(wave)

    if(plot)
        {
            alim<-max(abs(wave))
                                        # to get a single window view
            if(k==1 & j==1)
                {
                    if(!is.null(scroll))
                        {
                            if(!is.numeric(scroll)) stop("scroll has to a numeric")
                            if(length(scroll)>1) stop("length of scroll cannot be superior to 1")
                            if(zoom) stop("zoom and scroll cannot be used together")
                            if(identify) stop("identify and scroll cannot be used together")
                            step<-round(seq(0,n,length.out=scroll+1))
                            lstep<-length(step)
                            pos<-1:(lstep-1)
                            plot.dynosc<-function(panel)
                                {
                                    with(panel,
                                         {
                                             soscillo(wave = wave, f = f, from = step[pos]/f,
                                                      to=step[pos+1]/f,
                                                      colwave = colwave, collab = collab, tlab=tlab, alab=alab,
                                                      cexlab = cexlab, fontlab = fontlab, colline = colline,
                                                      colaxis = colaxis, cexaxis = cexaxis, fontaxis=fontaxis, coly0 = coly0, bty = bty,
                                                      tickup=max(abs(wave),na.rm=TRUE), ylim=c(-max(abs(wave)),max(abs(wave))))
                                             title(main=pos,col.main=coltitle,cex.main=cextitle,font.main=fonttitle)
                                         }
                                         )
                                    panel
                                }
                            osc.panel <- rpanel::rp.control("Window")
                            rpanel::rp.slider(osc.panel,pos,from=1,to=lstep-1,resolution=1,
                                              title = "Window", action=plot.dynosc)
                        }

                    else
                        {
                            if(zoom)
                                {
                                    par(tcl=0.5, col.axis=colaxis, cex.axis = cexaxis, font.axis=fontaxis, col=colline,las=0)
                                    plot(x=seq(from,to,length.out=n), y=wave,
                                         col=colwave, type=type, cex=cex,
                                         xaxs="i", yaxs="i",
                                         xlab="", ylab="", ylim=c(-alim,alim),
                                         xaxt=xaxt, yaxt=yaxt,
                                         cex.lab=0.8, font.lab=2,
                                         bty=bty
                                         )
                                    if(bty == "l" | bty == "o")
                                        {axis(side=1, col=colline,labels=FALSE)
                                         axis(side=2, at=max(abs(wave),na.rm=TRUE), col=colline,labels=FALSE)}
                                    mtext(tlab,col=collab, font=fontlab, cex=cexlab, side=1,line=3)
                                    mtext(alab,col=collab, font=fontlab, cex=cexlab,side=2,line=2.5)
                                    abline(h=0,col=coly0,lty=2)

                                    cat("choose start and end positions on the wave\n")
                                    if(.Platform$OS.type == "windows") flush.console()
                                    coord<-locator(n=2)
                                    from<-coord$x[1]; c<-from*f-a
                                    to<-coord$x[2]; d<-to*f-a
                                    if(d<c) {c<-d; d<-c}
                                    wave<-as.matrix(wave[c:d,1])
                                    n<-nrow(wave)
                                }

                            op<-par(tcl=tcl, col.axis=colaxis, cex.axis = cexaxis, font.axis=fontaxis, col=colline,las=0)

                            plot(x=seq(from,to,length.out=n), y=wave,
                                 col=colwave, type=type, cex=cex,
                                 xaxs="i", yaxs="i",
                                 xlab="", ylab="", ylim=c(-alim,alim),
                                 xaxt=xaxt, yaxt=yaxt,
                                 cex.lab=0.8, font.lab=2,
                                 bty=bty)

                            if(bty == "l" | bty == "o")
                                {
                                    axis(side=1, col=colline,labels=FALSE)
                                    axis(side=2, at=max(abs(wave),na.rm=TRUE), col=colline,labels=FALSE)
                                }

                            if(labels)
                                {
                                    mtext(tlab,col=collab, font=fontlab,side=1,line=3,cex=cexlab)
                                    mtext(alab,col=collab, font=fontlab, cex=cexlab,side=2,line=3)
                                }

                            abline(h=0,col=coly0,lty=2)

                            if(is.expression(title)) {title <- title}
                            else{
                                 if(title == FALSE) {title <- ""}
                                 else {
                                        if(is.character(title)) {title <- title}
                                        else {title <- paste("Total time =",as.character(round(n/f,3)), "s - f =",as.character(f),"Hz") }
                            }
                        }
                            title(main=title, col.main=coltitle, cex.main=cextitle, font.main=fonttitle)

                            if(identify)
                                {
                                    cat("choose points on the wave\n")
                                    if(.Platform$OS.type == "windows") flush.console()
                                    x<-seq(from=from,to=to,length.out=n)
                                    y<-wave
                                    if(is.null(nidentify)) {nidentify <- length(x)}
                                    id<-identify(x=x, y=y, labels=round(x,3), col="red", n=nidentify, plot=TRUE)
                                    time <- x[id]
                                    abline(v=time, col="red")
                                    amp <- y[id,1]
                                    res <- cbind(time, amp)
                                    return(res)
                                }
                            par(op)
                        }
                }

                                        # to get a multi-window view
            else
                {
                    if(!is.null(scroll)) stop("scroll cannot be used with a multi-frame window")
                    if(zoom) stop ("'zoom' does work with a single-frame window only ('k'=1 and 'j'=1)")
                    if(identify) stop ("'identify' does work with a single-frame window only ('k'=1 and 'j'=1)")
                    x<-n%/%p
                    def.par <- par(no.readonly = TRUE)
                    on.exit(par(def.par))
                    m<-matrix(1:p,k,j,byrow=byrow)
                    layout(m)
                    par(tcl=tcl,oma=c(3,2,2,0.5),
                        mar=rep(0,4)+0.8, mgp=c(0,0.15,0),
                        col.axis=colaxis, cex.axis = cexaxis, font.axis=fontaxis, col=colline, las=0)

                                        # plots the first window
                    wave1<-as.matrix(wave[0:x,]); n1<-nrow(wave1)
                    plot(x=seq(from,from+(x/f),length.out=n1), y=wave1,
                         col=colwave, type=type, cex=cex,
                         xaxs="i", yaxs="i",
                         xlab="", ylab="", ylim=c(-alim,alim),
                         xaxt=xaxt, yaxt=yaxt,
                         bty=bty)
                    axis(side=1, col=colline,labels=FALSE)
                    if(bty == "l" | bty == "o")
                        {axis(side=2, at=max(abs(wave)), col=colline,labels=FALSE)
                         axis(side=1, col=colline,labels=FALSE)}
                    abline(h=0,col=coly0,lty=2)

                                        # title
                    if(is.character(title)) title<-paste(title)
                    if(title == FALSE) {title <- paste("")}
                    else
                        {
                            title<-paste("Window time =",
                                         as.character(round(n/(p*f),3)),"s - Total time =",
                                         as.character(round(n/f,3)), "s - f =",
                                         as.character(f),"Hz")
                        }
                    mtext(paste(title),side=3,line=0.4,col=coltitle,cex=cextitle,font=fonttitle,outer=TRUE)

                                        # X-Y labels
                    if(labels)
                        {
                            mtext(tlab, col=collab, side=1,line=1.5, font=fontlab,cex=cexlab,outer=TRUE)
                            mtext(alab, col=collab, side=2, font=fontlab,cex=cexlab,
                                  line=0.4,outer=TRUE)
                        }

                                        # plots following windows
                    for(i in 1:(p-1))
                        {
                            xx<-((i*n)%/%p)+1
                            yy<-((i+1)*n)%/%p
                            wave2<-as.matrix(wave[xx:yy,]); n2<-nrow(wave2)
                            plot(x=seq(from+(xx/f),from+(yy/f),length.out=n2), y=wave2,
                                 col=colwave, type=type, cex=cex,
                                 xaxs="i", yaxs="i",
                                 xlab="", ylab="", ylim=c(-alim,alim),
                                 xaxt=xaxt, yaxt=yaxt,
                                 bty=bty)

                            if(bty == "l" | bty == "o")
                                {axis(side=2, at = max(abs(wave)), col=colline,labels=FALSE)
                                 axis(side=1, col=colline,labels=FALSE)}
                            abline(h=0,col=coly0,lty=2)
                        }
                }
            invisible(wave)
        }
    else return (wave)
}


################################################################################
##                                OSCILLOST
################################################################################

oscilloST <- function(
                    wave1,
                    wave2 = NULL,
                    f,
                    from = NULL,
                    to = NULL,
                    identify = FALSE,
                    plot = TRUE,
                    colwave1 = "black",
                    colwave2 = "blue",
                    coltitle = "black",
                    collab = "black",
                    cexlab = 1,
                    fontlab = 1,
                    colaxis = "black",
                    cexaxis = 1,
                    coly01 = "grey47",
                    coly02 = "black",
                    title = FALSE,
                    bty = "l"
                    )

{

  if(class(wave1)=="Wave" && wave1@stereo==TRUE) {
      wave2 <- mono(wave1, which="right")
      wave1 <- mono(wave1, which="left")
  }
  
  if(class(wave1)=="audioSample" && nrow(wave1 > 1)) {
      wave2 <- audio::audioSample(wave1[2,], rate=wave1$rate)
      wave1 <- audio::audioSample(wave1[1,], rate=wave1$rate)
  }
  
  input1 <- inputw(wave1, f=f, channel=1) ; wave1<-input1$w ; f<-input1$f ; rm(input1)
  wave2 <- inputw(wave2, f=f, channel=1)$w

  if(plot)
    {
      op<-par(mfrow=c(2,1),oma=c(5,3,2,2),mar=rep(0,4), cex.axis=cexaxis)
      
      oscillo(wave=wave1,f=f,
              from=from,to=to,zoom=FALSE,labels=FALSE,xaxt="n",
              colaxis=colaxis,colwave=colwave1,coly0=coly01,
              bty=bty)
      
      oscillo(wave=wave2,f=f,
              from=from,to=to,identify=identify,zoom=FALSE,labels=FALSE,
              colaxis=colaxis,colwave=colwave2,coly0=coly02,
              bty=bty)

      mtext("Time (s)",col=collab,font=fontlab,cex=cexlab,side=1,line=2.8,outer=TRUE)
      mtext("Amplitude",col=collab,font=fontlab,cex=cexlab,side=2,line=1.5,outer=TRUE)

      par(op)
      invisible(cbind(wave1,wave2))
    }

  else return(cbind(wave1,wave2))
}


################################################################################
##                                PASTEW
################################################################################

pastew <- function (wave1,
                    wave2,
                    f,
                    at = "end",
                    join = FALSE,
                    tjunction = 0,
                    choose = FALSE,
                    plot = FALSE,
                    marks = TRUE,
                    output = "matrix", 
                    ...
                    )
    
{
    input1 <- inputw(wave = wave1, f = f)
    wave1 <- input1$w
    f <- input1$f
    bit <- input1$bit
    rm(input1)
    wave2 <- inputw(wave = wave2, f = f)$w
    n <- nrow(wave2)
    if (choose) {
        cat("choose position on the wave\n")
        if (.Platform$OS.type == "windows") 
            flush.console()
        oscillo(wave2, f = f)
        coord <- locator(n = 1)
        at <- coord$x[1]
        abline(v = at, col = 2, lty = 2)
    }
    else {switch(at, start = {at <- 0}, middle = {at <- n/(2 * f)}, end = {at <- n/f})}
    pos <- round(at * f)
    if (at == 0) wave2a <- NULL
    else wave2a <- wave2[c(1:pos), 1]
    if (at == n/f) wave2b <- NULL
    else wave2b <- wave2[c(pos:n), 1]
    if (join) wave3 <- c(wave2a[-length(wave2a)], wave1, wave2b)
    else wave3 <- c(wave2a, wave1, wave2b)
    junction <- round(tjunction * f)
    if (junction != 0) {
        p <- length(wave2a)
        q <- length(wave2b)
        r <- floor(junction/2)
        s <- length(wave1)
        if (junction > s) {
            stop("t_junction must be largely shorter than wave1")
        }
        if (!is.null(wave2a)) {
            if (r > p) {stop("t_junction/2 must be largely shorter than the first part of wave2")}
            else {wave3[(p - r + 1):(p + junction - r)] <- seq(wave3[p - r + 1], wave3[p + junction - r], length.out = junction)}
        }
        if (!is.null(wave2b)) {
            if (junction - r > q) {stop("t_junction/2 must be largely shorter than the second part of wave2")}
            else {wave3[(p + s - r + 1):(p + s + junction - r)] <- seq(wave3[p + s - r + 1], wave3[p + s + junction - r], length.out = junction)}
        }
    }
    wave3 <- outputw(wave = wave3, f = f, format = output)
    if (plot) {
        def.par <- par(no.readonly = TRUE)
        on.exit(par(def.par))
        par(mfrow = c(3, 1), oma = c(0, 0.1, 0, 0))
        oscillo(wave1, f = f, k = 1, j = 1)
        title(main = "signal to be pasted")
        oscillo(wave2, f = f, k = 1, j = 1)
        title(main = "signal to be completed")
        oscillo(wave3, f = f, k = 1, j = 1)
        title(main = "resulting signal")
        if (marks) {
            abline(v = at, col = "red", lty = 2)
            abline(v = at + (nrow(wave1))/f, col = "red", lty = 2)
        }
        invisible(wave3)
    }
    else {return(wave3)}
}



################################################################################
##                                PHASEPLOT
################################################################################


phaseplot <- function(wave,
                      f,
                      dim = 3,
                      plot = TRUE,
                      type = "l",
                      xlab = "1st derivative",
                      ylab = "2nd derivative",
                      zlab = "3rd derivative",
                      ...)
{
  ## Error messages
  if(dim <= 1 | dim > 3) stop("'dim' has to be set to 2 or 3.")
  ## Input
  wave <- inputw(wave = wave, f = f)$w
  wave <- wave/max(wave) 
  ## Derivatives  
  z <- matrix(nrow = nrow(wave)-3, ncol = 3)
  z[,1] <- diff(wave, 1, 1)[-(1:2)]
  z[,2] <- diff(wave, 1, 2)[-1]
  if(dim==3) z[,3] <- diff(wave, 1, 3)
  ## Plot
  if(plot){
    if(dim==3) rgl::plot3d(z[,1], z[,2], z[,3], type=type, xlab=xlab, ylab=ylab, zlab=zlab,...)
    else plot(z[,1], z[,2], type=type, xlab=xlab, ylab=ylab,...)
    invisible(z)
  }
  else return(z)
}


################################################################################
##                                PLAYLIST
################################################################################

playlist <- function(directory, sample=FALSE, loop=1)
  {
    files <- dir(directory, pattern={".wav|.mp3|.ogg|.aiff"})
    n <- length(files)
    if(isTRUE(sample)) {files <- sample(files, size=n, replace=FALSE)}
    i <- 1
    while(i <= loop){
      for(j in 1:n) {play(shQuote(paste(directory, files[j], sep="")))}
      i <- i+1
    }
  }

################################################################################
##                                PULSEW
################################################################################

pulsew<-function(
                dbefore,
                dpulse,
                dafter,
                f,
                plot = FALSE,
                output = "matrix",
                ...
                )

{
  wave<-c(rep(0,dbefore*f),rep(1,dpulse*f),rep(0,dafter*f))

  wave <- outputw(wave=wave, f=f, format=output)

  if(plot)
    {
      oscillo(wave,f=f,...)
      invisible(wave)
    }
  else {return(wave)}
}




################################################################################
##                                Q
################################################################################

Q <- function(spec,
              f = NULL,
              level = -3,
              mel = FALSE,
              plot = TRUE,
              colval = "red",
              cexval = 1,
              fontval = 1,
              flab = NULL,
              alab = "Relative amplitude (dB)",
              type = "l",
              ...)

{ if (!is.null(f) & mel) {f <- 2*mel(f/2)}
  if(is.null(f))
      {
          if(is.vector(spec)) stop("'f' is missing")  
          else if(is.matrix(spec)) f<-spec[nrow(spec),1]*2000*nrow(spec)/(nrow(spec)-1)
      }

  if(is.matrix(spec)) spec <- spec[,2]

  range<-c(f/2000/length(spec),f/2000)

  if(max(spec) == 1) stop ("data must be in dB")
  if(which.max(spec) == 1) stop ("maximal peak cannot be the first value of the spectrum") 

  n0<-length(spec)

  spec1<-approx(spec,n=102400)$y
  n1<-length(spec1)
  level2<-round(max(spec1),1)+level

  f0<-which.max(spec1)
  f0khz<-((f0/n1)*(range[2]-range[1]))+range[1]
  specA<-as.matrix(spec1[1:f0])
  nA<-nrow(specA)
  specB<-as.matrix(spec1[f0:length(spec1)])
  f1<-which(round(specA,1) == level2)
  f1khz<-((f1[length(f1)]/n1)*(range[2]-range[1]))+range[1]
  f2<-which(round(specB,1) == level2)+(nA-1)
  f2khz<-((f2[1]/n1)*(range[2]-range[1]))+range[1]

  Q <- f0/(f2[1]-f1[length(f1)])

  results <- list(Q=Q, dfreq=f0khz, fmin=f1khz, fmax=f2khz, bdw=f2khz-f1khz)

  ## plot based on original data (=> spectrum)
  if(plot)
      { if(is.null(flab)) {if (mel) flab <- "Frequency (kmel)" else flab <- "Frequency (kHz)"}
        x<-seq(range[1],range[2],length.out=n0)
        plot(x=x,y=spec,xlab=flab,ylab=alab,type=type,...)
        arrows(f1khz,level2,f2khz,level2,length=0.1,col=colval,code=3,angle=15)
        text(paste("Q =",as.character(round(Q,2))),x=f2khz,y=level2,pos=4,
             col=colval, cex=cexval, font=fontval)
        invisible(results)
    }
  return(results)
}


################################################################################
##                                REPW
################################################################################

repw<-function(
               wave,
               f,
               times = 2,
               join = FALSE,
               plot = FALSE,
               output = "matrix",
               ...
               )

{
  input <- inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; bit <- input$bit ; rm(input)
  n <- nrow(wave)
  if(join) wave1 <- rep(wave[-n],times=times) else wave1 <- rep(wave,times=times)
  wave1 <- outputw(wave=wave1, f=f, format=output)

  if(plot)
    {
      oscillo(wave=wave1,f=f,...)
      invisible(wave1)
    }
  else {return(wave1)}
}


################################################################################
##                                REVW
################################################################################

revw<-function(
               wave,
               f,
               env = TRUE,
               ifreq = TRUE,
               plot = FALSE,
               output = "matrix",
               ...)

{
  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; bit <- input$bit ; rm(input)

  if(env == FALSE & ifreq == FALSE) stop ("Both arguments 'env' and 'ifreq' cannot be set to FALSE.")

  if(env & ifreq) {wave2<-as.matrix(rev(wave[,1]))}
  else
    {
      wave.e<-env(wave[,1],f=f,plot=FALSE)
      wave.p<-ifreq(wave[,1],f=f,plot=FALSE)$p[,2]
      if(env & ifreq== FALSE) {wave2<-as.matrix(rev(wave.e)*cos(wave.p))}
      if(env==FALSE & ifreq) {wave2<-as.matrix(wave.e*cos(rev(wave.p)))}
    }

  wave2 <- outputw(wave=wave2, f=f, format=output)
  
  if(plot)
    {
      oscillo(wave=wave2,f=f,...)
      invisible(wave2)
    }
  else {return(wave2)}
}



################################################################################
##                                RESAMP
################################################################################


resamp<-function(
                 wave,
                 f,
                 g,
                 output = "matrix"
                 )

{
  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; bit <- input$bit ; rm(input)$w

  n<-nrow(wave)
  if(g==f) stop ("'f' and 'g' must be different")
  if(g<f)  {r<-f/g; wave1<-wave[seq(1,n,by=r),1]}
  if(g>f)  {s<-(n*g)/f; wave1<-approx(wave,n=s)$y}

  wave1 <- outputw(wave=wave1, f=g, format=output)
  
  return(wave1)
}


################################################################################
##                                RMAM
################################################################################

rmam <- function(
               wave,
               f,
               plot = FALSE,
               listen = FALSE,
               output = "matrix",
               ...
               )

{
  input <- inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; bit <- input$bit ; rm(input)
  wave <- wave/Mod(hilbert(wave,f=f))
  wave <- outputw(wave=wave, f=f, format=output)

  if(plot)
    {
      oscillo(wave=wave,f=f,...)
      if(listen) {listen(wave,f=f)}
      invisible(wave)
    }
  else
    {
      if(listen) {listen(wave,f=f)}
      return(wave)
    }
}


################################################################################
##                                RMOFFSET
################################################################################

rmoffset <-function(
                   wave,
                   f,
                   FUN = mean,
                   plot = FALSE,
                   output = "matrix",
                   ...
                   )

{
  input <- inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; bit <- input$bit ; rm(input)
  m <- apply(wave, MARGIN=2, FUN=FUN)
  wave <- wave-m

  wave <- outputw(wave=wave, f=f, format=output)
  
  if(plot)
    {
      oscillo(wave=wave,f=f,...)
      invisible(wave)
    }
  else return(wave)
}


################################################################################
##                               RMNOISE
################################################################################


rmnoise<-function(
                  wave,
                  f,
                  output = "matrix",
                  ...                  
                  )

  {
    input <- inputw(wave=wave, f=f); wave<-input$w; f<-input$f; bit <- input$bit ; rm(input)
    wave2 <- smooth.spline(wave, all.knots=TRUE,...)$y  
    wave2 <- outputw(wave=wave2, f=f, format=output)
    return(wave2)
  }




################################################################################
##                               RMS
################################################################################

rms <- function(
                x,
                ...
                )
  {
    x <- sqrt(mean(x^2, ...))
    return(x)
  }


################################################################################
##                               ROUGHNESS
################################################################################

roughness <- function(x, std=FALSE){
  if(std) x <- x/max(x)
  deriv2 <- diff(x, 1, 2)
  roughness <- sum(deriv2^2, na.rm=TRUE)
  return(roughness)
} 
  

################################################################################
##                               RUGO
################################################################################

rugo <- function(
                x,
                ...
                )
  {
    n <- length(x)
    y <- numeric(n-1)
    for(i in 1:(n-1)) {y[i] <- (x[i+1]-x[i])^2 }
    rug <- sqrt(mean(y, ...))
    return(rug)
  }


################################################################################
##                               SAVEWAV
################################################################################

savewav <- function(
                  wave,
                  f = NULL,
                  filename = NULL,
                  rescale = NULL
                  )

{
  ## ERROR MESSAGE
  if(!is.null(rescale))
    {
      if(rescale[1] >= 0) stop("The first value of 'rescale' should not be >=0")
      if(rescale[2] <= 0) stop("The first value of 'rescale' should not be <=0")
    }
  ## FILENAME  
  if(is.null(filename)) filename <- paste(as.character(deparse(substitute(wave))),".wav",sep="")
  
  ## INPUT
  input <- inputw(wave=wave, f=f) ; wave <- input$w ; f <- input$f ; bit <- input$bit
  rm(input)

  ## OUTPUT
  if(!is.null(rescale))
       {
        wave <- rescale(wave, lower=rescale[1], upper=rescale[2])
        wave <- Wave(left=wave, samp.rate=f, bit=bit)
       }
  else {
        wave <- Wave(left=wave, samp.rate=f, bit=bit)
        max <- max(wave@left)
        if(max <= 1) {level <- max} else {level <- 1}
        wave <- normalize(wave, unit=as.character(bit), level=level)
       }
  
  ## WRITING TO FILE
  writeWave(wave, filename=filename)
}


################################################################################
##                               SAX
################################################################################


SAX <- function (
                 x,
                 alphabet_size,
                 PAA_number,
                 breakpoints = "gaussian",
                 collapse = NULL
                )
    
{
    ## INPUT
    if (class(x) == "Wave") {x <- soundscapespec(x)[,2]}
    
    pos <- function(t,v) {which.max(v[v<=t])}

    if (alphabet_size > PAA_number)
        {warning('alphabet_size must preferably be smaller than PAA_number', call. = FALSE)}

    l <- length(x)

    PAA <- array(0, PAA_number)
    for (i in 1:PAA_number)
        {PAA[i] <- mean(x[round((i-1)*l/PAA_number+1):round(i*l/PAA_number)])}

    if (class(breakpoints) == 'character')
        {
            if (breakpoints == 'gaussian') {q <- (qnorm(seq(0,1,1/alphabet_size)))*sd(x)+mean(x)}
            else if (breakpoints == 'quantiles') {q <- as.numeric(quantile(x, probs = seq(0, 1, 1/alphabet_size)))}
            
            else {stop ('this method is not implemented')}
        }

    else if (class(breakpoints) == 'numeric') {
        if(length(breakpoints) != alphabet_size-1) stop('The number of breakpoints should not be different from the (alphabet_size - 1)')
        q <- c(-Inf,breakpoints, Inf)}
    else {stop ('breakpoints class must be character or numeric')}

    res <- paste(letters[sapply(PAA, pos, v = q)], collapse = collapse)
    
    return(res)
}


################################################################################
##                               SDDB
################################################################################

sddB <- function(
                 x,
                 level = "IL"
                 ) 

{
  if(level == "IL") {a <- 10} else {a <-20}
  x <- 10^(x/a)
  res <- a * (sd(x) / (mean(x)*log(10)))
  return(res)
}


################################################################################
##                               SEEDATA
################################################################################

seedata <- function(
                    data,
                    na.rm = FALSE,
                    col= "grey"
                    )

{
                                        # na
  if(na.rm) data<-na.omit(data)

                                        # layout
  def.par <- par(no.readonly = TRUE)
  on.exit(par(def.par))
  m<-layout(matrix(c(1,1,2,2,3,4,5,5), nrow = 4, ncol = 2, byrow=TRUE), heights=c(3,2,3,2))
  par(mar=c(4.5,4,3,5)+0.1,oma=c(0,0.5,0,0))

                                        # histogram and density
  hist(data, main=NULL, col=col, xlab="", prob=TRUE)
  lines(density(data), lwd=1.75)
  x <- NULL
  curve(dnorm(x, mean(data), sd(data)), add=TRUE, col="red")
  title(main="Histogram and density")


                                        # stripchart
  stripchart(data)
  title(main="Stripchart")

                                        # boxplot
  boxplot(data, col=col, lwd=0.5)
  title(main="Boxplot")

                                        # qq plot
  qqnorm(data, main="Quantile-quantile plot")
  qqline(data, col="red")

  
                                        # normality tests
  n <- length(data)

  shapi <- shapiro.test(data)
  shapiW <- round(shapi$statistic, 3)
  shapiP <- round(shapi$p.value, 3)

  par(mar=rep(1,4))
  plot.new()
  text(x=0.5, y=0.5,
       paste(
             "NORMALITY TEST", 
             "\nSample size: n=", n,
             "\nShapiro-Wilk test:  W=", as.character(shapiW), ", p=",as.character(shapiP),
             sep=""
             )
       )
}


################################################################################
##                               SETENV
################################################################################

setenv<-function(
                 wave1,
                 wave2,
                 f,
                 envt="hil",
                 msmooth = NULL,
                 ksmooth = NULL,
                 plot = FALSE,
                 listen = FALSE,
                 output = "matrix",
                 ...
                 )

{
  input1<-inputw(wave=wave1,f=f) ; wave1<-input1$w ; f<-input1$f ; bit <- input1$bit ; rm(input1)
  wave2<-inputw(wave=wave2,f=f)$w

  wave1<-rmoffset(wave1,f=f)
  wave1<-rmam(wave1,f=f)
  wave1<-wave1/max(abs(wave1))

  wave2<-rmoffset(wave2,f=f)
  wave2.env<-env(wave2,f=f,envt=envt,msmooth=msmooth,ksmooth=ksmooth,plot=FALSE)
  wave2.env<-approx(wave2.env, n=nrow(wave1))$y

  wave3<-wave1*wave2.env
  wave3<-wave3/max(abs(wave3))

  wave3 <- outputw(wave=wave3, f=f, format=output)

  if(plot)
    {
      oscillo(wave=wave3,f=f,...)
      if(listen) {listen(wave3,f=f)}
      invisible(wave3)
    }
  else
    {
      if(listen) {listen(wave3,f=f)}
      return(wave3)
    }
}


################################################################################
##                               SFM
################################################################################

sfm <- function(spec)
{
    if(is.matrix(spec)) spec <- spec[, 2]
    if(any(spec<0)) stop("Data do not have to be in dB")
    if(sum(spec)==0) flat<-NA
                                        # undersample spec if too long because prod(spec) tends towards zero
    if(length(spec) > 4000)
        {
            step<-seq(1,length(spec),by=round(length(spec)/256))
            spec<-spec[step]
        }
    spec<-ifelse(spec==0,yes=1e-5,no=spec)
    n<-length(spec)
    geo<-prod(spec^(1/n))
    ari<-mean(spec)
    flat<-geo/ari
    return(flat)
}

################################################################################
##                               SH
################################################################################

sh <- function(    
               spec,
               alpha = 'shannon'
               )

{

  if(is.matrix(spec) && ncol(spec)==2) spec<-spec[,2]
  options(warn=-1)  # to ignore warning messages
  if(any(is.na(spec))) {options(warn=-1) ; z <- NA ; options(warn=0)}  # options(warn) to ignore temporarely warning messages
  else
    {
      options(warn=0)  # to authorize warning messages
      if(any(spec<0, na.rm=TRUE)) stop("Data do not have to be in dB.")
      if(sum(spec)==0) {warning("Caution! This is a null spectrum. The spectral entropy is null!", call.=FALSE) ; z <- 0}
      else
        {
          N<-length(spec)
          spec[spec==0]<-1e-7
          spec<-spec/sum(spec)  # PMF
          if(alpha == 'shannon') {z <- -sum(spec*log(spec))/log(N)}
	  else if (alpha == 'simpson') {z <- 1 - sum(spec^2)}
          else {
            if(alpha < 0) stop ("'alpha' cannot be negative.")
            if(alpha == 1) stop ("'alpha' cannot be set to 1.")
            z <- (1/(1-alpha))*log(sum(spec^alpha))/log(N)
          }
        }
    }
  return(z)
}


################################################################################
##                                SIMSPEC
################################################################################

simspec <- function(spec1,
                    spec2,
                    f = NULL,
                    mel = FALSE,
                    norm = FALSE,
                    PMF = FALSE,
                    plot = FALSE,
                    type = "l",
                    lty = c(1,2,3),
                    col = c(2,4,1),
                    flab = NULL,
                    alab = "Amplitude (percentage)",
                    flim = NULL,
                    alim = NULL,
                    title = TRUE,
                    legend = TRUE,
                    ...
                    )

{
    leg<-c(as.character(deparse(substitute(spec1))),as.character(deparse(substitute(spec2))))

    if(norm & PMF) stop("'norm' and 'PMF' should not be both set to TRUE")

    if(is.matrix(spec1) && ncol(spec1)==2) s1<-spec1[,2] else s1 <- spec1
    if(is.matrix(spec2) && ncol(spec2)==2) s2<-spec2[,2] else s2 <- spec2

    n1<-length(s1)
    n2<-length(s2)

    if(n1 != n2) stop("spec1 and spec2 must have the same length")

    if(any(is.na(s1)) | any(is.na(s2)))
        {
            S <- NA
            warning("The data set contains 'NA' values. The returned values have been set to NA.", call.=FALSE)
            return(S)
        } 
    else
	{

            if(any(s1 < 0) | any(s2 < 0)) stop("spectra (spec 1 and/or spec 2) do not have to be in dB")
            if(sum(s1)==0) {warning(paste("Caution!, spec1 is a null spectrum"), call.=FALSE)}
            if(sum(s2)==0) {warning(paste("Caution!, spec2 is a null spectrum"), call.=FALSE)}

            if (norm) {s1<-s1/max(s1)
                       s2<-s2/max(s2)}

            if (PMF) {s1<-s1/sum(s1)
                      s2<-s2/sum(s2)}

            s1[s1==0] <- 1e-7
            s2[s2==0] <- 1e-7
            S1<-100*(pmin(s1,s2)/pmax(s1,s2))
            S<-sum(S1)/n1

            if(plot)
                {  
                    if (!is.null(f) & mel) {f <- 2*mel(f/2)}
                    if(is.null(f))
                        {
                            if(is.vector(spec1) & is.vector(spec2)) stop("'f' is missing")  
                            else
                                {
                                    if(is.matrix(spec1)) f<-spec1[nrow(spec1),1]*2000*nrow(spec1)/(nrow(spec1)-1)
                                    else if(is.matrix(spec2)) f<-spec2[nrow(spec2),1]*2000*nrow(spec2)/(nrow(spec2)-1)
                                }
                        }
                    if(is.null(alim)) alim<-c(0, 100)
                    if(is.null(flim)) flim<-c(0,f/2000)
                    if(is.null(flab)) {if (mel) flab <- "Frequency (kmel)" else flab <- "Frequency (kHz)"}
                    x<-seq(0,f/2000*(n1-1)/n1,length.out=n1)
                    plot(x=x, y=s1*100/max(s1,s2), type=type, lty=lty[1], col=col[1],
                         xlim=flim, xaxs="i", xlab=flab,
                         ylim=alim, yaxs="i", ylab=alab,...)
                    lines(x=x, y=s2*100/max(s1,s2), type=type, lty=lty[2], col=col[2])
                    lines(x=x, y=S1, type=type, lty=lty[3], col=col[3])
                    if(legend) legend("topleft", col=c(col[1],col[2]),lty=c(lty[1],lty[2]),legend=leg, bty="n")
                    if(title) title(main=paste("S = ", round(S, 3)))
                    invisible(S)
                }

            return(S)
        }
}


################################################################################
##                                SMOOTHW
################################################################################

smoothw <- function(
                    wave,
                    f,
                    wl,
                    padding=TRUE,
                    output="matrix"
                    )
{
   input <- inputw(wave = wave, f = f)
   wave <- input$w
   f <- input$f
   rm(input)
   wave <- sumsmooth(wave, wl=wl, padding=padding)
   wave <- rmoffset(wave, f=f, output=output)
   return(wave)
}


################################################################################
##                                SONGMETER
################################################################################

songmeter <- function(
                      x             ## a character vector (not a Wave object)
                     )
    
    {
        ## ERROR MESSAGES
        if(!is.character(x)) stop("'x' should be of mode character.")

        ## PREPARE RESULTS
        options(stringsAsFactors = FALSE)
        res <- data.frame(model = character(0), prefix = character(0), mic = character(0),
                          year = numeric(0), month = numeric(0), day = numeric(0),
                          hour = numeric(0), min = numeric(0), sec = numeric(0),
                          time = numeric(0), geo = logical(0))

        N <- length(x)  ## number of files

        ## LOOP
        for (i in 1: N)
            {
                tmp <- x[i]
                n <- nchar(tmp)
                extension <- substr(tmp, start = n-3, stop = n)
                if(any(extension != ".wav" & extension != ".wac"))
                    {
                        warning(paste("File '", tmp, "' is neither a '.wac' nor a '.wav' file", sep=""))
                    }
                else
                    {
                        f <- unlist(strsplit(tmp, extension))
                        if(any(unlist(strsplit(f, "")) == "$")) {
                            geo <- TRUE
                            f <- sub(pattern="$", replacement="_", x=f, fixed=TRUE)
                        }
                        else {geo <- FALSE}
                        
                        items <- unlist(strsplit(f, "_", fixed=TRUE))
                        
                        if (length(items) == 3)
                            {
                                model <- "SM2"
                                prefix <- items[1]
                                mic <- NA
                                year <- substr(items[2], start = 1, stop = 4)
                                month <- substr(items[2], start = 5, stop = 6)
                                day <- substr(items[2], start = 7, stop = 8)
                                hour <- substr(items[3], start = 1, stop = 2)
                                min <- substr(items[3], start = 3, stop = 4)
                                sec <- substr(items[3], start = 5, stop = 6)
                                time <- strptime(paste(items[2], items[3]), "%Y%m%d%H%M%S")
                                geo <- NA
                                res <- rbind(res, data.frame(model, prefix, mic, year, month, day, hour, min, sec, time, geo))
                            }
                        else if (length(items)==4)
                            {
                                model <- "SM3"
                                prefix <- items[1]
                                mic <- items[2]
                                mic[mic=="-0-"] <- "monoL"
                                mic[mic=="-1-"] <- "monoR"
                                mic[mic=="0+1"] <- "stereo"
                                year <- substr(items[3], start = 1, stop = 4)
                                month <- substr(items[3], start = 5, stop = 6)
                                day <- substr(items[3], start = 7, stop = 8)
                                hour <- substr(items[4], start = 1, stop = 2)
                                min <- substr(items[4], start = 3, stop = 4)
                                sec <- substr(items[4], start = 5, stop = 6)
                                time <- strptime(paste(items[3], items[4]), "%Y%m%d%H%M%S")
                                res <- rbind(res, data.frame(model, prefix, mic, year, month, day, hour, min, sec, time, geo))
                            }
                        else
                            {
                                warning(paste("The name of the file '", tmp, "' does not seem to be correct.", sep=""))
                            }
                    }
            }
        res$year <- as.numeric(res$year)
        res$month <- as.numeric(res$month)
        res$day <- as.numeric(res$day)
        res$hour <- as.numeric(res$hour)
        res$min <- as.numeric(res$min)
        res$sec <- as.numeric(res$sec)
        
        return(res)

        options(stringsAsFactors = TRUE)
            }


################################################################################
##                                SOUNDSCAPESPEC
################################################################################

soundscapespec <- function(wave, 
                           f,    
                           wl = 1024,
                           wn = "hamming", 
                           ovlp = 50,
                           plot = TRUE, 
                           xlab = "Frequency (kHz)",
                           ylim = c(0,1),
                           ...     
                           )
    
    {
        ## INPUT
        input <- inputw(wave=wave,f=f) ; wave <- input$w ; f <- input$f ; rm(input)
        ## WELCH SPEC
        n <- nrow(wave)   ## number of samples in wave
        step <- seq(1, n - wl, wl - (ovlp * wl/100))   ## positions of the sliding window for the STFT
        n_recs <- length(step) # EPK: Number of records processed
        x <- seq(f/wl, (f/2) - (f/wl), length.out = wl%/%2)/1000  ## frequency values of the spectrum       
        y <- stft(wave = wave, f = f, wl = wl, zp = 0, step = step, wn = wn, fftw = FALSE, scale = FALSE) ## STFT with the parameters set above
        y <- y^2     # EPK: square of the summed spectrum
        y <- apply(y, MARGIN = 1, FUN = sum)     ## sum of successive spectra
        # y <- y^2     # square of the summed spectrum
        # y <- y*2 / (sum(y)*f)  # normalisation as indicated in an email THE PROBLEM MIGHT BE HERE
        y <- y / (n_recs*f) # EPK: normalize for number of records and frequency
        y <- y*2            # EPK: If half the reflected spectra then *2 
        ## frequency bining 
        freq <- trunc(x)  ## truncate frequency digits to get only 0, 1 etc values
        freq <- freq[freq!=0]  ## eliminate frequency values between 0 and 1 kHz
        spec <- y[trunc(x)!=0]  ## eliminate amplitude values between 0 and 1 kHz
        ## 1 kHz frequency bins
        bin  <- unique(freq) 
        bin.length <- numeric(length(bin))
        for(i in bin) {bin.length[i] <- length(freq[freq==i])}  ## width of each frequency bin
        ## remove small frequency bin at the extreme right of the spectrum
        bin.width <- 1000 * wl / f
        bin.out <- which(bin.length < trunc(bin.width))
        if(length(bin.out)!=0) {bin.new <- bin[-bin.out]} else {bin.new <- bin}
        ## compute the energy of each frequency bin
        val <- numeric(length(bin.new))
        for (i in bin.new) val[i] <- sum(spec[freq==bin.new[i]])

        vlen <- sqrt(sum(val^2)) # EPK: Compute the vector length
        val <- val / vlen        # EPK: Normalize to the range [0,1]

        ## RETURN
        res <- cbind(frequency = bin.new, amplitude = val, deparse.level=0)
        if(plot)
            {
            barplot(val, names=as.character(bin.new), xlab=xlab, ylim=ylim, ...)    
            invisible(res)    
            }
        else return(res)
    }


################################################################################
##                                SOX
################################################################################

sox <- function(command, exename = NULL, path2exe = NULL)
  
  {

    if(.Platform$OS.type == "unix")
      {
        if(missing(exename)) exename <- "sox"
        if(missing(path2exe)) {exe <- exename} else {exe <- paste(path2exe, exename ,sep="/")}
        system(paste(exename, command),ignore.stderr = TRUE)
      }

   if(.Platform$OS.type == "windows")
     {
      if(missing(exename)) exename <- "sox.exe"
      if(missing(path2exe)) {exe <- paste("c:\\sox-14-4-0\\", exename, sep="")} else {exe <- paste(path2exe,exename,sep="/")} 
      e <- system(paste(exe, command, sep=" "), ignore.stderr = TRUE) 
     }

  }


################################################################################
##                                SPEC
################################################################################

spec <- function(
               wave,
               f,
               wl = 512,
               wn = "hanning",
               fftw = FALSE,
               norm = TRUE,
               PSD = FALSE,
               PMF = FALSE,
               dB = NULL,
               dBref = NULL,
               at = NULL,
               from = NULL,
               to = NULL,
               identify = FALSE,
               col = "black",
               cex = 1,
               plot = 1,
               flab = "Frequency (kHz)",
               alab = "Amplitude",
               flim = NULL,
               alim = NULL,
               type ="l",
               ...)

{
                                        # STOP MESSAGES
  if(!isTRUE(norm) & PMF) stop ("'PMF' can be computed only if 'norm' is TRUE")
  if(!isTRUE(norm) & !is.null(dB)) stop ("dB are computed on normalised spectra only, 'norm' should be turned to TRUE")
  if(!is.null(dB) & PMF) stop("PMF cannot be in dB")
  if(!is.null(dB) & PSD) stop("PSD cannot be in dB")
  if(is.null(dB) & !is.null(dBref)) stop("'dB' cannot be NULL  when 'dBref' is not NULL")
  if(is.logical(dB)) stop("'dB' is no more a logical. Please see the documentation: help(spec).")
  if(!is.null(dB) && all(dB!=c("max0","A","B","C","D")))
    stop("'dB' has to be one of the following character strings: 'max0', 'A', 'B', 'C' or 'D'")
  if(!is.null(wl) & wl%%2 == 1) stop("'wl' has to be an even number.")

                                        # INPUT
  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; rm(input)

                                        # FROM-TO SELECTION
  if(!is.null(from)|!is.null(to))
    {
      if(is.null(from) && !is.null(to)) {a<-1; b<-round(to*f)}
      if(!is.null(from) && is.null(to)) {a<-round(from*f); b<-length(wave)}
      if(!is.null(from) && !is.null(to))
        {
          if(from>to) stop("'from' cannot be superior to 'to'")
          if(from==0) {a<-1} else a<-round(from*f)
          b<-round(to*f)
        }
      wl <- (b - a) + 1    ## fixed by Arvind Sowmyan, 2014-03-10
      wave<-as.matrix(wave[a:b,])
    }

                                        # AT SELECTION
  if(!is.null(at))
    {
      c<-round(at*f)
      wl2<-wl%/%2
      wave<-as.matrix(wave[(c-wl2):(c+wl2),])
    }

                                        # FFT
  n<-nrow(wave)
  W<-ftwindow(n,wn=wn)
  wave<-wave*W
  if(fftw == FALSE) {y<-Mod(fft(wave[,1]))}
  else {
    p <- fftw::planFFT(n)
    y <- Mod(fftw::FFT(wave[,1], plan=p))
  }
  y <- y[1:(n%/%2)]

                                        # NORM OPTION
  if(norm) {y <- y /max(y)}
     
  y <- ifelse(y==0, yes=1e-6, no=y)  # replaces 0 values in spectra that cannot be processed by log10()

                                        # FREQUENCY DATA
  x<-seq(0, (f/2)-(f/wl), length.out= n%/%2) / 1000

                                        # PSD and PMF OPTIONS
  if(PSD) y<-y^2
  if(PMF) y<-y/sum(y)
 
                                        # DB
  if(!is.null(dB))
    {
      if(is.null(dBref)) {y <- 20*log10(y)} else {y <- 20*log10(y/dBref)}
      if(dB!="max0")
          {
            if(dB == "A") y <- dBweight(x*1000, dBref = y)$A 
            if(dB == "B") y <- dBweight(x*1000, dBref = y)$B 
            if(dB == "C") y <- dBweight(x*1000, dBref = y)$C 
            if(dB == "D") y <- dBweight(x*1000, dBref = y)$D
          }
    }
  
                                        # LIMITS OF THE AMPLITUDE AXIS
  if(is.null(alim))
    {
      if(is.null(dB)) {alim<-c(0,1.1)} else {alim <- c(min(y, na.rm=TRUE), max(y, na.rm=TRUE) + 20)}
      if(PMF | !isTRUE(norm)) alim<-c(0,max(y, na.rm=TRUE))
    }

                                        # HORIZONTAL PLOT
  if(plot == 1)
    {
      if(!is.null(dB))
        {
          plot(x=x, y=y,
               xaxs = "i", xlab = flab, xlim = flim,
               yaxs = "i", yaxt = "s", ylab = alab, ylim = alim,
               col = col, cex = cex,
               type = type, las = 1,
               ...)
        }
      else
        {         
          if(isTRUE(norm))
            {
              yaxt<-"n"
              ylab<-alab
              if(isTRUE(PMF)) {yaxt="s"}
            }
          else
            {
              yaxt<-"s"
              ylab<-" "
            }

          plot(x=x,y=y,
               xaxs="i", xlab=flab, xlim = flim,
               yaxs="i", yaxt=yaxt, ylab = ylab, ylim=alim,
               col = col, cex = cex,
               type = type, las = 1,
               ...)
        }

      if(identify)
        {
          cat("Choose points on the spectrum\n")
          if(.Platform$OS.type == "windows") flush.console()
          id<-identify(x=x,y=y,labels=round(x,2),tolerance=0.15,col="red")
          id.freq<-x[id]
          id.amp<-y[id]
          coord<-list(freq = id.freq ,amp = id.amp)
          return(coord)
        }
    }

                                        # VERTICAL PLOT
  if(plot == 2)
    {
      if(!is.null(dB))
        {
          plot(x=y,y=x,
               xaxs = "i", xlab = alab, xlim = alim,
               yaxs = "i", yaxt = "s", ylab = flab, ylim = flim,
               col = col, cex = cex,
               type = type, las = 1,
               ...)
        }
      else
        {
          if(isTRUE(norm))
            {
              xaxt <- "n"
              xlab <- alab
              if(isTRUE(PMF)) {xaxt = "s"}
            }
          else
            {
              xaxt<-"s"
              xlab<-" "
            }
          
          plot(x=y,y=x,
               xaxs = "i", xaxt = xaxt, xlab = xlab, xlim = alim,
               yaxs = "i", ylab = flab, ylim = flim,
               col = col, cex = cex,
               type = type, las = 1,
               ...)
        }

      if(identify)
        {
          cat("choose points on the spectrum\n")
          if(.Platform$OS.type == "windows") flush.console()
          id<-identify(x=y,y=x,labels=round(x,2),tolerance=0.15,col="red")
          id.freq<-x[id]
          id.amp<-y[id]
          coord<-list(freq = id.freq ,amp = id.amp)
          return(coord)
        }
    }

                                        # INVISIBLE RETURN DATA
  if(plot == 1 | plot == 2)
    {
      spec<-cbind(x,y)	
      invisible(spec)
    }
                                        # DATA RETURN WHEN NO PLOT
  else if(plot == FALSE) 
    {
      spec<-cbind(x,y)	
      return(spec)
    }
}


################################################################################
##                                SPECPROP
################################################################################

specprop <- function (spec,
                      f = NULL,
                      str = FALSE,
                      flim = NULL,
                      mel = FALSE,
                      plot = FALSE,
                      type = "l",
                      xlab = NULL,
                      ylab = NULL,
                      col.mode = 2,
                      col.quartiles = 4,
                      ...
                      )
    
{
    ## DATA 
    ## Input
    fhz <- f
    if (!is.null(f) & mel) {f <- 2*mel(f/2)}
    if (is.null(f))
        {
            if (is.vector(spec)) stop("'f' is missing")
            else if (is.matrix(spec)) f <- spec[nrow(spec),1]*2000*nrow(spec)/(nrow(spec)-1)
        }

    if (is.matrix(spec))
        {
            freq <- spec[, 1]
            freq = freq * 1000
            spec <- spec[, 2]
        }
    L <- length(spec)
    wl <- L * 2
    if (any(spec < 0)) stop("The frequency spectrum to be analysed should not be in dB")
    if (!is.null(flim)){
        if(flim[1] < 0 || flim[2] > fhz/2) stop("'flim' should range between 0 and f/2")
        if(mel) flim <- mel(flim*1000)/1000
    }
    else {flim <- c(0, f/2000)}

    spec <- spec[(flim[1] * 1000 * wl/f):(flim[2] * 1000 * wl/f)]
    L <- length(spec)
    
    ## Amplitude
    amp <- spec/sum(spec)
    cumamp <- cumsum(amp)
    
    ## Frequency    
    freq <- seq(from = flim[1] * 1000, to = flim[2] * 1000, length.out = L)

    ## RESULTS
    mean <- sum(amp*freq)
    sd <- sqrt(sum(amp*((freq-mean)^2)))
    sem <- sd/sqrt(L)
    median <- freq[length(cumamp[cumamp <= 0.5])+1]
    ## mad deleted at version 1.6.1
    mode <- freq[which.max(amp)]
    Q25 <- freq[length(cumamp[cumamp <= 0.25])+1]
    Q75 <- freq[length(cumamp[cumamp <= 0.75])+1]
    IQR <- Q75 - Q25
    cent <- sum(freq * amp)
    z <- sum(amp - mean(amp))
    w <- sd(amp)
    skew <- (sum((amp - mean(amp))^3)/(L - 1))/w^3
    kurt <- (sum((amp - mean(amp))^4)/(L - 1))/w^4
    sfm <- sfm(amp)
    sh <- sh(amp)
    prec <- f/wl
    
    ## VALUE    
    results <- list(mean = mean, sd = sd, median = median, sem = sem, 
                    mode = mode, Q25 = Q25, Q75 = Q75, IQR = IQR, 
                    cent = cent, skewness = skew, kurtosis = kurt, sfm = sfm, 
                    sh = sh, prec = prec)
    if (str) {
        results <- str(results, digits.d = 5, give.head = FALSE)
    }

    ## PLOT
    if (plot == 1) {if(is.null(xlab)) {if (mel) xlab <- "Frequency (kmel)" else xlab <- "Frequency (kHz)"}
                    if(is.null(ylab)) {ylab <- "Probability"}
                    par(mar = c(5, 5, 4, 2) + 0.1)
                    plot(x = freq/1000, y = amp, type = type, xlab = xlab, 
                         xaxs = "i", ylab = "", yaxs = "i", las = 1, ...)
                    mtext(ylab, side = 2, line = 4)
                    segments(x0 = mode/1000, y0 = 0, x1 = mode/1000, y1 = amp[which(freq == mode)], col = col.mode)
                    segments(x0 = median/1000, y0 = 0, x1 = median/1000, y1 = amp[which(freq == median)], col = col.quartiles)
                    segments(x0 = Q25/1000, y0 = 0, x1 = Q25/1000, y1 = amp[which(freq == Q25)], col = col.quartiles, lty = 2)
                    segments(x0 = Q75/1000, y0 = 0, x1 = Q75/1000, y1 = amp[which(freq == Q75)], col = col.quartiles, lty = 3)
                    legend("topright", legend = c("Q25", "median", "Q75", "mode"), col = c(col.quartiles, col.quartiles, col.quartiles, col.mode), lty = c(2, 1, 3, 1), bty = "n")
                   }
    if (plot == 2) {if(is.null(xlab)) {if (mel) xlab <- "Frequency (kmel)" else xlab <- "Frequency (kHz)"}
                    if(is.null(ylab)) {ylab <- "Cumulated probability"}
                    plot(x = freq/1000, y = cumamp, type = type, xlab = xlab, xaxs = "i", ylab = ylab, yaxs = "i", las = 1, ...)
                    segments(x0 = mode/1000, y0 = 0, x1 = mode/1000, y1 = cumamp[which(freq == mode)], col = col.mode)
                    segments(x0 = 0, y0 = cumamp[which(freq == mode)], x1 = mode/1000, y1 = cumamp[which(freq == mode)], col = col.mode)
                    segments(x0 = median/1000, y0 = 0, x1 = median/1000, y1 = max(cumamp)/2, col = col.quartiles)
                    segments(x0 = 0, y0 = max(cumamp)/2, x1 = median/1000, y1 = max(cumamp)/2, col = col.quartiles)
                    segments(x0 = Q25/1000, y0 = 0, x1 = Q25/1000, y1 = max(cumamp)/4, col = col.quartiles, lty = 2)
                    segments(x0 = 0, y0 = max(cumamp)/4, x1 = Q25/1000, y1 = max(cumamp)/4, col = col.quartiles, lty = 2)
                    segments(x0 = Q75/1000, y0 = 0, x1 = Q75/1000, y1 = max(cumamp) * 3/4, col = col.quartiles, lty = 3)
                    segments(x0 = 0, y0 = max(cumamp) * 3/4, x1 = Q75/1000, y1 = max(cumamp) * 3/4, col = col.quartiles, lty = 3)
                    legend("bottomright", legend = c("Q25", "median", "Q75", "mode"), col =  c(col.quartiles, col.quartiles, col.quartiles, col.mode), lty = c(2, 1, 3, 1), bty = "n")
                }
    if (plot == 1 | plot == 2) {
        invisible(results)
    }
    else if (plot == FALSE) {
        return(results)
    }
}


################################################################################
##                                SPECTRO
################################################################################

spectro <- function(
                  wave,
                  f,
                  wl = 512,
                  wn = "hanning",
                  zp = 0,
                  ovlp = 0,
                  complex = FALSE,
                  norm = TRUE,
                  fftw= FALSE,
                  dB = "max0",
                  dBref = NULL,
                  plot = TRUE,
                  grid = TRUE,
                  osc = FALSE,
                  scale = TRUE,
                  cont = FALSE,
                  collevels = NULL,
                  palette = spectro.colors,
                  contlevels = NULL,
                  colcont = "black",
                  colbg = "white", 
                  colgrid = "black",
                  colaxis = "black",
                  collab = "black",
                  cexlab = 1,
                  cexaxis = 1,   
                  tlab = "Time (s)",
                  flab = "Frequency (kHz)",
                  alab = "Amplitude",
                  scalelab = "Amplitude\n(dB)",
                  main = NULL, 
                  scalefontlab = 1,
                  scalecexlab =0.75,
                  axisX = TRUE,
                  axisY = TRUE,
                  tlim = NULL,
                  trel = TRUE,
                  flim = NULL,
                  flimd = NULL,
                  widths = c(6,1),
                  heights = c(3,1),
                  oma = rep(0,4),
                  listen = FALSE, ...) 
{
    if (!is.null(dB) && all(dB != c("max0", "A", "B", "C", "D"))) 
        stop("'dB' has to be one of the following character strings: 'max0', 'A', 'B', 'C' or 'D'")
    if (complex & norm) {
        norm <- FALSE
        warning("\n'norm' was turned to 'FALSE'")
    }
    if (complex & !is.null(dB)) {
        dB <- NULL
        warning("\n'dB' was turned to 'NULL'")
    }

    input <- inputw(wave = wave, f = f)
    if(!is.null(tlim) && trel && osc) {wave <- wave0 <- input$w} else {wave <- input$w} # necessary to use 'from' and 'to' of soscillo()
    f <- input$f
    rm(input)
    
  if(!is.null(tlim)) wave<-cutw(wave,f=f,from=tlim[1],to=tlim[2])                                      
                                        # dynamic vertical zoom (modifications of analysis parameters)
  if(!is.null(flimd))                        
    {
                                        # zoom magnification
      mag<-round((f/2000)/(flimd[2]-flimd[1]))
                                        # new parameters
      wl<-wl*mag                              
      if(ovlp==0) ovlp<-100
      ovlp<-100-round(ovlp/mag)
                                        # use of normal flim to follow axis modifications
      flim<-flimd     
    }
  

                                        # STFT
  n<-nrow(wave)
  step<-seq(1,n-wl,wl-(ovlp*wl/100))
  z<-stft(wave=wave,f=f,wl=wl,zp=zp,step=step,wn=wn,fftw=fftw,scale=norm, complex=complex)

                                        # X axis settings
  if(!is.null(tlim) && trel) {X<-seq(tlim[1],tlim[2],length.out=length(step))}
  else {X<-seq(0,n/f,length.out=length(step))}

                                        # Y axis settings

  if(is.null(flim)) {Y<-seq(0, (f/2)-(f/wl), length.out=nrow(z))/1000}
  else
    {
      fl1<-flim[1]*nrow(z)*2000/f
      fl2<-flim[2]*nrow(z)*2000/f
      z<-z[(fl1:fl2)+1,]
      Y<-seq(flim[1],flim[2],length.out=nrow(z))
    }

                                        # DB
  if(!is.null(dB))
    {
      if(is.null(dBref)) {z <- 20*log10(z)} else {z <- 20*log10(z/dBref)}
      if(dB!="max0")
          {
            if(dB == "A") z <- dBweight(Y*1000, dBref = z)$A 
            if(dB == "B") z <- dBweight(Y*1000, dBref = z)$B 
            if(dB == "C") z <- dBweight(Y*1000, dBref = z)$C 
            if(dB == "D") z <- dBweight(Y*1000, dBref = z)$D
          }
    }
 
  Z<-t(z)
  
  if(plot)
    {
      if(!isTRUE(norm) && isTRUE(scale)) stop("dB colour scale cannot be plot when 'norm' is FALSE")
      maxz <- round(max(z, na.rm=TRUE))
      if(!is.null(dB))
        {
         if(is.null(collevels)) collevels <- seq(maxz-30, maxz, by = 1)
         if(is.null(contlevels)) contlevels <- seq(maxz-30, maxz, by = 10)
        }
      else
        {
         if(is.null(collevels)) collevels <- seq(0, maxz, length = 30)
         if(is.null(contlevels)) contlevels <- seq(0, maxz, length = 3)
        }        
      
      Zlim<-range(Z, finite = TRUE, na.rm=TRUE)
      # SPECTRO + OSC + SCALE
      if(osc & scale)
        {
          layout(matrix(c(3, 1 ,2, 0), ncol = 2, byrow=TRUE), widths = widths, heights = heights)
          par(las=0, oma=oma, col="white", col=colaxis, col.lab=collab, cex.lab=cexlab, cex.axis=cexaxis, bg=colbg)
          # SCALE
          par(mar=c(0,1,4.5,3))
          dBscale(collevels=collevels,palette=palette,fontlab=scalefontlab,
                  cexlab=scalecexlab,collab=collab,textlab=scalelab,colaxis=colaxis)
          # OSCILLO
          par(mar=c(5,4.1,0,0))
          if(!is.null(tlim) && trel) {wave<-wave0; from<-tlim[1]; to<-tlim[2]} else {from<-FALSE; to <- FALSE}
          soscillo(wave=wave,f=f,bty="u", from=from, to=to,
                   collab=collab,colaxis=colaxis,colline=colaxis,
                   ylim=c(-max(abs(wave)),max(abs(wave))),
                   tickup=max(abs(wave),na.rm=TRUE),
                   tlab=tlab, alab=alab,
                   cexlab=cexlab,
                   cexaxis=cexaxis,
                   ...)
          # SPECTRO
          par(mar=c(0,4.1,1,0), las=1, cex.lab=cexlab+0.2)
          filled.contour.modif2(x=X ,y=Y, z=Z, levels=collevels, nlevels=20,
                                plot.title=title(main=main,xlab="",ylab=flab),
                                color.palette=palette,
                                axisX=FALSE, axisY=axisY
                                )
          if(grid) grid(nx=NA, ny=NULL, col=colgrid)
          if(cont){contour(X,Y,Z,add=TRUE,levels=contlevels,nlevels=5,col=colcont,...)}
          if(colaxis != colgrid) abline(h=0,col=colaxis) else abline(h=0,col=colgrid)
        }
           
      # SPECTRO + SCALE
      if(osc==FALSE & scale)
        {
          layout(matrix(c(2, 1), ncol = 2, byrow=TRUE), widths = widths)
          # SCALE
          par(mar=c(5,1,4.5,3), oma=oma, las=0, bg=colbg)
          dBscale(collevels=collevels,palette=palette,fontlab=scalefontlab,
                  cexlab=scalecexlab,collab=collab,textlab=scalelab,colaxis=colaxis)
          # SPECTRO
          par(mar=c(5,4.1,1,0),las=1,cex=1,col=colaxis,col.axis=colaxis,col.lab=collab,bg=colbg,cex.lab=cexlab+0.2)
          filled.contour.modif2(x=X ,y=Y, z=Z, levels=collevels, nlevels=20,
                                plot.title=title(main=main,xlab=tlab,ylab=flab),
                                color.palette=palette,
                                axisX=axisX, axisY=axisY)
          if(grid) grid(nx=NA, ny=NULL, col=colgrid)
          if(colaxis!=colgrid) abline(h=0,col=colaxis) else abline(h=0,col=colgrid)
          if(cont){contour(X,Y,Z,add=TRUE,levels=contlevels,nlevels=5,col=colcont,...)}
        }

      # SPECTRO + OSCILLO
      if(osc & scale==FALSE)
        {
          layout(matrix(c(2,1), nrow = 2, byrow=TRUE), heights=heights) 
          par(mar=c(5.1,4.1,0,2.1), las=0, oma=oma, bg=colbg)
          # OSCILLO
          if(!is.null(tlim) && trel) {wave <-wave0; from<-tlim[1]; to<-tlim[2]} else {from<-FALSE ; to<-FALSE}
          soscillo(wave=wave,f=f,bty="u", from=from, to=to,  
                   collab=collab,colaxis=colaxis,colline=colaxis,
                   tickup=max(abs(wave),na.rm=TRUE),
                   ylim=c(-max(abs(wave)),max(abs(wave))),
                   tlab=tlab, alab=alab,
                   cexlab=cexlab, cexaxis=cexaxis,
                   ...)
          # SPECTRO
          par(mar=c(0,4.1,2.1,2.1), las=1, cex.lab=cexlab)
          filled.contour.modif2(x=X ,y=Y, z=Z, levels=collevels, nlevels=20,
                                plot.title=title(main=main,xlab="",ylab=flab), color.palette=palette, axisX=FALSE, axisY=axisY,
                                col.lab=collab,colaxis=colaxis,...)		
          if(grid) grid(nx=NA, ny=NULL, col=colgrid)
          if(cont){contour(X,Y,Z,add=TRUE,levels=contlevels,nlevels=5,col=colcont,...)}
          if(colaxis!=colgrid) abline(h=0,col=colaxis) else abline(h=0,col=colgrid)
        }

      # SPECTRO ONLY
      if(osc==FALSE & scale==FALSE)
        {
          par(las=1, col=colaxis, col.axis=colaxis, col.lab=collab, bg=colbg, cex.axis=cexaxis, cex.lab=cexlab,...)
          filled.contour.modif2(x=X ,y=Y, z=Z, levels=collevels, nlevels=20,
                                plot.title=title(main=main,xlab=tlab,ylab=flab), color.palette=palette, axisX=axisX, axisY=axisY,
                                col.lab=collab,colaxis=colaxis)		
          if(grid) grid(nx=NA, ny=NULL, col=colgrid)
          if(cont){contour(X,Y,Z,add=TRUE,levels=contlevels,nlevels=5,col=colcont,...)}
          if(colaxis!=colgrid) abline(h=0,col=colaxis) else abline(h=0,col=colgrid)
        } 
      if(listen) {listen(wave, f=f)}
      invisible(list(time=X, freq=Y, amp=z))
    }
  else return(list(time=X, freq=Y, amp=z))
}



################################################################################
##                                SPECTRO3D
################################################################################

spectro3D<-function(
                    wave,
                    f,
                    wl = 512,
                    wn = "hanning",
                    zp = 0,
                    ovlp = 0,
                    norm = TRUE,
                    fftw = FALSE,
                    dB = "max0",
                    dBref = NULL,
                    plot = TRUE,
                    magt = 10,
                    magf = 10,
                    maga = 2,
                    palette = reverse.terrain.colors
                    )

{
                                        # STOP MESSAGES
  if(!is.null(dB) && all(dB!=c("max0","A","B","C","D")))
    stop("'dB' has to be one of the following character strings: 'max0', 'A', 'B', 'C' or 'D'")
  if(!is.null(wl) & wl%%2 == 1) stop("'wl' has to be an even number.")
  
                                        # INPUT
  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; rm(input)
  
                                        # STFT
  n <- nrow(wave)
  step <- seq(1, n - wl, wl - (ovlp * wl/100))
  z <- stft(wave = wave, f = f, wl = wl, zp = zp, step = step, wn = wn, fftw=fftw, scale=norm)

                                       # FREQUENCY DATA
  F <- seq(0, (f/2)-(f/wl), length.out=nrow(z)) / 1000
  
                                        # DB
  if(!is.null(dB))
    {
      if(is.null(dBref)) {z <- 20*log10(z)} else {z <- 20*log10(z/dBref)}
      if(dB!="max0")
          {
            if(dB == "A") z <- dBweight(Y*1000, dBref = z)$A 
            if(dB == "B") z <- dBweight(Y*1000, dBref = z)$B 
            if(dB == "C") z <- dBweight(Y*1000, dBref = z)$C 
            if(dB == "D") z <- dBweight(Y*1000, dBref = z)$D
          }
    }
 
  if(plot)
    {        
      X <- magt * (1:ncol(z))
      Y <- magf * (1:nrow(z))
      Z <- maga * z
      Xat <- seq(magt, magt * ncol(z), by = (magt * ncol(z))/4)
      Yat <- seq(magf, magf * nrow(z), by = (magf * nrow(z))/4)
      if(is.null(dB))
        {
         Zat <- seq(min(Z, na.rm=TRUE), max(Z, na.rm=TRUE), length.out = 5)
         Zlab <- as.character(round(seq(min(Z, na.rm=TRUE)/maga, max(Z, na.rm=TRUE)/maga, length.out=5), 1))
        }
      else
        {
         Zat <- seq(min(Z, na.rm=TRUE), maga, by = abs(min(Z, na.rm=TRUE))/4)
         Zlab <- as.character(round(seq(min(Z, na.rm=TRUE)/maga, 0, by = abs(min(Z, na.rm=TRUE))/(4*maga)), 1))
        }
      Xlab <- as.character(round(seq(0, n/f, by = n/(4 * f)),1))
      Ylab <- as.character(round(seq((f/1000)/(wl + zp), f/2000, by = f/(4*2000)), 1))
      Zlim <- range(Z, na.rm=TRUE)
      Zlen <- Zlim[2] - Zlim[1] + 1
      colorlut <- palette(Zlen)
      col <- colorlut[Z - Zlim[1] + 1]
      rgl::rgl.clear()
      rgl::rgl.bbox(color = "white", emission = "gray8", specular = "gray",
               shininess = 50, alpha = 0.8, xat = Yat, xlab = Ylab,
               xunit = 0, yat = Zat, ylab = Zlab, yunit = 0, zat = Xat,
               zlab = Xlab, zunit = 0)
      rgl::rgl.texts(x = 1, z = magt * ncol(Z)/2, y = min(Z), text = "Time (s)", color = "white")
      rgl::rgl.texts(z = 1, x = magf * nrow(Z)/2, y = min(Z), text = "Frequency (kHz)", color = "white")
      rgl::rgl.texts(x = 1, z = 0, y = 0, text = "Amplitude (dB)", color = "white")
      rgl::rgl.surface(Y, X, Z, color = col, back = "lines")
      invisible(list(time=seq(0,n/f,length.out=length(step)), freq=F, amp=z))
    }
  else return(list(time=seq(0,n/f,length.out=length(step)), freq=F, amp=z))
}


################################################################################
##                                SQUAREFILTER
################################################################################

squarefilter <- function (f,
                          from = NULL,
                          to = NULL,
                          bandpass = TRUE,
                          wl = 1024 
                          )

{
    if (length(from) != length(to)) stop ("The arguments 'from' and 'to' should have the same length") 
    if (bandpass) {a <- rep(0, wl/2)}
    else {a <- rep(1, wl/2)}
    if (length(from) == 0) {return(a)}
    else for (i in 1:length(from))
             {if (from[i] < 0 | to[i] < 0 | from[i] > f | to[i] > f)
                  stop("For all i, from[i] and to[i] should range between 0 and f")
              if (from[i] >= to[i])
                  stop("For all i, from[i] should be inferior to to[i]")
              if (i < length(from))
                  {if (to[i] > from[i+1])
                       warning("For all i, to[i] should be inferior to from[i+1]")}
              F <- round(wl * (from[i]/f))
              T <- round(wl * (to[i]/f))
              a[T:F] <- abs(a[T:F]-1)}
    return(freq=cbind(seq(0,f/1000,length.out = wl/2+1)[1:(wl/2)], amp=a))
}


################################################################################
##                                STFT.EXT
################################################################################

stft.ext <- function(
                  file,
                  wl = 512,
                  ovlp = 0,
                  mean = FALSE,
                  norm = FALSE,
                  dB = FALSE,
                  verbose = FALSE)
{
  .Call("stft", file, ovlp, wl, mean, norm, dB, verbose, PACKAGE = "seewave") 
}


################################################################################
##                                SYMBA
################################################################################

symba<-function(
                x,
                y = NULL,
                symb = 5,
                collapse = TRUE,
                entropy = "abs",
                plot = FALSE,
                type = "l",
                lty1 = 1,
                lty2 = 2,
                col1 = 2,
                col2 = 4,
                cex1 = 0.75,
                cex2 = 0.75,
                xlab = "index",
                ylab = "Amplitude",
                legend = "TRUE",
                ...)

{
                                        # input x
  s1<-discrets(x=x, symb=symb, collapse=FALSE)
                                        # frequency of each symbols
  freq1a<-table(s1,dnn="symbol frequency in the sequence")
                                        # entropy of the sequence
  freq1b<-as.vector(freq1a)
  freq1c<-freq1b/sum(freq1b)
  if(entropy=="abs") h1<- -sum(freq1c*log(freq1c))
  else if(entropy=="rel") h1<- -sum(freq1c*log(freq1c))/log(length(freq1c))
  
  if(is.null(y))
    {
      if(plot)
        {
          if(symb==3) {s1<-c(NA,s1)} else if(symb==5) {s1<-c(NA,s1,NA)}
          plot(x,type=type,lty=lty1,xlab=xlab,ylab=ylab,...)
          text(x=x,labels=s1,col=col1, cex=cex1)
        }
      if(collapse) s1<-paste(s1,collapse="")
      results <- list(s1=s1,freq1=freq1c, h1=h1)
      if(plot) {invisible(results)} else return(results)
    }

  if(!is.null(y))
    {
      if(length(x)!=length(y)) {stop("x and y should have the same length")}
                                        # input y
      y<-inputw(wave = y,f = NULL)$w
      s2<-discrets(y, symb=symb, collapse=FALSE)
                                        # frequency of each symbols
      freq2a<-table(s2,dnn="symbol frequency in the sequence")
                                        # entropy of the sequence
      freq2b<-as.vector(freq2a)
      freq2c<-freq1b/sum(freq2b)
      if(entropy=="abs") h2<- -sum(freq2c*log(freq2c))
      else if(entropy=="rel") h2<- -sum(freq2c*log(freq2c))/log(length(freq2c))
                                        # joint entropy
                                        # frequency of each pair of symbols
      freq12<-table(paste(s1,s2,sep=""))
      freq12<-as.vector(freq12)
      freq12<-freq12/sum(freq12)
                                        # joint entropy
      if(entropy=="abs") h12<- -sum(freq12*log(freq12))
      else if(entropy=="rel") h12<- -sum(freq12*log(freq12))/log(length(freq12))
                                        # mutual information
      I<-h1+h2-h12
      if(plot)
        {
          if(symb==3) {s1<-c(NA,s1);s2<-c(NA,s2)} else {s1<-c(NA,s1,NA);s2<-c(s2,NA)}
          plot(x,type=type, col=col1, lty=lty1,ylim=c(min(c(x,y)),max(c(x,y))),xlab=xlab,ylab=ylab,...)
          text(x=x,labels=s1, col=col1, cex=cex1)
          lines(y,type=type, col=col2,lty=lty2)
          text(x=y,labels=s2, col=col2, cex=cex2)
        }
      if(collapse) {s1<-paste(s1,collapse=""); s2<-paste(s2,collapse="")} 
      results <- list(s1=s1,freq1=freq1c, h1=h1, s1=s2,freq2=freq2c, h2=h2, I=I)
      if(plot) {invisible(results)} else return(results)      
    }
}


################################################################################
##                                SYNTH
################################################################################

synth <- function (f,
                   d,
                   cf,
                   a = 1,
                   signal = "sine",
                   shape = NULL,
                   p = 0,
                   am = c(0, 0, 0), 
                   fm = c(0, 0, 0, 0),
                   harmonics = 1,
                   plot = FALSE,
                   listen = FALSE,
                   output = "matrix",
                   ...
                   )
    
{
    if(length(am) == 2) {am <- c(am, 0)}
    if(length(fm) == 3) {fm <- c(fm, 0)}
    if (harmonics[1] != 1) {stop("Harmonics must begin with a value equal to 1")}
    else if (!all(harmonics <= 1)) {warning ("One harmonic has a higher amplitude than the fundamental")}
    l1 <- length(harmonics)
    l2 <- round(f*d)
    t <- array(0, c(l1,l2))
    for (i in 1:l1)
        {t[i,] <- synth0(f, d, i*cf, a*harmonics[i], shape, i*p, am, fm, signal=signal)}
    sound <- apply(t, 2, sum)
    sound <- outputw(wave = sound, f = f, format = output)
    if (plot) {
        spectro(sound, f = f, ...)
        if (listen) {
            listen(sound, f = f)
        }
        invisible(sound)
    }
    else {
        if (listen) {
            listen(sound, f = f)
        }
        return(sound)
    }
}

aux <- function(x, signal){
 if (signal == 'sine')   {return (sin(x))}
 if (signal == 'square') {return ((-1)^floor(x/pi))}
 if (signal == 'tria')   {return ((2*x/pi-2*round(x/pi))*(-1)^round(x/pi))}
 if (signal == 'saw')    {return (2*(x/(2*pi)-floor((x/(2*pi))+0.5)))}
}

synth0 <- function (f, d, cf, a = 1, shape = NULL, p = 0, am = c(0, 0, 0), 
                    fm = c(0, 0, 0, 0), signal, ...) 

{
    n <- round(f * d)
    amp <- am[1]/100
    amf <- am[2]
    amP <- am[3]
    fme <- fm[1]
    fmf <- fm[2]
    fmE <- fm[3]
    fmP <- fm[4]
    t <- seq(0, d * 2 * pi, length.out = n)
    if (fme == 0 && fmf != 0) stop("FM sinusoidal excursion has to be set")
    if (fme != 0 && fmf == 0 && fmE == 0) stop("FM sinusoidal frequency or FM linear excursion has to be set")
    if (fme != 0 && fmf == 0 && fmE != 0) stop("FM sinusoidal frequency has to be set")
    if (fmE > 0) {freq <- seq(0, fmE/2, length.out = n)} else {freq <- rev(seq(fmE/2, 0, length.out = n))}
    if (fme == 0 & fmf == 0) {sound <- (1 + amp * cos(amf * t + amP)) * aux((cf * t) + (freq * t) + p, signal = signal)}
    if (fme != 0 & fmf != 0) {
        if (fmE == 0) {sound <- (1 + amp * cos(amf * t + amP)) * aux(cf * t + (fme/fmf) * sin(fmf * t + fmP) + p, signal = signal)}
        else sound <- {(1 + amp * cos(amf * t)) * aux(cf * t + (fme/fmf) * sin(fmf * t + fmP) + (freq * t) + p, signal = signal)}
    }
    if (!is.null(shape)) {
        if (shape == "incr") {S <- seq(0, 1, length.out = n)}
        if (shape == "decr") {S <- seq(1, 0, length.out = n)}
        if (shape == "sine") {S <- sin(seq(0, pi, length.out = n))}
        if (shape == "tria") {
            if (n%%2 == 1) S <- c(seq(0, 1, length.out = n%/%2), seq(1, 0, length.out = n%/%2 + 1))
            else S <- c(seq(0, 1, length.out = n%/%2), seq(1,0, length.out = n%/%2))
        }
        sound <- S * sound
    }
    sound <- a * (sound/max(abs(sound)))
    return(sound)
}



################################################################################
##                                SYNTH2
################################################################################

synth2 <- function(
    env = NULL,
    ifreq,    # in Hz
    f,   # in Hz
    plot = FALSE,
    listen = FALSE,
    output = "matrix",
    ...
    )

{
    ## input
    nfreq <- length(ifreq)
    if(is.null(env)) {env <- rep(1, times=nfreq)}    
    nenv <- length(env)
    if(nenv != nfreq) {stop("'wave' and 'freq' should have exactly the same length")}

    ## get the instantaneous phase from the instantaneous frequency
    ## equivalent to the integral (primitive) of the instantaneous phase along time
    iphase <- 2*pi*cumsum(ifreq)/f

    ## synthesis
    sound <- env*cos(iphase)   
    sound <- outputw(wave = sound, f = f, format = output)

    ##output
    if (plot) {
        spectro(sound, f = f, ...)
        if (listen) {
            listen(sound, f = f)
        }
        invisible(sound)
    }
    else {
        if (listen) {
            listen(sound, f = f)
        }
        return(sound)
    }
}

################################################################################
##                                TH
################################################################################

th <- function(
             env
             )

{
  options(warn=-1)  # to ignore warning messages
  if(is.na(env)) z <- 0 
  else{
    options(warn=0) # to authorize warning messages
    if(any(env<0, na.rm=TRUE))  stop ("data must be an envelope, i. e. a vector including positive values only.")
    N<-length(env)
    env<-env/sum(env) # PMF
    if(any(is.nan(env))) 
      {
        warning("Caution! There is no signal in this data set! The temporal entropy is null!", call.=FALSE)
        return(0)
      }
    if(sum(env)/(N*env[1]) == 1 | sum(env)/N == 1)
      {
        warning("Caution! This is a square signal. The temporal entropy is null!", call.=FALSE)
        return(0)
      }
    env[env==0]<-1e-7
    z<--sum(env*log(env))/log(N)
  }
  return(z)
}


################################################################################
##                                TIMER
################################################################################

timer <- function(wave,
                   f,
                   threshold = 5,
                   dmin = NULL, 
                   envt = "abs",
                   power = 1,
                   msmooth = NULL,
                   ksmooth = NULL,
                   ssmooth = NULL,
                   tlim = NULL,
                   plot = TRUE,
                   plotthreshold = TRUE,
                   col = "black",
                   colval = "red",
                   xlab = "Time (s)",
                   ylab = "Amplitude",
                   ...
                   )
{
    ## INPUT
    input <- inputw(wave = wave, f = f)
    wave <- input$w
    f <- input$f
    rm(input)
    n <- length(wave)
    thres <- threshold/100
  
    ## STOP MESSAGES
    if(power == 0) stop("'power' cannot equal to 0")
    if(!is.null(msmooth) && !is.null(ksmooth)) stop("'msmooth' and 'ksmooth' cannot be used together")
    if(!is.null(msmooth) && !is.null(ssmooth)) stop("'msmooth' and 'ssmooth' cannot be used together")
    if(!is.null(ksmooth) && !is.null(ssmooth)) stop("'ksmooth' and 'ssmooth' cannot be used together")
    if(!is.null(dmin)){
      if(length(dmin)!=1) stop("'dmin' should be a numeric vector of length 1") 
      if(dmin <= 0) stop("'dmin' cannot be negative or equal to 0")
      if(dmin >= n/f) stop("'dmin' cannot equal or be higher than the wave duration")
    } 
    ## TIME LIMITS
    if(!is.null(tlim)) {
      wave <- cutw(wave, f=f, from=tlim[1], to=tlim[2])
      n <- length(wave)
    }

    ## ENVELOPE
    wave1 <- env(wave = wave, f = f, msmooth = msmooth, ksmooth =
             ksmooth, ssmooth=ssmooth, envt = envt, norm = TRUE, plot = FALSE)
    n1 <- length(wave1)
    f1 <- f*(n1/n)

    ## SIGNAL DETECTION
    if(power!=1) wave1 <- wave1^power
    ## signal <= threshold = 1 and signal > threshold = 2
    wave2 <- ifelse(wave1 <= thres, yes = 1, no = 2)
    ## 2 = silence, 3 = change,  4 = signal
    n2 <- length(wave2)
    wave4 <- apply(as.matrix(1:(n2-1)), 1, function(x) wave2[x] + wave2[x+1])
    ## change at the start and end of the signal
    n4 <- length(wave4)
    wave4[c(1,n4)] <- 3
    ## look for change positions
    wave5 <- which(wave4 == 3)
    ## time threshold (dmin)
    if(!is.null(dmin)){
        ## look for events < dmin that have to be excluded
        event.dur <- diff(wave5)
        event.idx <- which(event.dur < dmin*f1)
        ## condition to be sure that there are numeric values in event.idx
        if(length(event.idx)!=0){
        ## replaces all events to be exluded by 2 (= pause)
        for(i in event.idx) {wave4[(wave5[i]):(wave5[i]+event.dur[i])] <- 2}
##        for(i in event.idx) {wave4[(wave5[i]-1):(wave5[i]+event.dur[i]+1)] <- 2}
        wave4[which(abs(diff(wave4))==2)] <- 3
        wave4[c(1,n4)] <- 3
        }
        wave5 <- which(wave4 == 3)
        if(length(wave5)==2){stop("'dmin' was set to a too high value, there are no signal longer than 'dmin'")}
    }
    
    wave5[-1] <- wave5[-1] + 1 # +1 for all positions except the first one (position 1,  0 s).
    
    ## DURATIONS OF THE DIFFERENT TIME ITEMS
    f4 <- f * (n4/n)
    wave4 <- ts(wave4, start = 0, end = n4/f4, frequency = f4)
    positions <- time(wave4)[wave5]
    durations <- diff(positions)
    npos <- length(positions) 
    ## the wave starts with a pause
    if(wave2[1] == 1) {
      first <- "pause"
      pause <- durations[seq(1, npos - 1, by = 2)]
      signal <- durations[seq(2, npos - 1, by = 2)]
      start.signal <- positions[seq(2, npos - 1, by = 2)]
      end.signal <- positions[seq(3, npos - 1, by = 2)]
    }
    ## the wave starts with a signal
    else {
      first <- "signal"
      pause <- durations[seq(2, npos - 1, by = 2)]
      signal <- durations[seq(1, npos - 1, by = 2)]
      start.signal <- positions[seq(1, npos - 1, by = 2)]
      end.signal <- positions[seq(2, npos - 1, by = 2)]
    }

    ratio <- sum(signal)/sum(pause)
    timer <- list(s = signal, p = pause, r = ratio,
                  s.start = start.signal, s.end = end.signal,
                  first = first)

    ## PLOT
    if(plot) {
        plot(x=seq(0, n1/f1, length.out=n1), y=wave1, xlab = xlab, ylab = ylab, yaxt = "n",
             ylim = c(0, 1 + 0.1), col = col, type = "l", xaxs ="i", ...)
        if(plotthreshold) {
            abline(h = thres, col = colval, lty = 2)
            mtext(paste(as.character(threshold), "%"), side = 2,
                line = 0.5, at = thres, las = 1, col = colval,
                cex = 0.8)
        }

        outline <- wave4-3
        outline[outline < 0] <- 0
        outline[1] <- NA
        lines(x=seq(0, n1/f1, length.out=n1), y=c(outline), col=colval)
        wave8 <- numeric(npos - 1)
        for (i in 2:npos) {
            wave8[i] <- ((wave5[i] - wave5[i - 1])/2) + wave5[i -
                1]
        }
        if(wave2[1] == 1) {
            wave8.1 <- wave8[seq(2, npos, by = 2)]/f1
            wave8.2 <- wave8[seq(3, npos, by = 2)]/f1
        }  else {
            wave8.2 <- wave8[seq(2, npos, by = 2)]/f1
            wave8.1 <- wave8[seq(3, npos, by = 2)]/f1
        }
        ypl <- as.character(round(pause, 2))
        ysl <- as.character(round(signal, 2))
        text(x = wave8.1, y = 0.075, ypl, col = colval, cex = 0.8)
        text(x = wave8.2, y = 1.075, ysl, col = colval, cex = 0.8)
        invisible(timer)
    }
    else {
        return(timer)
    }
}


################################################################################
##                                TKEO
################################################################################

TKEO <- function(
                 wave,
                 f,
                 m=1,
                 M=1,
                 plot = TRUE,
                 xlab = "Time (s)",
                 ylab = "Energy",
                 type = "l",
                 bty = "l",
                 ...
                 )
    
{
    input <- inputw(wave = wave, f = f)
    wave <- input$w
    f <- input$f
    n <- length(wave)
    rm(input)
    t <- seq(0, n/f, length.out=n)
    y <- wave^(2/m) - (c(wave[-(1:M)], rep(NA,M)) * c(rep(NA,M), wave[1:(length(wave)-M)]))^(1/m)
    results <- cbind(time=t, tkeo=y)
    if(plot) {
    plot(x=t, y=y, xlab=xlab, ylab=ylab, type=type, bty=bty, xaxs="i", ...)
    invisible(results)
    }
    else return(results)
}


################################################################################
##                                WASP
################################################################################

wasp<-function(
               f,
               t = 20,
               c = NULL,
               s = NULL,
               d = NULL,
               medium = "air"
               )

{
  if(medium == "air")
    {
      if(!is.null(d)) stop("Depth (d) is not a valuable argument for air medium")
      if(!is.null(s)) stop("Salinity (s) is not a valuable argument for air medium")
      if(!is.null(c)) C<-c
      else C<-331.4+0.6*t
    }
  
  if(medium == "sea")
    {
      if(!is.null(c)) C<-c
      if(is.null(s))  stop("Please specify a salinity value (parts per thousand) for sea medium")
      if(is.null(d))  stop("Please specify a depth value (m) for sea medium")
      else 
        {
          C<-1448.96+4.591*t-(5.304e-2)*t^2+(2.374e-4)*t^3+1.34*(s-35)+(1.63e-2)*d+(1.675e-7)*d^2-(1.025e-2)*t*(s-35)-(7.139e-13)*t*d^3
        }  
    }
  
  if(medium == "fresh")
    {
      if(!is.null(c)) C<-c
      if(!is.null(d)) stop("Depth (d) is not a valuable argument for freshwater medium")
      if(!is.null(s)) stop("Salinity (s) is not a valuable argument for freshwater medium")
      else 
        {
          C<-1.402385e3+5.038813*t-(5.799136e-2)*t^2+(3.287156e-4)*t^3-(1.398845e-6)*t^4+(2.787860e-9)*t^5
        }
    }

  lambda<-C/f  
  results<-list(l=lambda,c=C)
  return(results)
}


################################################################################
##                                WAV2FLAC
################################################################################

wav2flac <- function(file,
                     reverse = FALSE,
                     overwrite = FALSE,
                     exename = NULL,
                     path2exe = NULL
                     )
    
{
    ## UNIX
    if (.Platform$OS.type == "unix")
        {
        if (missing(exename)) {exename <- "flac"}
        if (missing(path2exe)) {exe <- exename}
        else {exe <- paste(path2exe, exename, sep = "/")}
        
                                        # Give specific error when flac is not found
        if (system(paste(exe, "-v"), ignore.stderr = TRUE)!=0){stop("FLAC program was not found.")}

                                        # Check for reverse
        if (reverse) {e <- system(paste(exe, "-d", file), ignore.stderr = TRUE)}
        else{e <- system(paste(exe, file), ignore.stderr = TRUE)}
        }
    
    ## WINDOWS
    if (.Platform$OS.type == "windows") {
        if (missing(exename)) {exename <- "flac.exe"}
        if (missing(path2exe)) {
                                        # Drive letter in caps, otherwise it causes an error
            exe <- paste("C:/Program Files/FLAC/", exename, sep = "")
            if (!file.exists(exe)){exe <- paste("C:/Program Files (x86)/FLAC/", exename, sep = "")} #For 64bit systems
        }
        else {exe <- paste(path2exe, exename, sep = "/")}
                                        # Give specific error when flac is not found
        if (!file.exists(exe)){stop("FLAC program was not found.")}
        
        if (reverse) {e <- system(paste(shQuote(exe), "-d", shQuote(file, type = "cmd"), sep = " "), ignore.stderr = TRUE)}
        else {e <- system(paste(shQuote(exe), shQuote(file, type = "cmd"), sep = " "), ignore.stderr = TRUE)}
    }

#    if (e > 0) {stop("File not found or wrong format/encoding")}   ## error should be e < 0
    if (overwrite) {unlink(file)}
}


################################################################################
##                                WF
################################################################################

wf <- function(
               wave,
               f = NULL,
               wl = 512,
               zp = 0,
               ovlp = 0,
               fftw = FALSE,
               dB = "max0",
               dBref = NULL,
               wn = "hanning",
               x = NULL,
               hoff = 1,
               voff = 1,
               col = heat.colors,
               xlab = "Frequency (kHz)",
               ylab = "Amplitude (dB)",
               xaxis = TRUE,
               yaxis = TRUE,
               density = NULL,
               border = NULL,
               lines = FALSE,
               lwd = NULL,
               ...)

{
                                        # STOP MESSAGES
  if(missing(wave) && is.null(x)) stop("'wave' or 'x' has to be set up")
  if(!missing(wave) && !is.null(x)) stop("'wave' or 'x' has to set up, not both of them!")
  if(lines && !is.null(border)) stop("when 'lines' is TRUE, changing 'border' has no effect")
  if(lines && !is.null(density)) stop("when 'lines' is TRUE, changing 'density' has no effect")
  if(hoff < 0) stop("'hoff' cannot be negative")
  if(voff < 0) stop("'voff' cannot be negative")
  if(!is.null(dB) && all(dB!=c("max0","A","B","C","D")))
    stop("'dB' has to be one of the following character strings: 'max0', 'A', 'B', 'C' or 'D'")

                                        # INPUT
  if(is.null(x))
    {
      input <- inputw(wave = wave, f = f)
      wave <- input$w
      f <- input$f
      rm(input)
      n <- nrow(wave)
      step <- seq(1, n - wl, wl - (ovlp * wl/100))
      data <- stft(wave = wave, f = f, wl = wl, zp = zp, step = step, 
                   wn = wn, fftw = fftw)
                                        # dB WEIGHTS
      F <- seq((f/1000)/(wl + zp), f/2000, length.out=nrow(data))
  if(!is.null(dB))
    {
      if(is.null(dBref)) data <- 20*log10(data) else data <-20*log10(data/dBref)
      if(dB == "max0") data <- data
      if(dB == "A") data <- dBweight(F*1000, dBref = data)$A 
      if(dB == "B") data <- dBweight(F*1000, dBref = data)$B 
      if(dB == "C") data <- dBweight(F*1000, dBref = data)$C 
      if(dB == "D") data <- dBweight(F*1000, dBref = data)$D
    }

    }
  else
    {
      data <- x
    }

                                        # colors
  if(!is.function(col)) {col <- rep(col,dim(data)[2])} else {col <- col(dim(data)[2])}

                                        # seek for axes ticks: look for the column including the maximum value
                                        # and get the values for both axes with a 'virtual' plot

  if(xaxis | yaxis)
    {
      maxi <- which(data==max(data),arr.ind=TRUE)[1,2]
      if(!is.null(f)) {x<-seq(0,f/2000,length.out=nrow(data))} else {x<-1:nrow(data)}
      plot(x=x,y=data[,maxi],
           type="n", xaxt="n", yaxt="n",xlab="",ylab="", bty="n")
      atX <- axis(side=1, labels=FALSE, tick=FALSE)
      atY <- axis(side=2, labels=FALSE, tick=FALSE)
    }

                                        # keep code readable
  pht <- dim(data)[2]
  pwid <- dim(data)[1]
  hoff <- floor(hoff)

                                        # create empty plot frame of useful size
  par(las=1)
  plot(c(1, (pwid+hoff*pht)), c((range(data)[1]-2), (pht*voff+range(data)[2])),
       type="n",
       xlab=xlab, xaxs="i", xaxt="n",
       ylab=ylab, yaxs="i", yaxt="n",...)

                                        # apply horiz, vert offsets to input data
  ywf<-(sapply(seq(1,pht),function(a) replace(rep(a*voff,(pwid+pht*hoff)),
                                              seq(((pht-a)*hoff+1),((pht-a)*hoff+pwid)), (data[,a]+a*voff))))
  xwf<-seq(1,pwid+hoff*pht)

                                        # lines
  if(lines)
    {
      ywf <- apply(ywf, MARGIN=2, function(a) replace(a,which(a == max(a)), NA))
      invisible(mapply(function(x,cols) lines(ywf[,x],col=cols,lwd=lwd), seq(pht,1,-1), col))
    }

                                        # polygons
  else
    {
      invisible(mapply(function(x,cols) polygon(c(xwf,rev(xwf)),
                                                c((ywf[,x]), rep((range(data)[1]-4), length(ywf[,x]))),
                                                col=cols, border=border, density=density, lwd=lwd), seq(pht,1,-1), col))
    }

                                        # axes
  if(xaxis)
    {
      from <- dim(ywf)[1]-dim(data)[1]
      len <- length(atX)
      if(!is.null(f)) {to <- 2000*atX[len]*dim(ywf)[1]/f} else {to <- dim(ywf)[1]}
      at <- seq(from = from, to = to, length.out = len)
      labels <- as.character(atX)
      axis(side = 1, at = at, labels = labels)
    }
  if(yaxis){axis(side=2, at=atY+1, labels=as.character(atY))}
  box()
  par(las=0)
}





################################################################################
##                                ZAPSILW
################################################################################

zapsilw<-function(
                  wave,
                  f,
                  threshold = 5,
                  plot = TRUE,                  
                  output = "matrix",
                  ...)

{
  input <- inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; bit <- input$bit ; rm(input)

  wave1 <- afilter(wave,f=f,threshold=threshold,plot=FALSE)
  wave2 <- wave1[wave1!=0]
  wave2 <- outputw(wave=wave2, f=f, format=output)

  if(plot)
    {
      def.par <- par(no.readonly = TRUE)
      on.exit(par(def.par))    
      par(mfrow=c(2,1),oma=c(0,0.1,0,0))
      oscillo(wave=wave, f=f, ...)
      title(main="original")
      oscillo(wave=wave2, f=f,...)
      title(main="silence removed")
      invisible(wave2)
    }
  else
    {
      return(wave2)
    }
}


################################################################################
##                                ZC
################################################################################

zc<-function(
             wave,
             f,
             plot = TRUE,
             interpol = 1,
             threshold = NULL,
             xlab = "Time (s)",
             ylab = "Frequency (kHz)",
             ylim = c(0,f/2000),
             warning = TRUE,
             ...)

{
  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; rm(input)

  if(isTRUE(warning) & interpol > 5)
    {
      cat("please wait...")
      if(.Platform$OS.type == "windows") flush.console()
    }

  n<-nrow(wave)

  if(!is.null(threshold)) wave<-afilter(wave=wave,f=f,threshold=threshold,plot=FALSE)
  
  if(interpol > 1)
    {
      waveinterpol<-approx(wave,n=n*interpol)
      wave<-as.matrix(waveinterpol$y)
      F<-f*interpol
    }
  else F<-f

                                        # replaces null or positive values by 1 and negative values by 2  
  wave1<-ifelse(wave>=0,1,2) 

                                        # adds successive values in wave1, values of 3 corresponds to ZC
  wave2<-numeric(n*interpol)  
  for (i in 1:((n*interpol)-1))  {wave2[i]<-wave1[i]+wave1[i+1]}
                                        # replaces 2 by 0
  wave3<-ifelse(wave2==2, yes=0, no=wave2)
                                        # replaces 4 by 0
  wave4<-ifelse(wave3==4, yes=0, no=wave3)
                                        # replaces 3 by their index
  wave5<-replace(wave4,which(wave4==3),which(wave4==3))

                                        # computes the period T between two successive zc 
  wave6<-which(wave2==3) 
  nn<-length(wave6) 
  wave7<-numeric(nn)
  for (i in 2:(nn-1)) {wave7[i]<-wave6[i+1]-wave6[i-1]}

                                        # replaces index by T
  wave8<-replace(wave5,which(wave5!=0),wave7)

                                        # calculates the frequency
  wave9<-F/(wave8)/1000 
  y<-replace(wave9,which(wave9==Inf),NA)

  x<-seq(0,n/f,length.out=n*interpol)
  if(plot)
    {
      plot(x = x, y = y, xlab=xlab, ylab=ylab, las=1, ylim = ylim,...)
      invisible(cbind(x,y))
    }
  else return(cbind(x,y))
}


################################################################################
##                                ZCR
################################################################################

zcr <- function(wave,
                f,
                wl = 512,
                ovlp = 0,
                plot = TRUE,
                type = "o",
                xlab = "Time (s)",
                ylab = "Zero crossing rate",
                ...)
    
    {
        
        ## INPUT
        input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; rm(input)
        n <- nrow(wave)

        ## function to compute the zero-crossing rate
        ZCR <- function(x)
            {
                sign <- ifelse(x>=0, 1, -1) 
                res <- sum(abs(diff(sign)))/(2*length(x))
                return(res)
            }
        ## first option (wl=NULL): computation of zcr on the whole signal
        if(is.null(wl)) {results <- ZCR(wave) ; return(results)}
        ## second option (wl = 2 exp(i)): computation of zcr on a sliding window and plot if requested
        else { 
            step <- seq(1, n-wl, wl-(ovlp*wl/100))
            res <- apply(as.matrix(step), MARGIN=1, function(x) ZCR(wave[x:(wl+x-1)]))
            t <- seq(0, n/f, length.out=length(step))
            results <- cbind(time=t, ZCR=res)
            if(plot)
                {
                    plot(x=t, y=res, type=type, xlab=xlab, ylab=ylab, ...)
                    invisible(results)
                }
            else return(results)                           
        }
    }




###########################
###########################
### ACCESSORY FUNCTIONS ###
###########################
###########################


################################################################################
##                                BARTLETT.W
################################################################################ 

bartlett.w<-function (n)
{
  if(n <= 0) stop("'n' must be a positive integer")

  n<-n-1
  m<-n%/%2
  w<-c((2*(0:(m-1)))/n, 2-((2*(m:n))/n))
  return(w)
}


################################################################################
##                                BLACKMAN.W
################################################################################ 

blackman.w<-function (n)
{
  if(n <= 0) stop("'n' must be a positive integer")
  n <- n-1
  w <- 0.42-0.5*cos(2*pi*(0:n)/n)+0.08*cos(4*pi*(0:n)/n)
  return(w)
}


################################################################################
##                       FILLED.CONTOUR.MODIF2
################################################################################ 
                                        # modification of filled.contour in graphics by Ross Ihaka

filled.contour.modif2<-function (x = seq(0, 1, len = nrow(z)),
                                 y = seq(0, 1, len = ncol(z)), z, xlim = range(x, finite = TRUE),
                                 ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE),
                                 levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors,
                                 col = color.palette(length(levels) - 1), plot.title, plot.axes, key.title,
                                 asp = NA, xaxs = "i", yaxs = "i", las = 1, axisX = TRUE, axisY = TRUE,...) 
{
  if(missing(z)) {
    if(!missing(x)) {
      if(is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq(0, 1, len = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if(is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if(any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  plot.new()
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  if(!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
    stop("no proper 'z' matrix specified")
  if(!is.double(z)) 
    storage.mode(z) <- "double"
  .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
                          col = col)
  if(missing(plot.axes))
    {
      if(axisX)
        {
          title(main="", xlab="",ylab="")
          axis(1)
        }
      if(axisY)
        {
          title(main="", xlab="",ylab="")
          axis(2)
        }
    }
  box()
  if(missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}


################################################################################
##                                FLATTOP.W
################################################################################ 

flattop.w<-function (n)
{
  if(n <= 0) stop("'n' must be a positive integer")
  n<-n-1
  w<-0.2156-0.4160*cos(2*pi*(0:n)/n)+0.2781*cos(4*pi*(0:n)/n)
  -0.0836*cos(6*pi*(0:n)/n)+0.0069*cos(8*pi*(0:n)/n)   
  return(w)
}


################################################################################
##                                HAMMING.W
################################################################################ 

hamming.w<-function (n)
{
  if(n <= 0) stop("'n' must be a positive integer")
  n<-n-1
  w<-0.54-0.46*cos(2*pi*(0:n)/n)
  return(w)
}


################################################################################
##                                HANNING.W
################################################################################ 

hanning.w<-function (n)
{
  if(n <= 0) stop("'n' must be a positive integer")
  n<-n-1
  w<-0.5-0.5*cos(2*pi*(0:n)/n)
  return(w)
}


################################################################################
##                                INPUTW
################################################################################ 

inputw <- function(
                  wave,
                  f,
                  channel = 1,
                  bit = 16
                  )
    
{
  if(is.data.frame(wave))   {f<-f ; wave <- as.matrix(wave[,channel])}
  if(is.vector(wave))       {f<-f ; wave <- as.matrix(wave)}
  ## WaveMC and mts objects are matrix by default, there is then a conflict between is.matrix and is.mts
  if(is.matrix(wave) && !is.mts(wave) && class(wave) != "WaveMC") {f<-f ; wave <- wave[,channel,drop=FALSE]}
  if(is.ts(wave))           {f<-frequency(wave) ; wave <- as.matrix(wave)} 
  if(is.mts(wave))          {f<-frequency(wave) ; wave <- as.matrix(wave[, channel])} 
  if(class(wave)=="Sample") {f<-wave$rate ; wave <- as.matrix(wave$sound[channel, ]) ; bit=wave@bits}
  if(class(wave)=="audioSample"){f<-wave$rate ; bit=wave$bits; wave <- as.matrix(wave)}
  if(class(wave)=="Wave")
    {
      f <- wave@samp.rate
      bit <- wave@bit
      if(channel==1) {wave <- as.matrix(wave@left)}   
      if(channel==2) {wave <- as.matrix(wave@right)}     
    }
  if(class(wave)=="WaveMC")
    {
      f <- wave@samp.rate  
      bit <- wave@bit
      wave <- as.matrix(wave@.Data[, channel])
    }
  return(list(w=wave, f=f, bit=bit))
}


################################################################################
##                                OUTPUTW
################################################################################ 

outputw <- function(
                    wave,
                    f,
                    bit = 16,
                    format
                    )
  {
    if(format == "matrix") {wave <- as.matrix(wave)}
    if(format == "Wave")  {wave <- Wave(left=wave, samp.rate=f, bit=bit)}
    if(format == "audioSample") {wave <- audio::audioSample(wave, rate=f, bits=bit, clip=FALSE)}
    if(format == "ts") {wave <- ts(wave, start=0, end=length(wave)/f, frequency=f)}
    return(wave)
  }



################################################################################
##                                RECTANGLE.W
################################################################################ 

rectangle.w<-function (n)
{
  if(n <= 0) stop("'n' must be a positive integer")
  w<-rep(1,n)
  return(w)
}




################################################################################
##                                RESCALE
################################################################################
## by Ethan Brown
## https://github.com/statisfactions/AutoregressionTones/blob/master/Autoregression.R
rescale <- function(x, lower, upper) {
  nrange <- upper-lower
  out <- ((x-min(x))*nrange/(max(x)-min(x)) - nrange/2)
}


################################################################################
##                                REVERSE.CM.COLORS
################################################################################
## reverse.cm.colors, reversion of cm.colors in grDevices package
## originally by R Development Core Team and contributors worldwide

reverse.cm.colors <- function (x)
{
  n <- x
  if((n <- as.integer(n[1])) > 0) {
    even.n <- n%%2 == 0
    k <- n%/%2
    l1 <- k + 1 - even.n
    l2 <- n - k + even.n
    rev(c(if(l1 > 0) hsv(h = 6/12, s = seq(0.5, ifelse(even.n, 
                                      0.5/k, 0), length = l1), v = 1), if(l2 > 1) hsv(h = 10/12, 
                                                                 s = seq(0, 0.5, length = l2)[-1], v = 1)))
  }
  else character(0)
}



################################################################################
##                                REVERSE.GRAY.COLORS.1
################################################################################ 
reverse.gray.colors.1 <- function (x) gray(seq(from = 1^1.7, to = 0, length = x)^(1/1.7))



################################################################################
##                                REVERSE.GRAY.COLORS.2
################################################################################ 
reverse.gray.colors.2 <- function (x)  gray(seq(from = 1, to = 0, length = x))



################################################################################
##                                REVERSE.HEAT.COLORS
################################################################################
## reverse.heat.colors, reversion of heat.colors in grDevices package
## originally by R Development Core Team and contributors worldwide 

reverse.heat.colors<-
  function (x) 
{
  n<-x
  if((n <- as.integer(n[1])) > 0) {
    j <- n%/%4
    i <- n - j
    rev(c(rainbow(i, start = 0, end = 1/6), if(j > 0) hsv(h = 1/6, 
                                  s = seq(from = 1 - 1/(2 * j), to = 1/(2 * j), length = j), 
                                  v = 1)))
  }
  else character(0)
}



################################################################################
##                        REVERSE.TERRAIN.COLORS
################################################################################
## reverse.terrain.colors, reversion of terrain.colors in grDevices package
## originally by R Development Core Team and contributors worldwide 

reverse.terrain.colors<-
  function (x)
{
  n<-x
  if((n <- as.integer(n[1])) > 0) {
    k <- n%/%2
    h <- c(4/12, 2/12, 0/12)
    s <- c(1, 1, 0)
    v <- c(0.65, 0.9, 0.95)
    rev(c(
          hsv(h = seq(h[1], h[2], length = k),
              s = seq(s[1], s[2], length = k),
              v = seq(v[1], v[2], length = k)),
          
          hsv(h = seq(h[2], h[3], length = n - k + 1)[-1],
              s = seq(s[2], s[3], length = n - k + 1)[-1],
              v = seq(v[2], v[3], length = n - k + 1)[-1])
          ))
  }
  else character(0)
}



################################################################################
##                                REVERSE.TOPO.COLORS
################################################################################
## reverse.topo.colors, reversion of topo.colors in grDevices package
## originally by R Development Core Team and contributors worldwide 

reverse.topo.colors<-
  function (x) 
{
  n<-x
  if((n <- as.integer(n[1])) > 0) {
    j <- n%/%3
    k <- n%/%3
    i <- n - j - k
    rev(c(if(i > 0) hsv(h = seq(from = 43/60, to = 31/60, length = i)), 
          if(j > 0) hsv(h = seq(from = 23/60, to = 11/60, 
                           length = j)), if(k > 0) hsv(h = seq(from = 10/60, 
                                                          to = 6/60, length = k), s = seq(from = 1, to = 0.3, 
                                                                                    length = k), v = 1)))
  }
  else character(0)
}



################################################################################
##                                SUMSMOOTH
################################################################################

sumsmooth <- function(x, wl, padding=TRUE, norm=FALSE){
  n <- length(x)
  WL <- wl-1
  x <- apply(as.matrix(1:(n-WL)), 1, function(y) sum(x[y:(y+WL)]))

  if(padding){
    x <- c(rep(0, times=wl/2), x, rep(0, times=wl/2))
    if(wl %% 2 == 0) {x <- x[-1]} # wl is even, remove one 0 at the start
  }

  if(norm) x <- x/max(x)
  return(x)
}


################################################################################
##                                SOSCILLO
################################################################################

soscillo<-function
(
 wave,
 f,
 from = FALSE,
 to =FALSE,
 colwave = "black",
 coltitle = "black",
 collab = "black",
 colline = "black",
 colaxis = "black",
 cexaxis = 1,
 coly0 = "grey47",
 tlab = "Times (s)",
 alab = "Amplitude",
 cexlab = 1,
 fontlab = 1,
 title = FALSE,
 xaxt="s",
 yaxt="n",
 tickup = NULL,
 ... 
 )

{
  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; rm(input)

  if(from|to)
    {
      if(from == 0) {a<-1; b<-round(to*f)}
      if(from == FALSE) {a<-1; b<-round(to*f);from<-0}
      if(to == FALSE) {a<-round(from*f); b<-nrow(wave);to<-nrow(wave)/f}
      else {a<-round(from*f); b<-round(to*f)}
      wave<-as.matrix(wave[a:b,])
      n<-nrow(wave)
    }
  else {n<-nrow(wave) ; from<-0 ; to<-n/f}

  par(tcl=0.5, col.axis=colaxis, cex.axis = cexaxis, col=colline, col.lab=collab,las=0)

  wave<-ts(wave[0:n,], start=from, end=to, frequency = f)

  plot(wave,
       col=colwave, type="l",
       xaxs="i", yaxs="i",
       xlab="", ylab="",
       xaxt=xaxt, yaxt=yaxt, bty="l",
       ...)
  axis(side=1, col=colline,labels=FALSE)
  axis(side=2, at=tickup, col=colline,labels=FALSE)

  mtext(tlab,col=collab,font=fontlab,cex=cexlab,side=1,line=3)
  mtext(alab,col=collab,font=fontlab,cex=cexlab,side=2,line=3)

  abline(h=0,col=coly0,lty=2)
}



################################################################################
##                                SSPECTRO
################################################################################

sspectro <- function
(
 wave,
 f,
 wl = 512,
 ovlp = 0,
 wn="hanning",
 norm=TRUE
 )

{
  input<-inputw(wave=wave,f=f) ; wave<-input$w ; f<-input$f ; rm(input)
  n<-nrow(wave)
  step<-seq(1,n-wl,wl-(ovlp*wl/100))
  W<-ftwindow(wl=wl,wn=wn)
  z <- apply(as.matrix(step), 1, function(x) Mod(fft(wave[x:(wl+x-1),]*W)))
  z<-z[2:(1+wl/2),]
  if(norm) {z <- z/max(z)}
  return(z)
}



################################################################################
##                                SPECTRO.COLORS
################################################################################

spectro.colors<-
  function (n)
{
  if((n <- as.integer(n[1])) > 0)
    {
      j <- n%/%3
      k <- n%/%3
      i <- n - j - k
      c(if(i > 0) hsv(h = seq(from = 31/60, to = 43/60, length = i), s = seq(0,1,length=i)),
        if(j > 0) hsv(h = seq(from = 21/60, to = 9/60, length = j), v = seq(0.5,0.8,length=j)),
        if(k > 0) hsv(h = seq(from = 8/60, to = 1/60, length = k), s = seq(from = 0.5, to = 1, length = k), v=1))
    }
  else character(0)
}

################################################################################
##                                 STFT
################################################################################

stft <- function(
    wave,
    f,
    wl,
    zp,
    step,
    wn,
    scale = TRUE,
    norm = FALSE,
    fftw = FALSE,
    complex = FALSE
    )

{
    zpl <- zp%/%2
    
    if(zpl==0)
        {
            W <- ftwindow(wl=wl,wn=wn)
            if(fftw)
                {
                    p <- fftw::planFFT(wl)
                    z <- apply(as.matrix(step), 1, function(x) fftw::FFT(wave[x:(wl+x-1)]*W, plan=p))
                }
            else {z <- apply(as.matrix(step), 1, function(x) fft(wave[x:(wl+x-1),]*W))}
        }
    
    else
        {
            W <- ftwindow(wl=wl+zp,wn=wn)
            if(fftw)
                {
                    p <- fftw::planFFT(wl+zp)
                    z <- apply(as.matrix(step), 1, function(x) fftw::FFT(c(1:zpl,wave[x:(wl+x-1),],1:zpl)*W, plan=p))
                }
            else {z <- apply(as.matrix(step), 1, function(x) fft(c(1:zpl,wave[x:(wl+x-1),],1:zpl)*W))}
        }

    z <- z[1:((wl+zp)%/%2), , drop=FALSE]	# to keep only the relevant frequencies (half of the FFT)

    if(complex == FALSE)
        {
            z <- Mod(z)
            
            if(scale)
                {
                    if(norm)  # normalise to 1 each column (ie each FFT), called by dynspec()
                        {
                            z1 <- matrix(numeric(length(z)))
                            dim(z1) <- dim(z)
                            for(i in 1:ncol(z1)) {z1[,i] <- z[,i]/max(z[,i])}
                            z <- z1
                            rm(z1)
                        }
                    else {z <- z/max(z)} # normalise to 1 the complete matrix
                }
        }
    return(z)
}




################################################################################
##                                 TEMP.COLORS
################################################################################

temp.colors<-function (n)
{
  if((n <- as.integer(n[1])) > 0) {
    j <- n%/%3
    k <- n%/%3
    i <- n - j - k
    c(
      if(i > 0) hsv(h=seq(from=44/60, to=31/60, length=i), s=seq(from=1, to=0.3, length=i), v=1),
      if(j > 0) hsv(h=seq(from=31/60, to=8/60,  length=j), s=seq(from=0.3, to=0.6, length=j), v=1),
      if(k > 0) hsv(h=seq(from= 8/60, to=1/60,  length=k), s=seq(from=0.6, to=1, length=k), v=1)
      )
  }
  else character(0)
}



################################################################################
##                                UNWRAP
################################################################################
## formerly in the signal package
## Original Octave version by Bill Lash. Conversion to R by Tom Short.

unwrap <- function (a, tol = pi, dim = 1) 
{
  sz = dim(a)
  nd = length(sz)
  if(nd == 0) {
    sz = length(a)
    nd = 1
  }
  if(!(length(dim) == 1 && dim == round(dim)) && dim > 0 && 
      dim < (nd + 1)) 
    stop("unwrap: dim must be an integer and valid dimension")
  while (dim < (nd + 1) && sz[dim] == 1) dim = dim + 1
  if(dim > nd) 
    dim = 1
  tol = abs(tol)
  rng = 2 * pi
  m = sz[dim]
  if(m == 1) 
    return(a)
  idx = list()
  for (i in 1:nd) idx[[i]] = 1:sz[i]
  idx[[dim]] = c(1, 1:(m - 1))
  d = a[unlist(idx)] - a
  p = rng * (((d > tol) > 0) - ((d < -tol) > 0))
  if(nd == 1) 
    r = cumsum(p)
  else r = apply(p, MARGIN = dim, FUN = cumsum)
  a + r
}
