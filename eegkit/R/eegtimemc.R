eegtimemc <-
  function(time,voltmat,channel,size=c(0.75,0.75),
           vadj=0.5,hadj=0.5,xlab="",ylab="",voltSE=NULL,...){
    ###### Plots Multi-Channel EEG Time Courses
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: February 16, 2015
    
    voltmat <- as.matrix(voltmat)
    vdim <- dim(voltmat)
    if(length(time)!=vdim[1]){stop("Number of rows of voltmat must match length of time.")}
    if(length(channel)!=vdim[2]){stop("Number of columns of voltmat must match length of channel.")}
    if(length(size)!=2L){stop("Size must give size of each subplot (in inches).")}
    if(!is.null(voltSE[1])){
      voltSE <- as.matrix(voltSE)
      sedim <- dim(voltSE)
      if(any(c(sedim[1]!=vdim[1],sedim[2]!=vdim[2]))){stop("Input 'voltSE' must be matrix same size as input 'voltmat'.")}
      if(any(c(voltSE)<0)){stop("Input 'voltSE' must contain non-negative numbers.")}
    }
    if(any(size<=0)){stop("Incorrect input for size.")}
    vadj <- vadj[1]
    hadj <- hadj[1]
    if(vadj<=0){stop("Incorrect input for vadj.")}
    if(hadj<=0){stop("Incorrect input for hadj.")}
    eegcoord <- NULL
    data(eegcoord,envir=environment())
    chnames <- rownames(eegcoord)
    cidx <- match(channel,chnames)
    naidx <- which(is.na(cidx))
    if(length(naidx)>0){
      cidx <- cidx[-naidx]
      channel <- channel[-naidx]
      voltmat <- voltmat[,-naidx]
      vdim[2] <- (vdim[2]-length(naidx))
    }
    eegcoord <- eegcoord[cidx,]
    chnames <- chnames[cidx]
    rx <- range(eegcoord[,4])
    xcoord <- (eegcoord[,4]-rx[1])/(rx[2]-rx[1])
    ry <- range(eegcoord[,5])
    ycoord <- (eegcoord[,5]-ry[1])/(ry[2]-ry[1])
    
    # adapted from Frank E Harrell Jr's subplot function (Hmisc package) 
    oldpar <- par(no.readonly=TRUE)
    on.exit(par(oldpar))
    pin <- par("pin")
    cusr <- par("usr")
    cplt <- par("plt")
    for(j in 1:vdim[2]){
      tx <- (xcoord[j]-cusr[1])/(cusr[2]-cusr[1])
      ty <- (ycoord[j]-cusr[3])/(cusr[4]-cusr[3])
      xx <- c(tx-hadj*size[1]/pin[1], tx+(1-hadj)*size[1]/pin[1])
      yy <- c(ty-vadj*size[2]/pin[2], ty+(1-vadj)*size[2]/pin[2])
      fx <- xx*(cplt[2]-cplt[1])+cplt[1]
      fy <- yy*(cplt[4]-cplt[3])+cplt[3]
      if(j==1L){
        par(plt=c(fx,fy))
      } else {
        par(plt=c(fx,fy),new=TRUE)
      }
      if(is.null(voltSE[1])){
        eegtime(time,voltmat[,j],xlab=xlab,ylab=ylab,...)
      } else {
        eegtime(time,voltmat[,j],xlab=xlab,ylab=ylab,voltageSE=voltSE[,j],...)
      }
      mtext(channel[j],line=-0.5)
    }
        
  }