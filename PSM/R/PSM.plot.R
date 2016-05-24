PSM.plot <- function(Data, Smooth = NULL, indiv = NULL, type = NULL) {
  isdata <- length(Data)
  for(i in 1:length(Data))
    isdata <- isdata && all(c('Y','Time') %in% names(Data[[i]]))
  if(!isdata)
    stop('Data does not appear to be valid.')
  if(!is.null(Smooth)) {
    issmooth <- length(Smooth)
    for(i in 1:length(Smooth))
      issmooth <- issmooth && all(c("Time","Xs","Ps","Xf","Pf",
                                    "Xp","Pp","Yp","R")
                                  %in% names(Smooth[[i]]))
    if(!issmooth)
      stop('Smooth does not appear to be valid.')
    if(length(Smooth)!=length(Data))
      stop('Smooth and Data has different length')
  }
  if(!is.null(indiv))
    if(!is.integer(indiv))
      stop('indiv must be a vector of integers')
  smoothplot <- function(type,obs,len,log) {
    if(is.null(Smooth))
       stop(paste('Smooth must be provied for plot type:',type))
    for(i in 1:len) {
      if(type!='Yp') {
        plot(Smooth[[j]]$Time,Smooth[[j]][[type]][i,],type="l",xlab="",ylab="",log=log)
      } else {
        plot(Data[[j]]$Time,na.omit(Smooth[[j]][[type]][i,]),type="l",
             xlab="",ylab="",log=log)
      }
      if(obs) {
        points(Data[[j]]$Time,Data[[j]]$Y[i,])
        rug(Data[[j]]$Time)
      }
      if(j==indiv[1])
        mtext(paste(type,i,sep=""),side=2,line=2.5)
      if(i==1 && k==1)
        title(main=paste('Subject',j))
    }
  }
  dataXplot <- function(type,log) {
    if(is.null(Data[[j]][[paste(type,'X',sep="")]]))
      stop(paste('Data does not contain ',type,'X',sep=""))
    for(i in 1:dimX) {
      plot(Data[[j]][[paste(type,'Time',sep="")]],
           Data[[j]][[paste(type,'X',sep="")]][i,],
           type="l",xlab="",ylab="",log=log)
      if(j==indiv[1])
        mtext(paste(type,'X',i,sep=""),side=2,line=2.5)
      if(i==1 && k==1)
        title(main=paste('Subject',j))
    }
  }
  dataplot <- function(log, type) { #plots Y and U
    if(is.null(Data[[j]][[type]]))
      stop(paste('Data does not contain ',type,sep=""))
    len <- dim(Data[[j]][[type]])[1]
    if(type=='Y') {
      len <- dimY
      plottype <- 'p'
    } else {
      len <- dimU
      plottype <- 's' #show zero-order hold on input.
    }
    for(i in 1:len) {
      plot(Data[[j]][['Time']], Data[[j]][[type]][i,], xlab="",
           ylab="", log=log, type=plottype)
      if(j==indiv[1])
        mtext(paste(type,i,sep=""),side=2,line=2.5)
      if(i==1 && k==1)
        title(main=paste('Subject',j))
    }
  }
  resfun <- function() {
    if(is.null(Smooth))
       stop(paste('Smooth must be provied for plot type: res and acf'))
    subs <- ceiling(length(Smooth[[j]]$Time)/length(Data[[j]]$Time)) - 1
    idx <- (1:dimT)*(subs+1)-subs
    Data[[j]]$Y-Smooth[[j]][['Yp']][,idx]
  }
  resplot <- function(log) {
    res <- resfun()
    for(i in 1:dimY) {
      plot(Data[[j]]$Time,res[i,],xlab="",ylab="",log=log)
      abline(h=0)
      if(j==indiv[1])
        mtext(paste('Y',i," - Yp",i,sep=""),side=2,line=2.5)
      if(i==1 && k==1)
        title(main=paste('Subject',j))
    }
  }
  acfplot <- function() {
    res <- resfun()
    for(i in 1:dimY) {
      tmpylab <- ifelse(j==indiv[1],paste('ACF (Y',i," - Yp",i,')',sep=""),'')
      acf(res[i,],main="",ylab=tmpylab)
      if(i==1 && k==1)
        title(main=paste('Subject',j))
    }
  }
  etaplot <- function() {
    sim <- ''
    if(is.null(Smooth))
      if(!is.null(Data[[j]]$eta)) {
        sim <- 'sim-'
      }
    if(sim=='') { #cannot use ifelse due to possible NULL
      tmp <- Smooth[[j]]$eta
    } else {
      tmp <- Data[[j]]$eta
    }
    plot.new()
    plot.window(xlim=c(0,10),ylim=c(0,10))
    len <- length(tmp)
    if(len==0) {
      mtext("No eta's",cex=.7)
    } else {
      for(m in 1:length(tmp))
        mtext(paste(sim,'eta',m,': ',signif(tmp[m],4),sep=""),line=-1*m+1,cex=.7)
    }
  }
  #Default plots
  if(is.null(type))
    if(is.null(Smooth)) {
      type <- c('Y')
    } else {
      type <- c('Xs','Ys.Y')
    }
  #Common vars
  dimS <- length(Data)
  if(is.null(Smooth)) {
    dimX <- ifelse(is.null(Data[[1]]$X),0,dim(Data[[1]]$X)[1])
  } else {
    dimX <- dim(Smooth[[1]]$Xs)[1]
  }
  dimY <- dim(Data[[1]]$Y)[1]
  dimU <- ifelse(is.null(Data[[1]]$U),0,dim(Data[[1]]$U)[1])
  if(is.null(indiv))
    indiv = 1:dimS
  numrows <- 0
  for(k in 1:length(type) )
    if('X' %in% unlist(strsplit(type[k],''))) {
      numrows <- numrows + dimX
    } else {
      if ('U' %in% unlist(strsplit(type[k],''))) {
        numrows <- numrows + dimU
      } else {
        numrows <- numrows + dimY
      }
    }
  #Setup plot
  par(mfcol=c(numrows,length(indiv)),mar=c(2,4,2,0)+.1)
  #Loop over individuals
  for(j in indiv) {
    dimT <- dim(Data[[j]]$Y)[2]
    for (k in 1:length(type)) {
      tk <- type[k]
      log <- ''
      for(i in 1:2) {
        if(substring(tk,1,5) %in% c('logx.','logy.')) {
          log <- paste(log,substring(tk,4,4),sep="")
          tk <- substring(tk,6)
        }
      }
      switch(tk,
             Xp = smoothplot(tk,obs=FALSE,len=dimX,log=log),
             Xf = smoothplot(tk,obs=FALSE,len=dimX,log=log),
             Xs = smoothplot(tk,obs=FALSE,len=dimX,log=log),
             Yp = smoothplot(tk,obs=FALSE,len=dimY,log=log),
             Ys = smoothplot(tk,obs=FALSE,len=dimY,log=log),
             Yp.Y = smoothplot('Yp',obs=TRUE,len=dimY,log=log),
             Ys.Y = smoothplot('Ys',obs=TRUE,len=dimY,log=log),
             X = dataXplot('',log=log),
             longX = dataXplot('long',log=log),
             Y = dataplot(log=log,type=tk),
             U = dataplot(log=log,type=tk),
             res = resplot(log=log),
             acf = acfplot(),
             eta = etaplot(),
             stop(paste('Unknown type:',tk))
             )       
    }
  }
}
