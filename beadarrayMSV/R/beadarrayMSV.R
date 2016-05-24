#New class-definitions
setClass('BeadSetIllumina',contains='eSet', prototype =
         prototype(new("VersionedBiobase",
                       versions = c(classVersion("eSet"), BeadSetIllumina = "1.0.0"))))
setMethod("initialize", "BeadSetIllumina",
          function(.Object,
                   assayData = assayDataNew(R = R, se.R = se.R, G = G, se.G = se.G,
                     no.beads = no.beads, storage.mode = 'list'),
                   R = new("matrix"), se.R = new("matrix"), G = new("matrix"),
                   se.G = new("matrix"), no.beads = new("matrix"), ...) {
            if (!missing(assayData) &&
                any(!missing(R),!missing(G),!missing(se.R),!missing(se.G),!missing(no.beads))){
              warning("Using 'assayData'; ignoring 'R', 'G', 'se.R', 'se.G', 'no.beads'")
            }
            callNextMethod(.Object, assayData = assayData, ...)
          })
setValidity("BeadSetIllumina", function(object) {
  val <- assayDataValidMembers(assayData(object), c("R", "se.R", "G", "se.G", "no.beads"))
  if (is.logical(val)){
    val <- storageMode(object) %in% "list"
    if (!val){
      val <- "\nRequired storageMode: 'list'"
    }
  }
  val
})

setClass('AlleleSetIllumina',contains='eSet', prototype =
         prototype(new("VersionedBiobase",
                       versions = c(classVersion("eSet"), AlleleSetIllumina = "1.0.0"))))
setMethod("initialize", "AlleleSetIllumina",
          function(.Object, intensity = new("matrix"), theta = new("matrix"),
                   SE = new("matrix"), ...) {
            callNextMethod(.Object, intensity = intensity, theta = theta, SE = SE,
                           storage.mode = 'list', ...)
          })
setValidity("AlleleSetIllumina", function(object) {
  val <- assayDataValidMembers(assayData(object), c("intensity", "theta", "SE"))
  if (is.logical(val)){
    val <- storageMode(object) %in% "list"
    if (!val){
      val <- "\nRequired storageMode: 'list'"
    }
  }
  val
})



#Reads iScan's summary files. Uses code from 'readIllumina()' {beadarray}
readBeadSummaryOutput <- function(arrayNames=NULL,path='.',pattern='beadTypeFile.txt',recursive=FALSE,sep=',',fullPaths=NULL,sepchar='_',prList=NULL){
  #Locate files
  if (is.null(fullPaths))
    fullPaths <- dir(path=path,pattern=pattern,recursive=recursive)
  if (length(fullPaths) == 0) 
    stop(paste("No files with extension", pattern, "found"))
  tmp <- unlist(strsplit(fullPaths,'/'))
  nRec <- length(tmp)/length(fullPaths)
  indNames <- vector(length=length(tmp))
  indNames[seq(nRec,length(tmp),nRec)] <- TRUE
  rFiles <- tmp[indNames]
  rPaths <- rep('',length(fullPaths))
  if (length(grep('/',fullPaths)) == length(fullPaths)){
    if (nRec>1){
      for (i in 1:(nRec-1)){
        indi <- seq(i,length(tmp),nRec)
        rPaths <- paste(rPaths,tmp[indi],'/',sep='')
      }
    }
  }
  arrays <- strtrim(rFiles, nchar(rFiles) - nchar(pattern) - 1)
  if (!is.null(arrayNames)) {
    arrayNames <- gsub(paste('_',pattern,sep=''), "", arrayNames)
    ind <- which(arrays %in% arrayNames)
    if (length(ind) == 0) 
      stop("'arrayNames' did not match with files in specified path(s)")
    arrays <- arrays[ind]
    rFiles <- rFiles[ind]
    rPaths <- rPaths[ind]
  }
  nArrays <- length(arrays)
  message('Found ', nArrays, ' arrays')
  
  #Prepare phenoData
  arrayInfo <- data.frame(arrayNames=as.character(arrays),no.beads=rep(0,nArrays),stringsAsFactors=FALSE,row.names=arrays)
  if (length(grep(sepchar,arrays)) == nArrays){
    tmp <- unlist(strsplit(arrays,sepchar))
    if (length(tmp) == nArrays*3){
      arrayInfo$chip <- tmp[seq(1, length(tmp), by=3)]
      arrayInfo$row <- tmp[seq(2, length(tmp), by=3)]
      arrayInfo$col <- tmp[seq(3, length(tmp), by=3)]
    }else{
      arrayInfo$chip <- tmp[seq(1, length(tmp), by=2)]
      arrayInfo$row <- tmp[seq(2, length(tmp), by=2)]
      arrayInfo$col <- rep(1, length(arrayInfo$chip))
    }
  }
  labelDesc <- c('Array ID','No. detected beads','Chip ID','Array ID within chip','Strip ID within array')
  phenoMetadata <- data.frame(labelDescription=labelDesc,row.names=colnames(arrayInfo))
  
  #Read assayData
  filesFull <- file.path(path, paste(rPaths,rFiles,sep=""))
  fc <- file(filesFull[1],open="r")
  header <- strsplit(readLines(fc,n=1),sep)
  if (!all(unlist(header)==c("Illumicode","N","Mean GRN","Dev GRN","Mean RED","Dev RED")))
    stop(paste('Unknown fields in header:',header))
  dat <- scan(file=fc,sep=sep,quiet=TRUE,
              what=list(ProbeID=integer(0),no.beads=integer(0),G=integer(0),sd.G=integer(0),
                R=integer(0),sd.R=integer(0)))
  close(fc)
  if (is.null(prList)){
    prList <- dat$ProbeID
  }
  if (length(prList)>length(dat$ProbeID)){
    stop(paste('The array',arrays[1],'contains fewer probes than specified in "prList"'))
  }else if (length(prList)<length(dat$ProbeID)){
    warning(paste('The array',arrays[1],'contains more probes than specified in "prList"'))
  }
  indP <- dat$ProbeID %in% prList #sapply(prList,function(x,ProbeID) which(ProbeID%in%x)[1],dat$ProbeID)
  if (!identical(prList,dat$ProbeID[indP])){
    stop(paste('The markers in array',arrays[1],'deviate from those in "prList"\n or are sorted differently.'))
  }
  nProbes <- length(prList)
  G <- se.G <- R <- se.R <- no.beads <-
    matrix(nrow=nProbes,ncol=nArrays,dimnames=list(prList,arrays))
  G[,1] <- dat$G[indP]
  se.G[,1] <- dat$sd.G[indP]/sqrt(dat$no.beads[indP])
  R[,1] <- dat$R[indP]
  se.R[,1] <- dat$sd.R[indP]/sqrt(dat$no.beads[indP])
  no.beads[,1] <- dat$no.beads[indP]
  if (nArrays>1){
    for (i in 2:nArrays){
      fc <- file(filesFull[i],open="r")
      dat <- scan(file=fc,sep=sep,quiet=TRUE,skip=1,
        what=list(ProbeID=integer(0),no.beads=integer(0),G=integer(0),sd.G=integer(0),
          R=integer(0),sd.R=integer(0)))
      close(fc)
      if (identical(prList,dat$ProbeID)){
        indP <- rep(TRUE,nProbes)
      }else{
        indP <- dat$ProbeID %in% prList
        if (!identical(prList,dat$ProbeID[indP])){
          stop(paste('Different probe-sets!\n',length(prList),'probes in "prList",',length(dat$ProbeID),'probes in array',arrays[i]))
        }
        warning(paste('The array',arrays[i],'contains more probes than specified in "prList"'))
      }
      G[,i] <- dat$G[indP]
      se.G[,i] <- dat$sd.G[indP]/sqrt(dat$no.beads[indP])
      R[,i] <- dat$R[indP]
      se.R[,i] <- dat$sd.R[indP]/sqrt(dat$no.beads[indP])
      no.beads[,i] <- dat$no.beads[indP]
    }
  }
  arrayInfo$no.beads <- apply(no.beads,2,sum)
  phenoData = new("AnnotatedDataFrame", data=arrayInfo, varMetadata=phenoMetadata)
  featureInfo <- data.frame(Address=prList,row.names=prList)
  featureMetadata <- data.frame(labelDescription="Unique strand identifier",row.names=colnames(featureInfo))
  featureData = new("AnnotatedDataFrame", data=featureInfo, varMetadata=featureMetadata)
  
  #Create summary data
  BSData = new('BeadSetIllumina', R=R, se.R=se.R, G=G, se.G=se.G, no.beads=no.beads,
    phenoData=phenoData, annotation='illuminaProbeIDs',
    featureData=featureData)
  sampleNames(BSData) <- make.names(sampleNames(BSData))
  featureNames(BSData) <- make.names(featureNames(BSData))
  rm(G,R,se.G,se.R,no.beads,arrayInfo,featureInfo,phenoData,featureData); gc()
  BSData
}



#Retrieves indexes to probes corresponding to specified Norm.ID's
getNormInd <- function(beadInfo,featureNames,normID=NULL,verbose=TRUE){
  if (is.null(normID))
    normID <- unique(beadInfo$Norm.ID)   #use all
  normInd <- matrix(FALSE,nrow=length(featureNames),ncol=length(normID),
                    dimnames=list(sub('X','',featureNames),normID))
  for (i in 1:length(normID)){
    indi <- beadInfo$Norm.ID==normID[i]
    probeIDs <- beadInfo$Address[indi]
    if (nchar(normID[i])>1)
      probeIDs <- c(probeIDs,beadInfo$Address2[indi])
    normInd[as.character(probeIDs),i] <- TRUE
  }
  rownames(normInd) <- featureNames
  if (verbose){
    cat('Sum Norm.IDs:\n')
    print(apply(normInd,2,sum))
  }
  normInd
}



#Plots A vs. B, or R vs. G, depending on BSData
scatterArrays <- function(BSData,arrays,markers=seq(1,nrow(BSData),5),smooth=TRUE,newFigure=TRUE,maxPlots=72,...){
  an <- sampleNames(BSData)[arrays]
  nSelect <- length(an)
  if (newFigure){
    if (nSelect>maxPlots)
      stop("Number of arrays for plotting exceeds 'maxPlots'")
    pp <- par()
    nC <- ceiling(sqrt(nSelect))
    nR <- ceiling(nSelect/nC)
    par(mfrow=c(nR,nC),mai=c(1.2,1,1,.5)/sqrt(nR),...)
  }else if (nSelect>1){
    stop("Only one array allowed when 'newFigure'==FALSE")
  }
  if (all(c('A','B') %in% names(assayData(BSData)))){
    x <- 'A'
    y <- 'B'
  }else if (all(c('R','G') %in% names(assayData(BSData)))){
    x <- 'R'
    y <- 'G'
  }else{
    stop('Cartesian coordinate intensities required in assayData.')
  }
  for (i in 1:nSelect){
    if (smooth){
      #require('geneplotter')
      smoothScatter(assayData(BSData)[[x]][markers,an[i]],assayData(BSData)[[y]][markers,an[i]],main=an[i],xlab=x,ylab=y)
    }else{
      plot(assayData(BSData)[[x]][markers,an[i]],assayData(BSData)[[y]][markers,an[i]],main=an[i],xlab=x,ylab=y,pch='.')
    }
    xy <- range(par('usr'))
    lines(xy,xy,col='green')
    lines(xy,c(0,0),col='green')
    lines(c(0,0),xy,col='green')
  }
  if (newFigure){
    par(pp[c('mfcol','mai')])
  }
}



#Rotate and shear BEFORE any other transf.
shearRawSignal <- function(BSData,plot=FALSE,newFigure=plot,normOpts=setNormOptions(),maxPlots=72,...){
  nArrays <- ncol(BSData)
  if (newFigure){
    if (nArrays>maxPlots)
      stop("Number of arrays for plotting exceeds 'maxPlots'")
    pp <- par()
    nC <- ceiling(sqrt(nArrays))
    nR <- ceiling(nArrays/nC)
    par(mfrow=c(nR,nC),mai=c(1,1,1,.5)/sqrt(nC))
  }else if (plot & nArrays>1){
    stop("Only one array allowed when 'plot'==TRUE and 'newFigure'==FALSE")
  }
  message('Rotating and shearing raw signal...')
  for (i in 1:nArrays){
    rawData <- list()
    rawData$x <- assayData(BSData)$R[,i]
    rawData$y <- assayData(BSData)$G[,i]
    trData <- normalizeIllumina(rawData,normOpts=normOpts,plot=plot,...)
    assayData(BSData)$R[,i] <- trData$x
    assayData(BSData)$G[,i] <- trData$y
    gc()
  }
  if (newFigure) par(pp[c('mfcol','mai')])
  validObject(BSData)
  BSData
}



#If 'indTrain' is a list containing indexes to red and green infinium-I
#probes, lines are fitted through these data-points only. Otherwise, same
#probes are used for finding both red and grean shear. Default: all probes.
#NB! Be aware of how xBase and yBase is calculated
normalizeIllumina <- function(rawData,indTrain=rep(TRUE,length(rawData$x)),normOpts=setNormOptions(),plot=FALSE,xylim=NULL,verbose=FALSE){
  if (is.list(indTrain)){ #2-element list with indexes to infinium-I probes
    xTrainR <- rawData$x[indTrain$'10r']
    xTrainG <- rawData$x[indTrain$'20g']
    yTrainR <- rawData$y[indTrain$'10r']
    yTrainG <- rawData$y[indTrain$'20g']
    iInorm <- TRUE
  }else{
    xTrainR <- xTrainG <- rawData$x[indTrain]
    yTrainR <- yTrainG <- rawData$y[indTrain]
    iInorm <- FALSE
  }
    
  #Fit line through AA
  rngx <- range(xTrainR,na.rm=TRUE)
  step <- diff(rngx)/(normOpts$nBins-1)
  breaks <- seq(rngx[1],rngx[2],step)
  midsx <- breaks[-normOpts$nBins] + step/2
  xBase <- rep(NA,length=normOpts$nBins-1)
  for (i in 1:(normOpts$nBins-1)){
    indi <- xTrainR>breaks[i] & xTrainR<breaks[i+1]
    if (sum(indi,na.rm=T)>normOpts$minSize){
      if (iInorm){
        xBase[i] <- quantile(yTrainR[indi],prob=normOpts$prob,na.rm=TRUE)  #Best for inf-I
      }else{
        xBase[i] <- quantile(sort(yTrainR[indi])[1:normOpts$minSize],prob=normOpts$prob,na.rm=TRUE)
      }
    }
  }
  if (sum(!is.na(xBase))>2){
    abx <- coef(lm(xBase~midsx))
  }else{
    abx <- rep(NA,2)
    names(abx) <- c('(Intercept)','midsx')
  }
  
  #Fit line through BB 
  rngy <- range(yTrainG,na.rm=TRUE)
  step <- diff(rngy)/(normOpts$nBins-1)
  breaks <- seq(rngy[1],rngy[2],step)
  midsy <- breaks[-normOpts$nBins] + step/2
  yBase <- rep(NA,length=normOpts$nBins-1)
  for (i in 1:(normOpts$nBins-1)){
    indi <- yTrainG>breaks[i] & yTrainG<breaks[i+1]
    if (sum(indi,na.rm=T)>normOpts$minSize){
      if (iInorm){
        yBase[i] <- quantile(xTrainG[indi],prob=normOpts$prob,na.rm=TRUE)
      }else{
        yBase[i] <- quantile(sort(xTrainG[indi])[1:normOpts$minSize],prob=normOpts$prob,na.rm=TRUE)
      }
    }
  }
  if (sum(!is.na(yBase))>2){
    aby <- coef(lm(yBase~midsy))
    aby <- c(-aby[1]/aby[2],1/aby[2])
  }else{
    aby <- rep(NA,2)
    names(aby) <- c('(Intercept)','midsy')
  }
  
  #Rotate and shear 
  normData <- rawData
  normData$x0 <- (aby[1]-abx[1])/(abx[2]-aby[2])
  normData$y0 <- abx[1] + abx[2]*normData$x0
  normData$x <- normData$x - normData$x0
  normData$y <- normData$y - normData$y0
  xx <- cbind(c(1,1),c(0,1))
  xAngle <- atan2(diff(xx%*%abx),1)
  yAngle <- atan2(diff(xx%*%aby),1)
  tmpX2 <- cos(xAngle)*normData$x + sin(xAngle)*normData$y
  normData$y <- -sin(xAngle)*normData$x + cos(xAngle)*normData$y
  normData$x <- tmpX2-normData$y/tan(yAngle-xAngle)
  
  #Plot input data with fitted lines
  if (plot){
    if (is.null(xylim)){
      plot(rawData$x,rawData$y,pch='.',main=paste('nBins=',normOpts$nBins,', prob quantile=',normOpts$prob),xlab='R',ylab='G')
      if (iInorm){
        points(xTrainR,yTrainR,pch='.',col='red')
        points(xTrainG,yTrainG,pch='.',col='green')
      }else{
        points(rawData$x[-which(indTrain)],rawData$y[-which(indTrain)],pch='.',col='red')
      }
      xylim <- par('usr')
    }else{
      plot.new()
      plot.window(xylim[1:2],xylim[3:4])
      axis(1); axis(2)
      points(rawData$x,rawData$y,pch='.',main=paste('nBins=',normOpts$nBins,', prob quantile=',normOpts$prob),xlab='R',ylab='G')
      if (iInorm){
        points(xTrainR,yTrainR,pch='.',col='red')
        points(xTrainG,yTrainG,pch='.',col='green')
      }else{
        points(rawData$x[-which(indTrain)],rawData$y[-which(indTrain)],pch='.',col='red')
      }
   }
    xy <- range(xylim)
    points(midsx,xBase,col='blue',pch=16)
    lines(xy,abx[1]+xylim[1:2]*abx[2],col='red')
    points(yBase,midsy,col='blue',pch=16)
    lines(xy,aby[1]+xy*aby[2],col='green')
    #points(x0,y0,col='blue',pch=16)
    lines(xy,rep(normData$y0,2),col='blue')
    lines(rep(normData$x0,2),xy,col='blue')
  }
  if (verbose){
    message(sprintf('\tOffset\tSlope'))
    message(sprintf('R:\t%.4g\t%.4g',abx[1],abx[2]))
    message(sprintf('G:\t%.4g\t%.4g',aby[1],aby[2]))
    #message(sprintf('%.4g\t',aby))
  }
  normData
}



#Returns a matrix of parametrized noise for each array
getNoiseDistributions <- function(BSData,subBeadPool=NULL,normInd,normOpts=setNormOptions(),plot=FALSE,newFigure=plot,maxPlots=72,xlim=NULL,...){
  nArrays <- ncol(BSData)
  if (newFigure){
    if (nArrays>maxPlots)
      stop("Number of arrays for plotting exceeds 'maxPlots'")
    pp <- par()
    nC <- ceiling(sqrt(nArrays))
    nR <- ceiling(nArrays/nC)
    par(mfrow=c(nR,nC),mai=c(1,1,1,.5)/sqrt(nC))
  }else if (plot & nArrays>1){
    stop("Only one array allowed when 'plot'==TRUE and 'newFigure'==FALSE")
  }
  an <- sampleNames(BSData)
  
  sPools <- colnames(normInd)
  if (is.null(subBeadPool))
    subBeadPool <- sPools[nchar(sPools)==1]
  #b10x <- normInd[,'101']
  #b20x <- normInd[,'201']
  #if (length(subBeadPool)>1){
    #for (i in 2:length(subBeadPool)){
      #b10x <- b10x | normInd[,as.character(100+i)]
      #b20x <- b20x | normInd[,as.character(200+i)]
    #}
  #}
  b10x <- normInd[,as.character(100+as.numeric(subBeadPool[1]))]
  b20x <- normInd[,as.character(200+as.numeric(subBeadPool[1]))]
  for (i in subBeadPool[-1]){
    b10x <- b10x | normInd[,as.character(100+as.numeric(i))]
    b20x <- b20x | normInd[,as.character(200+as.numeric(i))]
  }
  
  noiseDist <- matrix(nrow=nArrays,ncol=4,dimnames=list(an,c('medianR','medianG','madR','madG')))
  trChannel <- transformChannels(assayData(BSData)$R[b20x,],assayData(BSData)$G[b10x,],normOpts=normOpts)
  noiseDist[,'medianR'] <- apply(trChannel$X,2,median,na.rm=TRUE)
  noiseDist[,'medianG'] <- apply(trChannel$Y,2,median,na.rm=TRUE)
  noiseDist[,'madR'] <- apply(trChannel$X,2,mad,na.rm=TRUE)
  noiseDist[,'madG'] <- apply(trChannel$Y,2,mad,na.rm=TRUE)
  if (plot){
    if (is.null(xlim)){
      meanD <- apply(noiseDist,2,mean,na.rm=TRUE)
      xlim <- c(floor(min(meanD[1:2]-4*max(meanD[3:4]))),ceiling(max(meanD[1:2]+4*max(meanD[3:4]))))
    }
    for (i in 1:nArrays){
      hx <- hist(trChannel$X[,i],breaks=normOpts$breaks,plot=FALSE)
      hy <- hist(trChannel$Y[,i],breaks=normOpts$breaks,plot=FALSE)
      plot(hx$mids,hx$density,type='l',main=paste('Noise level,',an[[i]]),xlab='',ylab='',col='gray',xlim=xlim,...)
      lines(hx$mids,dnorm(hx$mids,noiseDist[i,1],noiseDist[i,3]),col='red')
      lines(hy$mids,hy$density,col='gray')
      lines(hy$mids,dnorm(hy$mids,noiseDist[i,2],noiseDist[i,4]),col='seagreen')
    }
  }
  if (newFigure){par(pp[c('mfcol','mai')])}
  noiseDist
}



#Transform BSData
transformChannels <- function(X,Y=NULL,normOpts=setNormOptions()){
  str1 <- '%d of %d %s-intensities below %d set to NA'
  switch(normOpts$transf,
         none={lstr=' '},
         log={
           X[X<=-normOpts$offset] <- NA
           message(sprintf(str1,sum(is.na(X)),length(X),'X',-normOpts$offset))
           X <- log2(X+normOpts$offset)
           lstr <- 'Log2'
           if (normOpts$offset!=0)
             lstr <- paste(lstr,normOpts$offset,'+')
           if (!is.null(Y)){
             Y[Y<=-normOpts$offset] <- NA
             message(sprintf(str1,sum(is.na(Y)),length(Y),'Y',-normOpts$offset))
             Y <- log2(Y+normOpts$offset)
           }},
         root={
           X[X<=-normOpts$offset] <- NA
           message(sprintf(str1,sum(is.na(X)),length(X),'X',-normOpts$offset))
           X <- (X+normOpts$offset)^(1/normOpts$nthRoot)
           lstr <- paste('nth-root (n=',normOpts$nthRoot,')',sep='')
           if (normOpts$offset!=0)
             lstr <- paste(lstr,', ',normOpts$offset,' +',sep='')
           if (!is.null(Y)){
             Y[Y<=-normOpts$offset] <- NA
             message(sprintf(str1,sum(is.na(Y)),length(Y),'Y',-normOpts$offset))
             Y <- (Y+normOpts$offset)^(1/normOpts$nthRoot)
           }}
         )
  trChannel <- list()
  trChannel$X <- as.matrix(X)
  trChannel$Y <- as.matrix(Y)
  trChannel$lstr <- lstr
  trChannel
}



#Scatter-plot with origin and noise indicated
#NB! 'noiseDist' must be transformed according to 'normOpts'. normInd = single index-vector
plotEstimatedNoise <- function(BSData,noiseDist,normInd=rep(TRUE,nrow(BSData)),normOpts=setNormOptions(),newFigure=TRUE,maxPlots=72,...){
  if (!identical(sampleNames(BSData),rownames(noiseDist)))
    stop("'BSData' and 'noiseDist' contain different samples")
  if (newFigure){
    if (ncol(BSData)>maxPlots)
      stop("Number of arrays for plotting exceeds 'maxPlots'")
    pp <- par()
    nC <- ceiling(sqrt(ncol(BSData)))
    nR <- ceiling(ncol(BSData)/nC)
    par(mfrow=c(nR,nC),mai=c(1,1,1,.5)/sqrt(nC))
  }else if (ncol(BSData)>1){
    stop("Only one array allowed when 'newFigure'==FALSE")
  }
  an <- sampleNames(BSData)
  for (i in 1:ncol(BSData)){
    trData <- transformChannels(assayData(BSData)$R[normInd,an[i]],assayData(BSData)$G[normInd,an[i]],normOpts=normOpts)
    plot(trData$X,trData$Y,pch='.',main=paste(an[i],', nSD =',normOpts$nSD),...)
    xy <- range(par('usr'))
    lines(rep(noiseDist[an[i],'medianR'],2),xy,col='green')
    lines(xy,rep(noiseDist[an[i],'medianG'],2),col='green')
    points(c(noiseDist[an[i],'medianR']-2*noiseDist[an[i],'madR'],
             noiseDist[an[i],'medianR']+2*noiseDist[an[i],'madR']),
           rep(noiseDist[an[i],'medianG'],2),type='b',pch=16,col='red')
    points(rep(noiseDist[an[i],'medianR'],2),
           c(noiseDist[an[i],'medianG']-2*noiseDist[an[i],'madG'],
             noiseDist[an[i],'medianG']+2*noiseDist[an[i],'madG']),type='b',pch=16,col='red')
  }
  if (newFigure) par(pp[c('mfcol','mai')])
}



#Transform SE's using a first order Taylor-expansion around the mean
transformSEs <- function(X,se.X,normOpts=setNormOptions()){
  str1 <- '%d of %d intensities below %d set to NA'
  switch(normOpts$transf,
         none={lstr=' '},
         log={
           se.X[X<=-normOpts$offset] <- NA
           message(sprintf(str1,sum(is.na(se.X)),length(se.X),-normOpts$offset))
           se.X <- se.X/(X+normOpts$offset)/log(2)
           lstr <- 'Log2'
           if (normOpts$offset!=0){
             warning('transformSEs() not properly tested for transf==\'log\' and offset!=0')
             lstr <- paste(lstr,normOpts$offset,'+')
           }
         },
         root={
           se.X[X<=-normOpts$offset] <- NA
           message(sprintf(str1,sum(is.na(se.X)),length(se.X),-normOpts$offset))
           se.X <- se.X*(X+normOpts$offset)^((1-normOpts$nthRoot)/normOpts$nthRoot)/normOpts$nthRoot
           lstr <- paste('nth-root (n=',normOpts$nthRoot,')',sep='')
           if (normOpts$offset!=0)
             lstr <- paste(lstr,', ',normOpts$offset,' +',sep='')
         })
  trSE <- list()
  trSE$SE <- se.X
  trSE$lstr <- lstr
  trSE
}



#Centres and scales based on percentile, assumes already rotated and sheared
normalizeShearedChannels <- function(trChannel,noiseDist,normOpts=setNormOptions()){#,b101=NULL,b201=NULL){
  nProbes <- nrow(trChannel$X)
  trChannel$X0 <- noiseDist[,'medianR']
  trChannel$Y0 <- noiseDist[,'medianG']
  trChannel$X.SE <- noiseDist[,'madR']
  trChannel$Y.SE <- noiseDist[,'madG']
  trChannel$X <- trChannel$X - matrix(1,nProbes,1)%*%trChannel$X0
  trChannel$Y <- trChannel$Y - matrix(1,nProbes,1)%*%trChannel$Y0
  
  switch(normOpts$method,
         none={trChannel$method <- 'The channels are not normalized'},
         quantNorm={   #Quantile normalization
           #require('limma')
           #iX <-  apply(trChannel$X,2,function(x) all(!is.na(x)))
           #iY <-  apply(trChannel$Y,2,function(x) all(!is.na(x)))
           iX <-  apply(trChannel$X,2,function(x) sum(!is.na(x))>1)
           iY <-  apply(trChannel$Y,2,function(x) sum(!is.na(x))>1)
           for (i in which(iX & iY)){
             normXY <- normalizeQuantiles(cbind(trChannel$X[,i],trChannel$Y[,i]))
             trChannel$X[,i] <- normXY[,1]
             trChannel$Y[,i] <- normXY[,2]
           }
           #trChannel$X.SE <- trChannel$Y.SE <- sqrt(trChannel$X.SE^2+trChannel$Y.SE^2)/2   #NB!
           trChannel$method <- 'The data have been quantile normalized (SE\'s not modified)'
         },
         medianAF={  #Scales red channel such that median(AF) is close to 0.5
           message('Median(AF)-estimates for a few arrays:')
           nDisp <- min(4,ncol(trChannel$X))
           message(sprintf('%s\t',colnames(trChannel$X)[1:nDisp]))
           for (rep in 1:8){
             PAF <- trChannel$X/(trChannel$X+trChannel$Y)
             medianAF <- apply(PAF,2,median,na.rm=TRUE)
             message(sprintf('%.6f\t',medianAF[1:nDisp]))
             trChannel$X <- trChannel$X/matrix(1,nProbes,1)%*%medianAF/2
             trChannel$X.SE <- trChannel$X.SE/medianAF/2
             #gc()
           }
           trChannel$method <- 'The red channel scaled such that median(AF) is close to 0.5'
         },
         linPeak={  #Scales each channel by its scale-th quantile
           trChannel$XScale <- apply(trChannel$X,2,quantile,normOpts$scale,na.rm=TRUE)
           trChannel$YScale <- apply(trChannel$Y,2,quantile,normOpts$scale,na.rm=TRUE)
           trChannel$X <- trChannel$X/matrix(1,nProbes,1)%*%trChannel$XScale
           trChannel$Y <- trChannel$Y/matrix(1,nProbes,1)%*%trChannel$YScale
           trChannel$X.SE <- trChannel$X.SE/trChannel$XScale
           trChannel$Y.SE <- trChannel$Y.SE/trChannel$YScale
           trChannel$method <- paste('Both channels scaled using the ',normOpts$scale,'-th quantile',sep='')
         },
         #linHomo={  #Scales each channel using quantile of homozygotes
         #  stop('linHomo not implemented for matrices')
         #  trChannel$XScale <- quantile(trChannel$X[b101],normOpts$scale,na.rm=TRUE)
         #  trChannel$YScale <- quantile(trChannel$Y[b201],normOpts$scale,na.rm=TRUE)
         #  trChannel$X <- trChannel$X/trChannel$XScale
         #  trChannel$Y <- trChannel$Y/trChannel$YScale
         #  trChannel$X.SE <- trChannel$X.SE/trChannel$XScale
         #  trChannel$Y.SE <- trChannel$Y.SE/trChannel$YScale
         #  trChannel$method <- paste('Each channel scaled using the ',normOpts$scale,'-th quantile of its infinium I homozygotes',sep='')
         #}
  )
  trChannel
}



#Merge BSData (R/G) into BSRed (A/B)
#beadInfo must contain 'Name', 'SNP', 'ILMN.Strand', 'Address', 'Address2', 'Norm.ID'
createAlleleSet <- function(BSData,beadInfo,normOpts,includeAB=FALSE){
  nSNP <- nrow(beadInfo)
  sPools <- unique(beadInfo$Norm.ID)
  subBeadPool <- sPools[nchar(sPools)==1]
  A <- B <- Var.pooled <- Min.Beads <- 
    matrix(nrow=nSNP,ncol=ncol(BSData),dimnames=list(rownames(beadInfo),sampleNames(BSData)))
  
  #Addresses are sorted ascendingly
  allProbes <- featureData(BSData)[[1]]
  indA1 <- which(is.element(allProbes,beadInfo$Address))
  indA2 <- which(is.element(allProbes,beadInfo$Address2))
  #Addresses ordered to follow beadInfo
  sortAdd1 <- sort(beadInfo$Address,index.return=TRUE)
  unSortAdd1 <- sort(sortAdd1$ix,index.return=TRUE)  #same as 'order()'
  Add2 <- beadInfo$Address2
  Add2[Add2==0] <- NA
  sortAdd2 <- sort(Add2,index.return=TRUE)
  unSortAdd2 <- sort(sortAdd2$ix[1:length(indA2)],index.return=TRUE)
  #Indexes used with 'Address'
  norm1 <- norm101 <- norm201 <- rep(FALSE,nSNP)
  for (i in subBeadPool){
    norm1 <- norm1 | beadInfo$Norm.ID==as.character(i)
    norm101 <- norm101 | beadInfo$Norm.ID==as.character(100+i)
    norm201 <- norm201 | beadInfo$Norm.ID==as.character(200+i)
  }
  #Indexes used with 'Address2'
  norm101_2 <- norm201_2 <- rep(FALSE,length(indA2))
  for (i in subBeadPool){
    norm101_2 <- norm101_2 | beadInfo[norm201|norm101,'Norm.ID']==as.character(100+i)
    norm201_2 <- norm201_2 | beadInfo[norm201|norm101,'Norm.ID']==as.character(200+i)
  }

  A[norm1,] <- assayData(BSData)$R[indA1,][unSortAdd1$ix,][norm1,]
  B[norm1,] <- assayData(BSData)$G[indA1,][unSortAdd1$ix,][norm1,]
  A[norm101,] <- assayData(BSData)$R[indA1,][unSortAdd1$ix,][norm101,]
  A[norm201,] <- assayData(BSData)$G[indA1,][unSortAdd1$ix,][norm201,]
  B[norm101,] <- assayData(BSData)$R[indA2,][unSortAdd2$ix,][norm101_2,]
  B[norm201,] <- assayData(BSData)$G[indA2,][unSortAdd2$ix,][norm201_2,]
  Min.Beads[norm1,] <- assayData(BSData)$no.beads[indA1,][unSortAdd1$ix,][norm1,]
  Min.Beads[norm101,] <- 
    pmin(assayData(BSData)$no.beads[indA1,][unSortAdd1$ix,][norm101,],
         assayData(BSData)$no.beads[indA2,][unSortAdd2$ix,][norm101_2,])
  Min.Beads[norm201,] <- 
    pmin(assayData(BSData)$no.beads[indA1,][unSortAdd1$ix,][norm201,],
         assayData(BSData)$no.beads[indA2,][unSortAdd2$ix,][norm201_2,])
  pol <- cart2pol(A,B,dist=normOpts$dist,pNorm=normOpts$pNorm)
  Var.pooled[norm1,] <- 
    (assayData(BSData)$se.R[indA1,][unSortAdd1$ix,][norm1,]^2 +
     assayData(BSData)$se.G[indA1,][unSortAdd1$ix,][norm1,]^2)/2
  Var.pooled[norm101,] <- 
    (assayData(BSData)$se.R[indA1,][unSortAdd1$ix,][norm101,]^2 +
     assayData(BSData)$se.R[indA2,][unSortAdd2$ix,][norm101_2,]^2)/2
  Var.pooled[norm201,] <- 
    (assayData(BSData)$se.G[indA1,][unSortAdd1$ix,][norm201,]^2 +
     assayData(BSData)$se.G[indA2,][unSortAdd2$ix,][norm201_2,]^2)/2
  SE <- findSeTheta(sqrt(Var.pooled),pol$r,dist=normOpts$dist,pNorm=normOpts$pNorm)
  
  beadMetadata = data.frame(labelDescription=character(ncol(beadInfo)),
    row.names=colnames(beadInfo),stringsAsFactors=FALSE)
  beadMetadata[c('Name','SNP','ILMN.Strand','Address','Address2','Norm.ID'),1] <- 
    c('Marker ID', 'Polymorphism', 'TOP/BOT-category', 'Illumina strand ID',
      'Second strand ID (Infinium I)', 'Sub-bead pool ID')
  featureData <- new("AnnotatedDataFrame", data=beadInfo, varMetadata=beadMetadata,
                     dimLabels=c('featureNames','featureColumns'))
  BSRed <- new('AlleleSetIllumina', intensity=pol$r, theta=pol$th*2/pi,
               SE=SE, min.beads=Min.Beads, phenoData=phenoData(BSData),
               featureData=featureData)
  if (includeAB){
    assayData(BSRed)$A <- A
    assayData(BSRed)$B <- B
  } #end if
  validObject(BSRed)
  rm(pol,A,B,Var.pooled,SE,Min.Beads,featureData,beadInfo,BSData); gc()
  BSRed
}



#Transform cartesian coordinates to polar
cart2pol <- function(x,y,dist='euclidean',pNorm=NULL){
  pol <- list()
  switch(dist,
         euclidean={pol$r <- sqrt(x^2 + y^2)},
         manhattan={pol$r <- x+y},
         minkowski={pol$r <- (abs(x)^pNorm + abs(y)^pNorm)^(1/pNorm)}
         )
  pol$th <- atan2(y,x)
  pol
}



#Make sure 'dataFiles' points to the relevant data. NB! 'beadInfo' MUST contain
#total number of SNPs in the dataFiles, and 'markers' contains row-indexes to
#beadInfo. The
#vectors 'markers' and 'arrays' are indexes to rows and columns, i.e. they do NOT
#depend on the order in which the markers and arrays are given in the files.
createAlleleSetFromFiles <- function(dataFiles,markers,arrays,phenoInfo=NULL,beadInfo=NULL,sep='\t',quote=''){
  if (!all(c('intFile','thFile','seFile') %in% names(dataFiles)))
    stop("The filenames 'intFile', 'thFile', and 'seFile' must all be provided")
  if ('resFile' %in% names(dataFiles)){
    message("Loading featureData (ignoring 'beadInfo'-argument) ...")
    beadInfo <- read.table(dataFiles[['resFile']],header=TRUE,sep=sep,quote=quote,as.is=TRUE)
    fMetadata <- data.frame(labelDescription=character(ncol(beadInfo)), row.names=colnames(beadInfo), stringsAsFactors=FALSE)
    fMetadata[c('Name','SNP','ILMN.Strand','Address','Address2','Norm.ID','Classification',  #try()?
                'Cent.Deviation','Within.SD','HW.Chi2','HW.P','BAF.Locus1','BAF.Locus2',
                'Call.Rate'),1] <- 
                  c('Marker ID','Polymorphism','TOP/BOT-category','Illumina strand ID',
                    'Second strand ID (Infinium I)','Sub-bead pool ID',
                    'Genotype call from automatic clustering',
                    'Largest distance from cluster-centre to its ideal position',
                    'Largest within-cluster spread',
                    'Chi-squared statistic from test of Hardy-Weinberg equilibrium',
                    'Probability of Hardy-Weinberg equilibrium',
                    'Estimated B-allele frequency',
                    'Estimated BAF of second paralogue (if exists)',
                    'Ratio of samples assigned to clusters')
  }else if(!is.null(beadInfo)){
    fMetadata <- data.frame(labelDescription=character(ncol(beadInfo)), row.names=colnames(beadInfo),stringsAsFactors=FALSE)
    fMetadata[c('Name','SNP','ILMN.Strand','Address','Address2','Norm.ID'),1] <- 
      c('Marker ID', 'Polymorphism', 'TOP/BOT-category', 'Illumina strand ID',
        'Second strand ID (Infinium I)', 'Sub-bead pool ID')
  }else{
    stop("No available featureData. Provide 'beadInfo'-table or 'resFile'-filename")
  }
  allMarkers <- unlist(read.table(dataFiles[['intFile']],nrows=1,as.is=TRUE,sep=sep))
  if (missing(markers)){
    markers <- 1:min(nrow(beadInfo),144)
  }else if (is.character(markers)){
    markers <- which(allMarkers %in% markers)
  }else{
    markers <- (1:length(allMarkers))[markers]
  }
  featureData <- new("AnnotatedDataFrame", data=beadInfo[markers,], varMetadata=fMetadata, dimLabels=c('featureNames','featureColumns'))
  
  if ('phFile' %in% names(dataFiles)){
    message("Loading phenoData (ignoring 'phenoInfo'-argument) ...")
    phenoInfo <- read.table(dataFiles[['phFile']],header=TRUE,sep=sep,quote=quote,as.is=TRUE)
    phenoInfo$noiseIntensity <- as.numeric(phenoInfo$noiseIntensity)
  }else if(is.null(phenoInfo)){
    stop("No available phenoData. Provide 'phenoInfo'-table or 'phFile'-filename")
  }
  pMetadata <- data.frame(labelDescription=character(ncol(phenoInfo)), row.names=colnames(phenoInfo), stringsAsFactors=FALSE)
  pMetadata[c('arrayNames','no.beads','chip','row','col','noiseIntensity'),1] <-
    c('Array ID','No. detected beads','Chip ID','Array ID within chip',
      'Strip ID within array','Estimated signal detection-limit')
  if (missing(arrays)){
    #arrays <- 1:min(nrow(phenoInfo),500)
    arrays <- 1:nrow(phenoInfo)
  }else if(is.character(arrays)){
    arrays <- which(rownames(phenoInfo) %in% arrays)
  }else{
    arrays <- (1:nrow(phenoInfo))[arrays]
  }
  phenoData <- new("AnnotatedDataFrame", data=phenoInfo[arrays,], varMetadata=pMetadata, dimLabels=c('sampleNames','sampleColumns'))

  colClasses <- rep('NULL',length(allMarkers)+1)
  colClasses[1] <- 'character'
  colClasses[markers+1] <- 'numeric'
  orderM <- order(sort(markers,index.return=TRUE)$ix)
  message('Loading intensities...')
  intensity <- t(read.table(dataFiles[['intFile']],header=TRUE,sep=sep,quote=quote,colClasses=colClasses,nrows=max(arrays)))
  intensity <- intensity[orderM,arrays]
  message('Loading thetas...')
  theta <- t(read.table(dataFiles[['thFile']],header=TRUE,sep=sep,quote=quote,colClasses=colClasses,nrows=max(arrays)))
  theta <- theta[orderM,arrays]
  message('Loading standard errors...')
  SE <- t(read.table(dataFiles[['seFile']],header=TRUE,sep=sep,quote=quote,colClasses=colClasses,nrows=max(arrays)))
  SE <- SE[orderM,arrays]
  BSRed <- new('AlleleSetIllumina', intensity=intensity, theta=theta, SE=SE, phenoData=phenoData, featureData=featureData)
  rm('intensity','theta','SE','phenoData','featureData'); gc()

  if ('callFile' %in% names(dataFiles)){
    callHeader <- c('',unlist(read.table(dataFiles[['callFile']],nrows=1,as.is=TRUE,sep=sep)))
    callCols <- rep('NULL',length(callHeader))
    callCols[1] <- 'character'
    callCols[callHeader %in% sampleNames(BSRed)] <- 'numeric'
    indFound <- sampleNames(BSRed) %in% callHeader
    orderA <- order(sort(arrays[indFound],index.return=TRUE)$ix)
    call <- matrix(nrow=length(markers),ncol=length(arrays),dimnames=list(featureNames(BSRed),sampleNames(BSRed)))
    message('Loading calls...')
    skip <- min(markers)
    callFrame <- read.table(dataFiles[['callFile']],header=FALSE,sep=sep,quote=quote,nrows=max(markers)+1-skip,colClasses=callCols,skip=skip,row.names=1,col.names=callHeader)
    call[,indFound] <- as.matrix(callFrame[markers+1-skip,orderA])
    assayData(BSRed)$call <- call
    validObject(BSRed)
  }
  if ('genoFile' %in% names(dataFiles)){
    genoHeader <- c('',unlist(read.table(dataFiles[['genoFile']],nrows=1,as.is=TRUE,sep=sep)))
    genoCols <- rep('NULL',length(genoHeader))
    genoCols[1] <- 'character'
    genoCols[genoHeader %in% sampleNames(BSRed)] <- 'character'
    indFound <- sampleNames(BSRed) %in% genoHeader
    orderA <- order(sort(arrays[indFound],index.return=TRUE)$ix)
    genotype <- matrix('',nrow=length(markers),ncol=length(arrays),dimnames=list(featureNames(BSRed),sampleNames(BSRed)))
    message('Loading genotypes...')
    skip <- min(markers)
    genoFrame <- read.table(dataFiles[['genoFile']],header=FALSE,sep=sep,quote=quote,nrows=max(markers)+1-skip,colClasses=genoCols,skip=skip,row.names=1,col.names=genoHeader)
    genotype[,indFound] <- gsub('\"','',as.matrix(genoFrame[markers+1-skip,orderA]))
    assayData(BSRed)$genotype <- genotype
    validObject(BSRed)
  }
  BSRed
}



createMultiSetFromFiles <- function(dataFiles,markers,arrays,phenoInfo=NULL,beadInfo=NULL,sep='\t',quote=''){
  if ('resFile' %in% names(dataFiles)){
    message("Loading featureData (ignoring 'beadInfo'-argument) ...")
    beadInfo <- read.table(dataFiles[['resFile']],header=TRUE,sep=sep,quote=quote,as.is=TRUE)
    fMetadata <- data.frame(labelDescription=character(ncol(beadInfo)), row.names=colnames(beadInfo), stringsAsFactors=FALSE)
    fMetadata[c('Name','SNP','ILMN.Strand','Address','Address2','Norm.ID','Classification',  #try()?
                'Cent.Deviation','Within.SD','HW.Chi2','HW.P','BAF.Locus1','BAF.Locus2',
                'Call.Rate'),1] <- 
                  c('Marker ID','Polymorphism','TOP/BOT-category','Illumina strand ID',
                    'Second strand ID (Infinium I)','Sub-bead pool ID',
                    'Genotype call from automatic clustering',
                    'Largest distance from cluster-centre to its ideal position',
                    'Largest within-cluster spread',
                    'Chi-squared statistic from test of Hardy-Weinberg equilibrium',
                    'Probability of Hardy-Weinberg equilibrium',
                    'Estimated B-allele frequency',
                    'Estimated BAF of second paralogue (if exists)',
                    'Ratio of samples assigned to clusters')
  }else if(!is.null(beadInfo)){
    fMetadata <- data.frame(labelDescription=character(ncol(beadInfo)), row.names=colnames(beadInfo),stringsAsFactors=FALSE)
    fMetadata[c('Name','SNP','ILMN.Strand','Address','Address2','Norm.ID'),1] <- 
      c('Marker ID', 'Polymorphism', 'TOP/BOT-category', 'Illumina strand ID',
        'Second strand ID (Infinium I)', 'Sub-bead pool ID')
  }else{
    stop("No available featureData. Provide 'beadInfo'-table or 'resFile'-filename")
  }
  allMarkers <- rownames(beadInfo)
  if (missing(markers)){
    markers <- 1:min(nrow(beadInfo),128)
  }else if (is.character(markers)){
    markers <- which(allMarkers %in% markers)
  }else{
    markers <- (1:length(allMarkers))[markers]
  }
  featureData <- new("AnnotatedDataFrame", data=beadInfo[markers,], varMetadata=fMetadata, dimLabels=c('featureNames','featureColumns'))
  
  if ('phFile' %in% names(dataFiles)){
    message("Loading phenoData (ignoring 'phenoInfo'-argument) ...")
    phenoInfo <- read.table(dataFiles[['phFile']],header=TRUE,sep=sep,quote=quote,as.is=TRUE)
  }else if(is.null(phenoInfo)){
    stop("No available phenoData. Provide 'phenoInfo'-table or 'phFile'-filename")
  }
  pMetadata <- data.frame(labelDescription=character(ncol(phenoInfo)), row.names=colnames(phenoInfo), stringsAsFactors=FALSE)
  pMetadata[c('arrayNames','no.beads','chip','row','col','noiseIntensity'),1] <-
    c('Array ID','No. detected beads','Chip ID','Array ID within chip',
      'Strip ID within array','Estimated signal detection-limit')
  if (missing(arrays)){
    #arrays <- 1:min(nrow(phenoInfo),500)
    arrays <- 1:nrow(phenoInfo)
  }else if(is.character(arrays)){
    arrays <- which(rownames(phenoInfo) %in% arrays)
  }else{
    arrays <- (1:nrow(phenoInfo))[arrays]
  }
  phenoData <- new("AnnotatedDataFrame", data=phenoInfo[arrays,], varMetadata=pMetadata, dimLabels=c('sampleNames','sampleColumns'))

  colClasses <- rep('NULL',length(allMarkers)+1)
  colClasses[1] <- 'character'
  colClasses[markers+1] <- 'numeric'
  orderM <- order(sort(markers,index.return=TRUE)$ix)
  if ('intFile' %in% names(dataFiles)){
    message('Loading intensities...')
    intensity <- t(read.table(dataFiles[['intFile']],header=TRUE,sep=sep,quote=quote,colClasses=colClasses,nrows=max(arrays)))
    intensity <- intensity[orderM,arrays]
  }else{
    intensity <- NULL
  }
  if ('thFile' %in% names(dataFiles)){
    message('Loading thetas...')
    theta <- t(read.table(dataFiles[['thFile']],header=TRUE,sep=sep,quote=quote,colClasses=colClasses,nrows=max(arrays)))
    theta <- theta[orderM,arrays]
  }else{
    theta <- NULL
  }
  if ('seFile' %in% names(dataFiles)){
    message('Loading standard errors...')
    SE <- t(read.table(dataFiles[['seFile']],header=TRUE,sep=sep,quote=quote,colClasses=colClasses,nrows=max(arrays)))
    SE <- SE[orderM,arrays]
  }else{
    SE <- NULL
  }
  BSRed <- new('MultiSet', intensity=intensity, theta=theta, SE=SE, phenoData=phenoData, featureData=featureData, storage.mode='list')
  sampleNames <- rownames(phenoData@data)
  featureNames <- rownames(featureData@data)
  rm('intensity','theta','SE','phenoData','featureData'); gc()
  
  if ('callFile' %in% names(dataFiles)){
    callHeader <- c('',unlist(read.table(dataFiles[['callFile']],nrows=1,as.is=TRUE,sep=sep)))
    callCols <- rep('NULL',length(callHeader))
    callCols[1] <- 'character'
    callCols[callHeader %in% sampleNames] <- 'numeric'
    indFound <- sampleNames %in% callHeader
    orderA <- order(sort(arrays[indFound],index.return=TRUE)$ix)
    call <- matrix(nrow=length(markers),ncol=length(arrays),dimnames=list(featureNames,sampleNames))
    message('Loading calls...')
    skip <- min(markers)
    callFrame <- read.table(dataFiles[['callFile']],header=FALSE,sep=sep,quote=quote,nrows=max(markers)+1-skip,colClasses=callCols,skip=skip,row.names=1,col.names=callHeader)
    call[,indFound] <- as.matrix(callFrame[markers+1-skip,orderA])
    assayData(BSRed)$call <- call
    validObject(BSRed)
  }
  if ('genoFile' %in% names(dataFiles)){
    genoHeader <- c('',unlist(read.table(dataFiles[['genoFile']],nrows=1,as.is=TRUE,sep=sep)))
    genoCols <- rep('NULL',length(genoHeader))
    genoCols[1] <- 'character'
    genoCols[genoHeader %in% sampleNames] <- 'character'
    indFound <- sampleNames %in% genoHeader
    orderA <- order(sort(arrays[indFound],index.return=TRUE)$ix)
    genotype <- matrix('',nrow=length(markers),ncol=length(arrays),dimnames=list(featureNames,sampleNames))
    message('Loading genotypes...')
    skip <- min(markers)
    genoFrame <- read.table(dataFiles[['genoFile']],header=FALSE,sep=sep,quote=quote,nrows=max(markers)+1-skip,colClasses=genoCols,skip=skip,row.names=1,col.names=genoHeader)
    genotype[,indFound] <- gsub('\"','',as.matrix(genoFrame[markers+1-skip,orderA]))
    assayData(BSRed)$genotype <- genotype
    validObject(BSRed)
  }
  BSRed
}



callGenotypes <- function(BSRed,gO=setGenoOptions(largeSample=ncol(BSRed)>250)){
  #Quality control, SNP- and array-filtering
  if ('call' %in% names(assayData(BSRed))){
    OK = FALSE
    while (!OK){
      ovrwrt <- readline(prompt='Warning! Call-matrix already exists. Overwrite? (yes/no) ')
      OK <- ovrwrt %in% c('yes','no')
    }
    if (ovrwrt %in% 'yes'){
      message("Matrix 'assayData(BSRed)$call' will be replaced!")
      pData(BSRed)$passRatio <- NULL
      fData(BSRed)$Ped.Errors <- NULL
      assayData(BSRed)[grep('ped.check',names(assayData(BSRed)))] <- NULL
    }else{
      stop("Exits, previous calls are retained")
    }
  }
  descNames <- c('Classification','Cent.Deviation','Within.SD','HW.Chi2','HW.P','BAF.Locus1',
                 'BAF.Locus2','Call.Rate')
  if (any(descNames %in% colnames(fData(BSRed)))){
    OK = FALSE
    while (!OK){
      ovrwrt <- readline(prompt='Warning! Overlapping column-names in featureData. Overwrite? (yes/no) ')
      OK <- ovrwrt %in% c('yes','no')
    }
    if (ovrwrt %in% 'yes'){
      message("Overlapping descriptors will be replaced!")
      featureData(BSRed) <- featureData(BSRed)[,-which(colnames(fData(BSRed)) %in% descNames)]
    }else{
      stop("Exits, previous descriptors are retained")
    }
  }
  iSnps <- apply(assayData(BSRed)$theta,1,function(x,cc) sum(!is.na(x))>cc,cc=gO$arrayPerSnpLim*ncol(BSRed))
  message(sum(!iSnps),' SNP(s) with more than ',100-gO$arrayPerSnpLim*100,'% missings disregarded')
  iArrays <- apply(assayData(BSRed)$theta,2,function(x,cc) sum(!is.na(x))>cc,cc=gO$snpPerArrayLim*nrow(BSRed))
  message(sum(!iArrays),' array(s) with more than ',100-gO$snpPerArrayLim*100,'% missings disregarded')
  mNoise <- mean(pData(BSRed)$noiseIntensity[iArrays],na.rm=TRUE)
  if (sum(iSnps)==1){
    R <- t(as.matrix(assayData(BSRed)$intensity[iSnps,iArrays]))
    Theta <- t(as.matrix(assayData(BSRed)$theta[iSnps,iArrays]))
    SE <- t(as.matrix(assayData(BSRed)$SE[iSnps,iArrays]))
    dimnames(R)[1] <- dimnames(Theta)[1] <- dimnames(SE)[1] <- names(iSnps)
  }else{
    R <- assayData(BSRed)$intensity[iSnps,iArrays]
    Theta <- assayData(BSRed)$theta[iSnps,iArrays]
    SE <- assayData(BSRed)$SE[iSnps,iArrays]
  }
  iFilter <- diff(apply(Theta,1,range,na.rm=TRUE))>gO$filterLim
  nKeep <- sum(iFilter)
  message(sum(!iFilter),' monomorphs with range(theta)<',gO$filterLim,' are detected')
  
  polyCent <- generatePolyCenters(ploidy=gO$ploidy)
  assayData(BSRed)$call <- matrix(nrow=nrow(BSRed),ncol=ncol(BSRed),dimnames=list(featureNames(BSRed),sampleNames(BSRed)))
  classVec <- character(nrow(BSRed))
  classVec[!iSnps] <- 'NA'
  classVec[iSnps][!iFilter] <- 'MONO-filt'
  fData <- data.frame(Classification=classVec,Cent.Deviation=rep(NA,nrow(BSRed)),Within.SD=rep(NA,nrow(BSRed)),HW.Chi2=rep(NA,nrow(BSRed)),HW.P=rep(NA,nrow(BSRed)),BAF.Locus1=rep(NA,nrow(BSRed)),BAF.Locus2=rep(NA,nrow(BSRed)),Call.Rate=rep(NA,nrow(BSRed)),stringsAsFactors=FALSE,row.names=featureNames(BSRed))
  fDesc <-  c('Genotype call from automatic clustering',
              'Largest distance from cluster-centre to its ideal position',
              'Largest within-cluster spread',
              'Chi-squared statistic from test of Hardy-Weinberg equilibrium',
              'Probability of Hardy-Weinberg equilibrium',
              'Estimated B-allele frequency',
              'Estimated BAF of second paralogue (if exists)',
              'Ratio of samples assigned to clusters')
  fMetadata <- data.frame(labelDescription=fDesc, row.names=colnames(fData), stringsAsFactors=FALSE)
  count <- 0
  for (i in which(iFilter)){ 
    count <- count + 1
    if (count %in% seq(min(100,nKeep),nKeep,100))
      message('Called ',count,' of ',nKeep,' genotypes')
    iMiss <- is.na(R[i,])|is.na(SE[i,])|is.na(Theta[i,])
    iDet <- R[i,]>mNoise & !iMiss & Theta[i,]>-.5 & Theta[i,]<1.5    #Filter bad samples
    #if (sum(iDet)<=sum(iArrays)*gO$detectLim){
    if (sum(iDet)<=sum(!iMiss)*gO$detectLim){
      fData$Classification[iSnps][[i]] <- 'FAIL'
    }else{
      sclR <- median(R[i,iDet],na.rm=TRUE)*2*gO$rPenalty                   #Factor by which R is divided before clustering
      X <- matrix(cbind(Theta[i,],R[i,]/sclR,SE[i,])[iDet,],ncol=3)
      if (is.null(gO$probsIndSE)){
        indSE <- rep(TRUE,sum(iDet))
      }else{
        indSE <- X[,3]<quantile(X[,3],probs=gO$probsIndSE,na.rm=TRUE)      #Uncertain points removed from clustering
      }
      #qRange <- getQRange(X[indSE,1],minScl=3,probs=gO$probsQRange)     #Initial cluster-centers
      #sConf <- distLikelyConfiguration(X[indSE,1:2],polyCent,qRange=qRange)  #Ordered list of likely classifications
      #maxConf <- sum(sConf$x<gO$mseLim,na.rm=TRUE)
      sConf <- getCenters(X[indSE,1],gO=gO,polyCent=polyCent)
      maxConf <- length(sConf$ix)
      OK <- FALSE
      iConf <- 1
      while (!OK){
        testCallRate <- FALSE
        goodSNP <- FALSE
        nCl <- polyCent$size[sConf$ix[iConf]]
        cntIdeal <- polyCent$centers[[sConf$ix[iConf]]]
        centers <- sConf$centers[[iConf]]
        clObj <- findPolyploidClusters(X,indSE=indSE,centers=centers,wss.update=FALSE)
        tmpRes <- rep(NA,7)
        if (any(clObj$size>0)){
          iCl <- which(clObj$size>0)
          tmpRes[1] <- max(abs(cntIdeal[iCl]-clObj$centers[iCl,1]))   #Max cluster-centre deviation (Theta)
          #tmpRes[1] <- mean(abs(cntIdeal[iCl]-clObj$centers[iCl,1]))   #Max cluster-centre deviation (Theta)
          if (tmpRes[1]<gO$devCentLim & sum(clObj$size<1)<2){
            wSpread <- tapply(X[,1],clObj$cluster,sd)   #NB! MAD?
            wSpread[clObj$size[iCl]<2] <- 0   #use gO$minClLim instead?
            tmpRes[2] <- max(wSpread)  #Max within-cluster spread (Theta)
            if (tmpRes[2]<gO$wSpreadLim){
              if (nCl==1){
                tmpRes[5:6] <- as.numeric(testHardyWeinberg(sizes=clObj$size,bestConf=sConf$ix[iConf],polyCent=polyCent,afList=gO$afList)[4:5])
                testCallRate <- TRUE
              }else{
                tmpRes[3:6] <- as.numeric(testHardyWeinberg(sizes=clObj$size,bestConf=sConf$ix[iConf],polyCent=polyCent,afList=gO$afList)[c(1,3:5)])
                minB <- (clObj$centers[iCl,1] - gO$nSdOverlap*wSpread)[-1]
                maxB <- (clObj$centers[iCl,1] + gO$nSdOverlap*wSpread)[-length(iCl)]
                ovrlp <- any(round(minB*100)<round(maxB*100))   #Overlapping clusters?
                                        #goodSNP <- tmpRes[1]+2*tmpRes[2]<gO$goodSnpLim  & sConf$ix[iConf]==which(polyCent$classification=='SNP') & sum(clObj$size>=min(gO$minClLim,ncol(BSRed)/100))==3
                                        #if (tmpRes[4]>gO$hwAlpha | goodSNP){
                if (tmpRes[4]>gO$hwAlpha & !ovrlp){
                  testCallRate <- TRUE
                }else{
                  iConf <- iConf + 1
                }
              }
              if (testCallRate){
                nClRed <- length(iCl)
                clCtrs <- matrix(clObj$centers[iCl,],nrow=nClRed)
                indOut <- matrix(FALSE,nrow=nrow(X),ncol=nCl)
                indKeep <- rep(FALSE,nrow(X))
                indOverlap <- rep(FALSE,nrow(X))
                Sratio <- rep(1,nCl)
                for (j in iCl){
                  indj <- clObj$cluster==j
                  nj <- sum(indj)
                  if (nj>=gO$minClLim){    #Hotellings T^2-ellipse superimposed on clusters
                    ctrX <- apply(X[indj,1:2],2,median,na.rm=TRUE)
                    Xj <- X[indj,1:2]-matrix(ctrX,nrow=nj,ncol=2,byrow=TRUE)
                    Sj <- t(Xj)%*%Xj/(nj-1)
                    Sratio[j] <- Sj[1,1]/(Sj[2,2]*gO$rPenalty)  #Cluster 'angle', given comparable axes
                    Xall <- X[,1:2]-matrix(ctrX,nrow=nrow(X),ncol=2,byrow=TRUE)
                    T2j <- diag(Xall%*%solve(Sj)%*%t(Xall))
                    f.alpha <- T2j*(nj-2)*nj/(2*(nj-1)*(nj+1))
                    pVal <- 1 - pf(f.alpha,2,nj-2)
                  }else{
                    pVal <- rep(0,nrow(X))
                    pVal[indj] <- 1
                  }
                  indOut[,j] <- pVal < gO$clAlpha
                  indKeep[indj] <- !indOut[indj,j]
                  indOverlap[!indOut[,j]][clObj$cluster[!indOut[,j]]!=j] <- TRUE
                }
                plot <- FALSE
                if (plot){
                  plot.new()
                  plot.window(c(-.2,1.2),c(0,max(max(X[,2]),1/gO$rPenalty)))
                  axis(1); axis(2)
                  points(X[,1],X[,2],pch=16)
                  points(X[!indKeep,1],X[!indKeep,2],pch=16,col='red')   #outside ellipses
                  points(X[indOverlap,1],X[indOverlap,2],pch=16,col='orange')  #overlapping ellipses
                  points(clCtrs[,1],clCtrs[,2],col='green',pch=16)    #cluster-centres
                }
                indKeep <- indKeep &! indOverlap
                #tmpRes[7] <- sum(indKeep)/sum(iArrays)
                tmpRes[7] <- sum(indKeep)/sum(!iMiss)
                                        #psRat <- prod(Sratio[clObj$size>ncol(BSRed)/10])
                psRat <- sum(Sratio*clObj$size)/(ncol(BSRed)*nClRed)
                                        #goodSNP <- goodSNP & sum(!indOverlap)==sum(iDet)
                
                if ((tmpRes[7]>gO$detectLim & psRat<=gO$rotationLim)){# | goodSNP){   #Success!
                  if (sum(iSnps)==1){
                    assayData(BSRed)$call[iSnps,iArrays][iDet][indKeep] <- cntIdeal[clObj$cluster][indKeep]
                  }else{
                    assayData(BSRed)$call[iSnps,iArrays][i,iDet][indKeep] <- cntIdeal[clObj$cluster][indKeep]
                  }
                  fData$Classification[iSnps][[i]] <- polyCent$classification[[sConf$ix[iConf]]]
                  fData$Cent.Deviation[iSnps][[i]] <- tmpRes[1]
                  fData$Within.SD[iSnps][i] <- tmpRes[2]
                  fData$HW.Chi2[iSnps][i] <- tmpRes[3]
                  fData$HW.P[iSnps][i] <- tmpRes[4]
                  fData$BAF.Locus1[iSnps][i] <- 1-tmpRes[5]  #NB! B-af, not A-af
                  fData$BAF.Locus2[iSnps][i] <- 1-tmpRes[6]  #NB!
                  fData$Call.Rate[iSnps][i] <- tmpRes[7]
                  #fData$Max.Rotation[iSnps][i] <- max(Sratio)
                  OK <- TRUE
                }else{
                  iConf <- iConf + 1
                } #end if success
              } #end if callrate
            }else{
              iConf <- iConf + 1
            } #end if wSpread
          }else{
            iConf <- iConf + 1
          } #end if devCent
          if (iConf>maxConf){
            fData$Classification[iSnps][[i]] <- 'FAIL'
            OK <- TRUE
          }
        }else{
          fData$Classification[iSnps][[i]] <- 'FAIL'
          OK <- TRUE
        } #end if any(clObj$size)
      } #end while !OK
    } #end if gO$detectLim
  } #end for
  featureData <- new('AnnotatedDataFrame', data=fData, varMetadata=fMetadata, dimLabels=c('featureNames','featureColumns'))
  featureData(BSRed) <- combine(featureData(BSRed), featureData)
  validObject(BSRed)
  BSRed
}


#Tests automatic genocalling for a single marker, plots different outcomes. Outputs all results from single marker
callGenotypes.verboseTest <- function(BSRed,singleMarker=1,gO=setGenoOptions(largeSample=ncol(BSRed)>250)){
  message(sprintf("Verbosely analyzing marker '%s'",featureNames(BSRed)[singleMarker]))
  iSnps <- apply(assayData(BSRed)$theta,1,function(x,cc) sum(!is.na(x))>cc,cc=gO$arrayPerSnpLim*ncol(BSRed))
  if (!iSnps[singleMarker]){
    message('Failure as the ratio of missing arrays exceeds (1 - "arrayPerSnpLim")')
  }
  iSnps <- singleMarker
  iArrays <- apply(assayData(BSRed)$theta,2,function(x,cc) sum(!is.na(x))>cc,cc=gO$snpPerArrayLim*nrow(BSRed))
  message(sum(!iArrays),' array(s) with more than ',100-gO$snpPerArrayLim*100,'% missings found')
  mNoise <- mean(pData(BSRed)$noiseIntensity[iArrays],na.rm=TRUE)
  R <- as.matrix(assayData(BSRed)$intensity[singleMarker,iArrays])
  Theta <- as.matrix(assayData(BSRed)$theta[singleMarker,iArrays])
  SE <- as.matrix(assayData(BSRed)$SE[singleMarker,iArrays])
  dimnames(R)[2] <- dimnames(Theta)[2] <- dimnames(SE)[2] <- featureNames(BSRed)[singleMarker]
  iFilter <- diff(range(Theta,na.rm=TRUE))>gO$filterLim
  if (!iFilter){
    message('The selected marker would be filtered away as a monomorph based on range(Theta)')
    return()
  }
  BSRed <- BSRed[singleMarker,]
  polyCent <- generatePolyCenters(ploidy=gO$ploidy)
  #nTest <- length(polyCent$size)
  #message(sprintf('%d allowed cluster-combinations',nTest))
  #call <- matrix(nrow=nTest,ncol=ncol(BSRed),dimnames=list(polyCent$classification,sampleNames(BSRed)))
  #fData <- data.frame(Cent.Deviation=rep(NA,nTest),Within.SD=rep(NA,nTest),HW.Chi2=rep(NA,nTest),HW.P=rep(NA,nTest),BAF.Locus1=rep(NA,nTest),BAF.Locus2=rep(NA,nTest),Call.Rate=rep(NA,nTest),Overlap=rep(NA,nTest),Rotation=rep(NA,nTest),stringsAsFactors=FALSE,row.names=polyCent$classification)
  #fDesc <-  c('Largest distance from cluster-centre to its ideal position',
  #            'Largest within-cluster spread',
  #            'Chi-squared statistic from test of Hardy-Weinberg equilibrium',
  #            'Probability of Hardy-Weinberg equilibrium',
  #            'Estimated B-allele frequency',
  #            'Estimated BAF of second paralogue (if exists)',
  #            'Ratio of samples assigned to clusters',
  #            'TRUE if there is an initial cluster overlap (along \'Theta\')',
  #            'Large value (>1) means increasingly slanting clusters')
  #fMetadata <- data.frame(labelDescription=fDesc, row.names=colnames(fData), stringsAsFactors=FALSE)

  iMiss <- is.na(R)|is.na(SE)|is.na(Theta)
  iDet <- R>mNoise & !iMiss & Theta>-.5 & Theta<1.5    #Filter bad arrays
  message(sprintf('%d low-quality or missing array(s) filtered away for this marker',sum(!iDet)))
  #if (sum(iDet)<=sum(iArrays)*gO$detectLim){
  if (sum(iDet)<=sum(!iMiss)*gO$detectLim){
    str <- 'Failure as rate of high quality to non-missing arrays (%.2f) smaller than "detectLim"'
    message(sprintf(str,sum(iDet)/sum(!iMiss)))
  }else{
    sclR <- median(R[iDet],na.rm=TRUE)*2*gO$rPenalty                   #Factor by which R is divided before clustering
    X <- matrix(cbind(Theta,R/sclR,SE)[iDet,],ncol=3)
    if (is.null(gO$probsIndSE)){
      indSE <- rep(TRUE,sum(iDet))
    }else{
      indSE <- X[,3]<quantile(X[,3],probs=gO$probsIndSE,na.rm=TRUE)      #Uncertain points removed from clustering
    }
    #qRange <- getQRange(X[indSE,1],minScl=3,probs=gO$probsQRange)     #Initial cluster-centers
    #sConf <- distLikelyConfiguration(X[indSE,1:2],polyCent,qRange=qRange)  #Ordered list of likely classifications
    sConf <- getCenters(X[indSE,1],gO=gO,polyCent=polyCent)
    print(sConf)
    nTest <- length(sConf$ix)
    message(sprintf('%d cluster-combinations tested',nTest))
    rNames <- unlist(polyCent$classification[sConf$ix])
    if (nTest>1){
      for (i in 2:length(rNames)){
        sumRep <- sum(rNames[1:(i-1)]%in%rNames[i])
        if (sumRep>0){
          rNames[i] <- paste(rNames[i],sumRep,sep='.')
        }
      }
    }
    call <- matrix(nrow=nTest,ncol=ncol(BSRed),dimnames=list(rNames,sampleNames(BSRed)))
    fData <- data.frame(Cent.Deviation=rep(NA,nTest),Within.SD=rep(NA,nTest),HW.Chi2=rep(NA,nTest),HW.P=rep(NA,nTest),BAF.Locus1=rep(NA,nTest),BAF.Locus2=rep(NA,nTest),Call.Rate=rep(NA,nTest),Overlap=rep(NA,nTest),Rotation=rep(NA,nTest),stringsAsFactors=FALSE,row.names=rNames)
    fDesc <-  c('Largest distance from cluster-centre to its ideal position',
                'Largest within-cluster spread',
                'Chi-squared statistic from test of Hardy-Weinberg equilibrium',
                'Probability of Hardy-Weinberg equilibrium',
                'Estimated B-allele frequency',
                'Estimated BAF of second paralogue (if exists)',
                'Ratio of samples assigned to clusters',
                'TRUE if there is an initial cluster overlap (along \'Theta\')',
                'Large value (>1) means increasingly slanting clusters')
    fMetadata <- data.frame(labelDescription=fDesc, row.names=colnames(fData), stringsAsFactors=FALSE)

    pp <- par()
    nC <- ceiling(sqrt(nTest))
    nR <- ceiling(nTest/nC)
    par(mfrow=c(nR,nC),mai=c(1,1,1,.5)/sqrt(nC))
    
    for (iConf in 1:nTest){
      current.i <- sConf$ix[iConf]
      nCl <- polyCent$size[current.i]
      cntIdeal <- polyCent$centers[[current.i]]
      centers <- sConf$centers[[iConf]]
      clObj <- findPolyploidClusters(X,indSE=indSE,centers=centers,wss.update=FALSE)
      if (!any(clObj$size>0)){
        message(sprintf("No clusters assigned for '%s'",polyCent$classification[current.i]))
      }else{
        iCl <- which(clObj$size>0)
        #fData[current.i,1] <- max(abs(cntIdeal[iCl]-clObj$centers[iCl,1]))   #Max cluster-centre deviation (Theta)
        fData[iConf,1] <- max(abs(cntIdeal[iCl]-clObj$centers[iCl,1]))   #Max cluster-centre deviation (Theta)
        wSpread <- tapply(X[,1],clObj$cluster,sd)   #NB! MAD?
        wSpread[clObj$size[iCl]<2] <- 0   #use gO$minClLim instead?
        #fData[current.i,2] <- max(wSpread)  #Max within-cluster spread (Theta)
        fData[iConf,2] <- max(wSpread)  #Max within-cluster spread (Theta)
        if (nCl==1){
          #fData[current.i,5:6] <- as.numeric(testHardyWeinberg(sizes=clObj$size,bestConf=current.i,polyCent=polyCent,afList=gO$afList)[4:5])
          fData[iConf,5:6] <- as.numeric(testHardyWeinberg(sizes=clObj$size,bestConf=current.i,polyCent=polyCent,afList=gO$afList)[4:5])
        }else{
          #fData[current.i,3:6] <- as.numeric(testHardyWeinberg(sizes=clObj$size,bestConf=current.i,polyCent=polyCent,afList=gO$afList)[c(1,3:5)])
          fData[iConf,3:6] <- as.numeric(testHardyWeinberg(sizes=clObj$size,bestConf=current.i,polyCent=polyCent,afList=gO$afList)[c(1,3:5)])
          minB <- (clObj$centers[iCl,1] - gO$nSdOverlap*wSpread)[-1]
          maxB <- (clObj$centers[iCl,1] + gO$nSdOverlap*wSpread)[-length(iCl)]
          #fData[current.i,8] <- any(round(minB*100)<round(maxB*100))
          fData[iConf,8] <- any(round(minB*100)<round(maxB*100))
        }
        nClRed <- length(iCl)
        clCtrs <- matrix(clObj$centers[iCl,],nrow=nClRed)
        indOut <- matrix(FALSE,nrow=nrow(X),ncol=nCl)
        indKeep <- rep(FALSE,nrow(X))
        indOverlap <- rep(FALSE,nrow(X))
        Sratio <- rep(1,nCl)
        for (j in iCl){
          indj <- clObj$cluster==j
          nj <- sum(indj)
          if (nj>=gO$minClLim){    #Hotellings T^2-ellipse superimposed on clusters
            ctrX <- apply(X[indj,1:2],2,median,na.rm=TRUE)
            Xj <- X[indj,1:2]-matrix(ctrX,nrow=nj,ncol=2,byrow=TRUE)
            Sj <- t(Xj)%*%Xj/(nj-1)
            Sratio[j] <- Sj[1,1]/(Sj[2,2]*gO$rPenalty)  #Cluster 'angle', given comparable axes
            Xall <- X[,1:2]-matrix(ctrX,nrow=nrow(X),ncol=2,byrow=TRUE)
            T2j <- diag(Xall%*%solve(Sj)%*%t(Xall))
            f.alpha <- T2j*(nj-2)*nj/(2*(nj-1)*(nj+1))
            pVal <- 1 - pf(f.alpha,2,nj-2)
          }else{
            pVal <- rep(0,nrow(X))
            pVal[indj] <- 1
          }
          indOut[,j] <- pVal < gO$clAlpha
          indKeep[indj] <- !indOut[indj,j]
          indOverlap[!indOut[,j]][clObj$cluster[!indOut[,j]]!=j] <- TRUE
        }
        
        plot.new()
        plot.window(c(-.2,1.2),c(0,max(max(X[,2]),1/gO$rPenalty)))
        axis(1); axis(2)
        title(main=polyCent$classification[current.i],xlab='Theta',ylab='R')
        points(X[,1],X[,2],pch=16)
        points(X[!indKeep,1],X[!indKeep,2],pch=16,col='red')   #outside ellipses
        points(X[indOverlap,1],X[indOverlap,2],pch=16,col='orange')  #overlapping ellipses
        points(clCtrs[,1],clCtrs[,2],col='green',pch=16)    #cluster-centres
        
        indKeep <- indKeep &! indOverlap
        #fData[current.i,7] <- sum(indKeep)/sum(iArrays)
        #fData[current.i,9] <- sum(Sratio*clObj$size)/(ncol(BSRed)*nClRed)
        #call[current.i,iArrays][iDet][indKeep] <- cntIdeal[clObj$cluster][indKeep]
        #fData[iConf,7] <- sum(indKeep)/sum(iArrays)
        fData[iConf,7] <- sum(indKeep)/sum(!iMiss)
        fData[iConf,9] <- sum(Sratio*clObj$size)/(ncol(BSRed)*nClRed)
        call[iConf,iArrays][iDet][indKeep] <- cntIdeal[clObj$cluster][indKeep]
      } #end if (!any(clObj$size))
    } #for
    par(pp[c('mfcol','mai')])
    verboseResults <- list()
    verboseResults$call <- call
    verboseResults$fData <- fData
    verboseResults$fData[,5:6] <- 1-verboseResults$fData[,5:6]  #NB! B-af, not A-af
    verboseResults$fMetadata <- fMetadata
    test <- data.frame(Cent.Deviation=fData[,1]<gO$devCentLim,
                       WithinSD=fData[,2]<gO$wSpreadLim,
                       HW.P=fData[,4]>gO$hwAlpha,Call.Rate=fData[,7]>gO$detectLim,
                       Overlap=!fData[,8],Rotation=fData[,9]<=gO$rotationLim,
                       row.names=rNames)
    verboseResults$test <- test
    verboseResults$test[test==TRUE] <- 'PASS'
    verboseResults$test[test==FALSE] <- 'FAIL'
    verboseResults$test$sumFAIL <- apply(!test,1,sum,na.rm=TRUE)
    #verboseResults$test <- verboseResults$test[sConf$ix,]
    indPass <- verboseResults$test$sumFAIL==0
    if (any(indPass)){
      clss <- rownames(verboseResults$test)[which(indPass)[1]]
    }else{
      clss <- 'FAIL'
    }
    message(sprintf("Marker would be classified as '%s'",clss))
    verboseResults
  } #end if low qual.
}



#Manual genotype calling
callGenotypes.interactive <- function(BSRed,gO=setGenoOptions(largeSample=ncol(BSRed)>250)){
  #require('rggobi')
  mNoise <- mean(pData(BSRed)$noiseIntensity, na.rm=TRUE)
  R <- assayData(BSRed)$intensity
  Theta <- assayData(BSRed)$theta
  SE <- assayData(BSRed)$SE
  polyCent <- generatePolyCenters(ploidy=gO$ploidy)
  if (is.null(assayData(BSRed)$call)){
    assayData(BSRed)$call <-  matrix(nrow=nrow(BSRed), ncol=ncol(BSRed), dimnames=list(featureNames(BSRed),sampleNames(BSRed)))
    isCalled <- FALSE
  }else{
    isCalled <- TRUE
  }
  if (is.null(assayData(BSRed)$ped.check.parents) | !isCalled)
    assayData(BSRed)$ped.check.parents <- matrix(nrow=nrow(BSRed), ncol=ncol(BSRed), dimnames=list(featureNames(BSRed),sampleNames(BSRed)))
  rNames <- c("Classification", "Cent.Deviation", "Within.SD", "HW.Chi2", "HW.P", "BAF.Locus1", "BAF.Locus2", "Call.Rate", "Ped.Errors")   
  if (all(rNames[1:8] %in% colnames(fData(BSRed))) & isCalled){
    prevCalls <- fData(BSRed)$Classification
    if (!rNames[9] %in% colnames(fData(BSRed))){
      fData <- data.frame(Ped.Errors=rep(NA,nrow(BSRed)), stringsAsFactors=FALSE, row.names=featureNames(BSRed))
      fMetadata <- data.frame(labelDescription='Number of pedigree-errors, summed over offspring', row.names=colnames(fData), stringsAsFactors=FALSE)
      featureData <- new('AnnotatedDataFrame', data=fData, varMetadata=fMetadata, dimLabels=c('featureNames','featureColumns'))
      featureData(BSRed) <- combine(featureData(BSRed), featureData)
    }
    if ('Manual.Calls.R' %in% colnames(fData(BSRed))){
      indMiss <- fData(BSRed)$Manual.Calls.R %in% ''
      if (any(!indMiss)){
        message('Column \'Manual.Calls.R\' will be overwritten!')
        prevCalls[!indMiss] <- fData(BSRed)$Manual.Calls.R[!indMiss]
      }
    }else{
      fData <- data.frame(Manual.Calls.R=character(nrow(BSRed)), stringsAsFactors=FALSE, row.names=featureNames(BSRed))
      fMetadata <- data.frame(labelDescription='Genotype call from interactive clustering', row.names=colnames(fData), stringsAsFactors=FALSE)
      featureData <- new('AnnotatedDataFrame', data=fData, varMetadata=fMetadata, dimLabels=c('featureNames','featureColumns'))
      featureData(BSRed) <- combine(featureData(BSRed), featureData)
    }
  }else{
    iKeep <- !colnames(fData(BSRed)) %in% c(rNames,'Manual.Calls.R')
    fData <- data.frame(Cent.Deviation=rep(NA,nrow(BSRed)), Within.SD=rep(NA,nrow(BSRed)), HW.Chi2=rep(NA,nrow(BSRed)), HW.P=rep(NA,nrow(BSRed)), BAF.Locus1=rep(NA,nrow(BSRed)), BAF.Locus2=rep(NA,nrow(BSRed)), Call.Rate=rep(NA,nrow(BSRed)), Ped.Errors=rep(NA,nrow(BSRed)), Manual.Calls.R=character(nrow(BSRed)), stringsAsFactors=FALSE, row.names=featureNames(BSRed))
    fDesc <-  c('Largest distance from cluster-centre to its ideal position',
                'Largest within-cluster spread',
                'Chi-squared statistic from test of Hardy-Weinberg equilibrium',
                'Probability of Hardy-Weinberg equilibrium',
                'Estimated B-allele frequency',
                'Estimated BAF of second paralogue (if exists)',
                'Ratio of samples assigned to clusters',
                'Number of pedigree-errors, summed over offspring',
                'Genotype call from interactive clustering')
    fMetadata <- data.frame(labelDescription=fDesc, row.names=colnames(fData), stringsAsFactors=FALSE)
    featureData <- new('AnnotatedDataFrame', data=fData, varMetadata=fMetadata, dimLabels=c('featureNames','featureColumns'))
    featureData(BSRed) <- combine(featureData(BSRed)[,iKeep], featureData)
    prevCalls <- 'FAIL'
  }
  str0 <- 'Calling marker no. %d of %d...'
  for (i in 1:nrow(BSRed)){
    message(sprintf(str0,i,nrow(BSRed)))
    if (isCalled & !prevCalls[i]%in%'FAIL'){
      marker <- data.frame(Theta=Theta[i,], R=R[i,], Call=assayData(BSRed)$call[i,], PedigreeID=pData(BSRed)$PedigreeID, stringsAsFactors=FALSE)
      gg <- ggobi(marker[,c('Theta','R')],name='marker')
      glc <- marker$Call*4 + 1
      glc[is.na(glc)] <- 8
      marker$PedCheck <- validateSingleCall(marker,prevCalls[i])
      iErr <- marker$PedCheck>0
      glc[iErr] <- 5 + unlist(marker$PedCheck[iErr])
      glyph_colour(gg$marker) <- glc
      glyph_type(gg$marker) <-  rep(6,nrow(marker))
      glyph_type(gg$marker)[iErr] <- 7
      message(sprintf('Number of pedigree-errors: %d',sum(iErr)))
      finished <- readline(prompt='Keep current clusters? (y/n) ')
      OK <- finished %in% c('y','Y','yes','Yes','YES','1','T','TRUE')
      close(gg)
    }else{
      OK <- FALSE
    } #end if isCalled
    if (!OK){
      iMiss <- is.na(R[i,])|is.na(SE[i,])
      iDet <- R[i,]>mNoise & !iMiss & Theta[i,]>-.5 & Theta[i,]<1.5    #Filter bad markers
      if (sum(iDet)<=ncol(BSRed)*gO$detectLim){
        message(sprintf('Marker %s failed due to many missing samples',featureNames(BSRed)[i]))
        featureData(BSRed)@data$Manual.Calls.R[[i]] <- 'FAIL'
        assayData(BSRed)$call[i,] <- NA
        assayData(BSRed)$ped.check.parents[i,] <- NA
        assayData(BSRed)$ped.check[i,] <- NA
        featureData(BSRed)@data$Cent.Deviation[[i]] <- NA
        featureData(BSRed)@data$Within.SD[i] <- NA
        featureData(BSRed)@data$HW.Chi2[i] <- NA
        featureData(BSRed)@data$HW.P[i] <- NA
        featureData(BSRed)@data$BAF.Locus1[i] <- NA
        featureData(BSRed)@data$BAF.Locus2[i] <- NA
        featureData(BSRed)@data$Call.Rate[i] <- 0
        #featureData(BSRed)@data$Max.Rotation[i] <- NA
        featureData(BSRed)@data$Ped.Errors[i] <- 0
      }else{
        sclR <- median(R[i,iDet],na.rm=TRUE)*2*gO$rPenalty                   #Factor by which R is divided before clustering
        X <- matrix(cbind(Theta[i,],R[i,]/sclR,SE[i,])[iDet,],ncol=3)
        if (is.null(gO$probsIndSE)){
          indSE <- rep(TRUE,sum(iDet))
        }else{
          indSE <- X[,3]<quantile(X[,3],probs=gO$probsIndSE,na.rm=TRUE)      #Uncertain points removed from clustering
        } #end if is.null
        #qRange <- getQRange(X[indSE,1],minScl=3,probs=gO$probsQRange)     #Initial cluster-centers
        #sConf <- distLikelyConfiguration(X[indSE,1:2],polyCent)  #Ordered list of likely classifications
        #sConf <- getCenters(X[indSE,1],gO=gO,polyCent=polyCent)
        #maxConf <- length(sConf$ix)
        while (!OK){
          marker <- data.frame(Theta=Theta[i,],R=R[i,],PedigreeID=pData(BSRed)$PedigreeID,stringsAsFactors=FALSE)[iDet,]
          gg <- ggobi(marker[,c('Theta','R')],name='marker')
          glyph_colour(gg$marker) <- 8
          clList <- unlist(polyCent$classification)#[sConf$ix])
          nCand <- length(clList)
          message('Genotype-calling categories: ')
          for (k in 1:nCand)
            message(sprintf('%d: %s',k,clList[k]))
          message(sprintf('%d: FAIL',k+1))
          rngT <- range(marker$Theta)
          message(sprintf('(Range Theta: [%.2f, %.2f])',rngT[1],rngT[2]))
          iConf <- NA
          while (!iConf %in% 1:(nCand+1)){
            str1 <- 'Please indicate the correct marker-type: (in 1:8) '
            iConf <- as.numeric(readline(prompt=str1))
          }
          if (iConf!=(nCand+1)){
            #Suggest clusters
            nCl <- polyCent$size[iConf]#[sConf$ix[iConf]]
            cntIdeal <- polyCent$centers[[iConf]]#[[sConf$ix[iConf]]]
            centers <- getSpecificCenters(X[indSE,1], classification=polyCent$classification[[iConf]], gO=gO, breaks=seq(-.25,1.25,gO$binWidth), polyCent=polyCent)
            clObj <- findPolyploidClusters(X,indSE=indSE,centers=centers,wss.update=FALSE)
            
            #Quality-check and update clusters
            marker$Call <- cntIdeal[clObj$cluster]
            gg <- manualCall(marker,cntIdeal,clList[iConf],gg,close.gg=FALSE)
            marker$Call <- gg$marker$Call
            marker$PedCheck <- gg$marker$PedCheck
            if (nCl==1){
              clObj$cluster <- gg$marker$Call-cntIdeal+1
            }else{
              clObj$cluster <- (gg$marker$Call-cntIdeal[1])/diff(range(cntIdeal))*(nCl-1)+1
            }
            glc <- gg$marker$Call*4 + 1
            glc[is.na(glc)] <- 8
            glt <- rep(6,nrow(marker))
            newCent1 <- tapply(X[,1],clObj$cluster,median,na.rm=TRUE) #NB! median
            newCent2 <- tapply(X[,2],clObj$cluster,median,na.rm=TRUE) #NB! median
            clObj$centers[as.numeric(names(newCent1)),] <- cbind(newCent1,newCent2)
            clObj$centers[-as.numeric(names(newCent1)),] <- NA
            clObj$size <- rep(0,nCl)
            newSizes <- tapply(rep(1,nrow(X)),clObj$cluster,sum)
            clObj$size[as.numeric(names(newSizes))] <- newSizes
            
            tmpRes <- rep(NA,7)
            iCl <- which(clObj$size>0)
            tmpRes[1] <- max(abs(cntIdeal[iCl]-clObj$centers[iCl,1]))   #Max cluster-centre deviation (Theta)
            tmpRes[2] <- max(tapply(X[,1],clObj$cluster,sd))  #Max within-cluster spread (Theta)
            if (nCl>1){
              tmpRes[3:6] <- as.numeric(testHardyWeinberg(sizes=clObj$size,bestConf=iConf,polyCent=polyCent,afList=gO$afList)[c(1,3:5)])
            } #end if nCl
            nClRed <- length(iCl)
            clCtrs <- matrix(clObj$centers[iCl,],nrow=nClRed)
            indOut <- matrix(FALSE,nrow=nrow(X),ncol=nCl)
            indKeep <- rep(FALSE,nrow(X))
            indOverlap <- rep(FALSE,nrow(X))
            Sratio <- rep(1,nCl)
            for (j in iCl){
              indj <- clObj$cluster%in%j
              nj <- sum(indj)
              if (nj>gO$minClLim){    #Hotellings T^2-ellipse superimposed on clusters
                ctrX <- apply(X[indj,1:2],2,median,na.rm=TRUE)
                Xj <- X[indj,1:2]-matrix(ctrX,nrow=nj,ncol=2,byrow=TRUE)
                Sj <- t(Xj)%*%Xj/(nj-1)
                Sratio[j] <- Sj[1,1]/(Sj[2,2]*gO$rPenalty)  #Cluster 'angle', given comparable axes
                Xall <- X[,1:2]-matrix(ctrX,nrow=nrow(X),ncol=2,byrow=TRUE)
                T2j <- diag(Xall%*%solve(Sj)%*%t(Xall))
                f.alpha <- T2j*(nj-2)*nj/(2*(nj-1)*(nj+1))
                pVal <- 1 - pf(f.alpha,2,nj-2)
              }else{
                pVal <- rep(0,nrow(X))
                pVal[indj] <- 1
              }
              indOut[,j] <- pVal < gO$clAlpha
              indKeep[indj] <- !indOut[indj,j]
              indOverlap[!indOut[,j]][clObj$cluster[!indOut[,j]]!=j] <- TRUE
            } #end for j
            indKeep <- indKeep &! indOverlap
            tmpRes[7] <- sum(indKeep)/ncol(BSRed)
            
            #Update calls, pedCheck, and scatterplot 
            glc[!indKeep] <- rep(8,sum(!indKeep))
            glyph_colour(gg$marker) <- glc
            marker$Call[!indKeep] <- gg$marker$Call[!indKeep] <- clObj$cluster[!indKeep] <- NA
            indSel <- marker$PedCheck>0
            marker$PedCheck[indSel] <- gg$marker$PedCheck[indSel] <- validateSingleCall(marker[indSel,],clList[iConf])
            iErr <- marker$PedCheck>0
            glyph_colour(gg$marker)[iErr] <- 5 + unlist(marker$PedCheck[iErr])
            glyph_type(gg$marker) <- glt
            glyph_type(gg$marker)[iErr] <- 7
            
            message(sprintf('Hardy-Weinberg p-value: %.2g',tmpRes[4]))
            message(sprintf('Call-rate: %.2g',tmpRes[7]))
            message(sprintf('Number of pedigree-errors: %d',sum(iErr)))
            finished <- readline(prompt='Save results and move on to next marker? (y/n) ')
          }else{
            finished <- readline(prompt='Fail marker and move on? (y/n) ')
          } #end if iConf!=8
          OK <- finished %in% c('y','Y','yes','Yes','YES','1','T','TRUE')
          close(gg)
        } #end while !OK
        if (iConf!=8){
          assayData(BSRed)$call[i,iDet][indKeep] <- cntIdeal[clObj$cluster][indKeep]
          assayData(BSRed)$call[i,!iDet] <- NA
          assayData(BSRed)$call[i,iDet][!indKeep] <- NA
          assayData(BSRed)$ped.check.parents[i,iDet] <- marker$PedCheck
          assayData(BSRed)$ped.check.parents[i,!iDet] <- NA
          assayData(BSRed)$ped.check[i,] <- NA
          assayData(BSRed)$ped.check[i,iDet][marker$PedCheck%in%0] <- TRUE
          assayData(BSRed)$ped.check[i,iDet][marker$PedCheck%in%1] <- FALSE
          featureData(BSRed)@data$Manual.Calls.R[[i]] <- polyCent$classification[[iConf]]#[[sConf$ix[iConf]]]
          featureData(BSRed)@data$Cent.Deviation[[i]] <- tmpRes[1]
          featureData(BSRed)@data$Within.SD[i] <- tmpRes[2]
          featureData(BSRed)@data$HW.Chi2[i] <- tmpRes[3]
          featureData(BSRed)@data$HW.P[i] <- tmpRes[4]
          featureData(BSRed)@data$BAF.Locus1[i] <- tmpRes[5]
          featureData(BSRed)@data$BAF.Locus2[i] <- tmpRes[6]
          featureData(BSRed)@data$Call.Rate[i] <- tmpRes[7]
          #featureData(BSRed)@data$Max.Rotation[i] <- max(Sratio)
          featureData(BSRed)@data$Ped.Errors[i] <- sum(iErr)
        }else{
          featureData(BSRed)@data$Manual.Calls.R[[i]] <- 'FAIL'
          assayData(BSRed)$call[i,] <- NA
          assayData(BSRed)$ped.check.parents[i,] <- NA
          assayData(BSRed)$ped.check[i,] <- NA
          featureData(BSRed)@data$Cent.Deviation[[i]] <- NA
          featureData(BSRed)@data$Within.SD[i] <- NA
          featureData(BSRed)@data$HW.Chi2[i] <- NA
          featureData(BSRed)@data$HW.P[i] <- NA
          featureData(BSRed)@data$BAF.Locus1[i] <- NA
          featureData(BSRed)@data$BAF.Locus2[i] <- NA
          featureData(BSRed)@data$Call.Rate[i] <- 0
          #featureData(BSRed)@data$Max.Rotation[i] <- NA
          featureData(BSRed)@data$Ped.Errors[i] <- 0
        } #end if iConf!=8
      } #end if sum(iDet)
    }else{
      featureData(BSRed)@data$Manual.Calls.R[[i]] <- featureData(BSRed)@data$Classification[[i]]
      assayData(BSRed)$ped.check.parents[i,] <- marker$PedCheck
    } #end if !OK
  } #end for i
  validObject(BSRed)
  BSRed
}



#Scales each value of the pooled SE with the arclength of the first quadrant
#semi-circle of radius 'R', as the arclength is one for all values of
#'R' in the theta-R plot. The arclength is defined by 'p' <- the norm of the
#distance measure used. If p=1 (manhattan) a=sqrt(2)*r, if p=2 (euclidean)
#a=pi/2*r, if p=inf (infinity-norm) a=2r.
findSeTheta <- function(pooledSE.raw,R,dist='manhattan',pNorm=NULL){
  integrand <- function(th,pNorm){
    2/pNorm*sqrt(cos(th)^(4/pNorm-2)*sin(th)^2+sin(th)^(4/pNorm-2)*cos(th)^2)
  }
  switch(dist,
         euclidean={scale <- pi*R/2},
         manhattan={scale <- sqrt(2)*R},
         minkowski={
           unitLength <- integrate(integrand,lower=0,upper=pi/2,pNorm=pNorm)
           scale <- unitLength$value*R
         })
  pooledSE.theta <- pooledSE.raw/scale
  pooledSE.theta
}
  


#Generates a list of possible cluster-alternatives
generatePolyCenters <- function(ploidy){
  centers <- list()
  classification <- list()
  if (ploidy=='di'){
    centers[[1]] <- 0               #Only A's
    centers[[2]] <- 1               #Only B's
    centers[[3]] <- c(0,.5,1)            #Regular SNP
    classification[1] <- 'MONO-a'
    classification[2] <- 'MONO-b'
    classification[3] <- 'SNP'
  }else if(ploidy=='tetra'){
    centers[[1]] <- 0               #Only A's
    centers[[2]] <- .5              #PSV or MSV with small minor AF
    centers[[3]] <- 1               #Only B's
    centers[[4]] <- c(0,.5,1)            #Regular SNP
    centers[[5]] <- c(0,.25,.5)          #MSV with one paralogue fixed AA
    centers[[6]] <- c(.5,.75,1)          #MSV with one paralogue fixed BB
    centers[[7]] <- c(0,.25,.5,.75,1)    #MSV with segregation in both paralogues
    classification[1] <- 'MONO-a'
    classification[2] <- 'PSV'
    classification[3] <- 'MONO-b'
    classification[4] <- 'SNP'
    classification[5] <- 'MSV-a'
    classification[6] <- 'MSV-b'
    classification[7] <- 'MSV-5'
  }else if(ploidy=='tetra.red'){
    centers[[1]] <- 0               #Only A's
    centers[[2]] <- .5              #PSV or MSV with small minor AF
    centers[[3]] <- 1               #Only B's
    centers[[4]] <- c(0,.5,1)            #Regular SNP
    centers[[5]] <- c(0,.25,.5)          #MSV with one paralogue fixed AA
    centers[[6]] <- c(.5,.75,1)          #MSV with one paralogue fixed BB
    classification[1] <- 'MONO-a'
    classification[2] <- 'PSV'
    classification[3] <- 'MONO-b'
    classification[4] <- 'SNP'
    classification[5] <- 'MSV-a'
    classification[6] <- 'MSV-b'
  }else{
    stop("Allowed values for ploidy are 'di', 'tetra', or 'tetra.red'")
  }
  polyCent <- list()
  polyCent$centers <- centers
  polyCent$classification <- classification
  polyCent$size <- sapply(polyCent$centers,length)
  polyCent
}



#Estimates likely range and cluster-centres
getCenters <- function(theta,gO=setGenoOptions(),breaks=seq(-.25,1.25,gO$binWidth),polyCent=generatePolyCenters(ploidy=gO$ploidy)){
##### for testing:
  #breaks <- seq(-.25,1.25,.1)
  #theta <- assayData(BSRed)$theta[125,]; theta <- theta[theta>breaks[1]&theta<max(breaks)]
  #hX <- hist(theta,breaks=breaks,plot=TRUE,ylim=c(0,50));abline(h=gO$minBin-1)
#####
  
  theta <- theta[theta>breaks[1] & theta<max(breaks)]
  nSamp <- length(theta)
  polyCl <- findClusters(theta=theta,breaks=breaks,minBin=gO$minBin)
  if (polyCl$nCl<3){
    breaks1 <- seq(breaks[1],max(breaks),gO$binWidth/2)
    polyCl1 <- findClusters(theta=theta,breaks=breaks1,minBin=gO$minBin)
    if (all(polyCl1$clSizes>nSamp/100)){
      polyCl <- polyCl1
    }
  }
  centers <- list()
  if (polyCl$nCl==0){
    ix <- NULL
  }else if (gO$ploidy%in%'tetra'){
    i1 <- polyCent$size==1
    i3 <- polyCent$size==3
    i5 <- polyCent$size==5
    if (polyCl$nCl==1){
      vec1 <- unlist(polyCent$centers[i1])
      centers[[1]] <- sum(polyCl$clPeaks*polyCl$clSizes)/sum(polyCl$clSizes)
      sqDev1 <- (vec1 - centers[[1]])^2/c(2,1,2)   #Disfavour PSVs
      cent3 <- polyCent$centers[[4]]
      indClose <- which(cent3 %in% (round(polyCl$clPeaks*2)/2))
      cent3[indClose] <- polyCl$clPeaks
      centers[[2]] <- cent3
      ix <- c(which(i1)[order(sqDev1)[1]],4)
    }else if (polyCl$nCl==2){
      ix <- numeric(3)
      cent5 <- polyCent$centers[[7]]
      indClose <- cent5 %in% (round(polyCl$clPeaks*4)/4)
      if (sum(indClose)==1){
        cDist <- sapply(polyCl$clPeaks,function(x,cent5) (x-cent5)^2, cent5)
        orderC <- order(cDist[indClose,])
        cDist[indClose,orderC[2]] <- 10
        indClose[order(cDist[,orderC[2]])[1]] <- TRUE
      }
      indClose <- which(indClose)
      closerAA <- sum((polyCl$clPeaks)^2) < sum((polyCl$clPeaks-1)^2)
      
      if (sum(c(1,3,5)%in%indClose)==2){  #1-3, 3-5, 1-5
        cent5[indClose] <- polyCl$clPeaks
        if (closerAA){
          ii2 <- 5
        }else{
          ii2 <- 6
        }
        cent3 <- polyCent$centers[[ii2]]
        indClose1 <- which(rank((cent3-polyCl$clPeaks[order(polyCl$clSizes)[2]])^2)==1)
        cent3[indClose1] <- polyCl$clPeaks[order(polyCl$clSizes)[2]]
        centers[1:2] <- list(cent5[c(1,3,5)],cent3)
        ix[1:2] <- c(4,ii2)
      }else if(sum(1:3%in%indClose)==2){   #1-2, 2-3
        cent5[indClose] <- polyCl$clPeaks
        cent3 <- polyCent$centers[[4]]
        cent3[1:2] <- polyCl$clPeaks
        centers[1:2] <- list(cent5[1:3],cent3)
        ix[1:2] <- c(5,4)
      }else if(sum(3:5%in%indClose)==2){   #4-5, 3-4
        cent5[indClose] <- polyCl$clPeaks
        cent3 <- polyCent$centers[[4]]
        cent3[2:3] <- polyCl$clPeaks
        centers[1:2] <- list(cent5[3:5],cent3)
        ix[1:2] <- c(6,4)
      }else{                               #1-4, 2-4, 2-5
        cent3 <- polyCent$centers[[4]]
        indClose <- which(cent3 %in% (round(polyCl$clPeaks*2)/2))
        if (length(indClose)==1){
          if (closerAA){                                #clusters closer to theta=0
            indClose <- which(cent3 %in% (floor(pmax(polyCl$clPeaks,0)*2)/2))
          }else{                                        #clusters closer to theta=1
            indClose <- which(cent3 %in% (ceiling(pmin(polyCl$clPeaks,1)*2)/2))
          }
        }
        cent3[indClose] <- polyCl$clPeaks
        if (closerAA){
          ii2 <- 5
        }else{
          ii2 <- 6
        }
        cent3a <- polyCent$centers[[ii2]]
        indClose1 <- which(rank((cent3a-polyCl$clPeaks[order(polyCl$clSizes)[2]])^2)==1)
        cent3a[indClose1] <- polyCl$clPeaks[order(polyCl$clSizes)[2]]
        centers[1:2] <- list(cent3,cent3a)
        ix[1:2] <- c(4,ii2)
      }
      vec1 <- unlist(polyCent$centers[i1])
      centers[[3]] <- sum(polyCl$clPeaks*polyCl$clSizes)/sum(polyCl$clSizes)
      sqDev1 <- (vec1 - centers[[3]])^2
      ix[3] <- which(i1)[order(sqDev1)[1]]
    }else if (polyCl$nCl==3){
      ix <- numeric(3)
      arr3 <- t(array(unlist(polyCent$centers[i3]),dim=c(sum(i3),3)))
      indx <- apply(arr3,1,function(x,clPeaks,lim) max(abs(x-clPeaks))<=lim,polyCl$clPeaks,gO$devCentLim)  #Don't attempt if devCentLim is exceeded
      if (any(indx)){
        arr3 <- matrix(arr3[indx,],nrow=sum(indx),byrow=TRUE)
        sqDev3 <- apply(arr3,1,function(x) sum((x-polyCl$clPeaks)^2))/c(2,1,1)[indx] #alt: apply(arr3,1,function(x) diff(range(x)))
        ix[1] <- which(i3)[indx][order(sqDev3)[1]]
        centers[[1]] <- polyCl$clPeaks
      }else if (diff(range(theta))<.25){
        vec1 <- unlist(polyCent$centers[i1])
        centers[[1]] <- sum(polyCl$clPeaks*polyCl$clSizes)/sum(polyCl$clSizes)
        sqDev1 <- (vec1 - centers[[1]])^2/c(2,1,2)   #Disfavour PSVs
        ix[1] <- which(i1)[order(sqDev1)[1]]
      }else{
        ix[1] <- 7
        centers[[1]] <- c(0,polyCl$clPeaks,1)
      }
      breaks1 <- seq(breaks[1],max(breaks),gO$binWidth/2)
      polyCl1 <- findClusters(theta=theta,breaks=breaks1,minBin=gO$minBin)
      if (polyCl1$nCl==5){
        centers[[2]] <- polyCl1$clPeaks
      }else{
        centers[[2]] <- seq(polyCl$clPeaks[1],polyCl$clPeaks[3],diff(range(polyCl$clPeaks))/4)
      }
      ix[2] <- 7
      
      diffPeaks <- diff(polyCl$clPeaks)
      iSame <- c(diffPeaks<=min(diffPeaks),FALSE)
      if (sum(iSame)==1){
        cent3 <- polyCent$centers[[4]]
        clPeaksRed <- polyCl$clPeaks[!iSame]
        indx <- which(iSame)
        clPeaksRed[indx] <- mean(polyCl$clPeaks[indx:(indx+1)])
        indClose <- which(cent3 %in% (round(clPeaksRed*2)/2))
        if (length(indClose)==1){
          closerAA <- sum((clPeaksRed)^2) < sum((clPeaksRed-1)^2)
          if (closerAA){                                #clusters closer to theta=0
            indClose <- c(1,2)
          }else{                                        #clusters closer to theta=1
            indClose <- c(2,3)
          }
        }
        cent3[indClose] <- clPeaksRed
        centers[[3]] <- cent3
        ix[3] <- 4
      }else{
        vec1 <- unlist(polyCent$centers[i1])
        centers[[3]] <- sum(polyCl$clPeaks*polyCl$clSizes)/sum(polyCl$clSizes)
        sqDev1 <- (vec1 -  centers[[3]])^2
        ix[3] <- which(i1)[order(sqDev1)[1]]
      }
    }else if (polyCl$nCl==4){
      cent5 <- polyCent$centers[[7]]
      #minDist <- tapply(cent5,1:5,function(x,clPeaks) min((x-clPeaks)^2),polyCl$clPeaks)
      minDist <- sapply(cent5,function(x,clPeaks) min((x-clPeaks)^2),polyCl$clPeaks)
      indClose <- !rank(minDist,ties.method='max')==5
      if (sum(indClose)<4){
        if (sum((polyCl$clPeaks)^2) < sum((polyCl$clPeaks-1)^2)){   #clusters closer to theta=0
          indClose <- c(rep(TRUE,4),FALSE)
        }else{                                        #clusters closer to theta=1
          indClose <- c(FALSE,rep(TRUE,4))
        }
      }
      cent5[indClose] <- polyCl$clPeaks
      centers[[1]] <- cent5
      breaks1 <- seq(breaks[1],max(breaks),gO$binWidth/2)
      polyCl1 <- findClusters(theta=theta,breaks=breaks1,minBin=gO$minBin)
      if (polyCl1$nCl==5){
        centers[[2]] <- polyCl1$clPeaks
      }else{
        centers[[2]] <- seq(polyCl$clPeaks[1],polyCl$clPeaks[4],diff(range(polyCl$clPeaks))/4)
      }
      diffPeaks <- diff(polyCl$clPeaks)
      iSame <- c(diffPeaks<=min(diffPeaks),FALSE)
      if (sum(iSame)==1){
        clPeaksRed <- polyCl$clPeaks[!iSame]
        indx <- which(iSame)
        clPeaksRed[indx] <- mean(polyCl$clPeaks[indx:(indx+1)])
      }else{
        clPeaksRed <- c(polyCl$clPeaks[1],mean(polyCl$clPeaks[2:3]),polyCl$clPeaks[4])
      }
      centers[[3]] <- clPeaksRed
      arr3 <- t(array(unlist(polyCent$centers[i3]),dim=c(sum(i3),3)))
      sqDev3 <- apply(arr3,1,function(x) sum((x-clPeaksRed)^2))/c(2,1,1) #alt: apply(arr3,1,function(x) diff(range(x)))
      ix <- c(7,7,which(i3)[order(sqDev3)[1]])
    }else if (polyCl$nCl==5){
      centers[[1]] <- polyCl$clPeaks
      #centers[[2]] <- seq(0,1,.25)
      centers[[2]] <- c(polyCl$clPeaks[1],mean(polyCl$clPeaks[2:4]),polyCl$clPeaks[5])
      ix <- c(7,4)
    }else{
      ix <- c(4,7)
      centers <- polyCent$centers[ix]
    } #end if polyCl$nCl
  }else{
    ix <- 1:length(polyCent$classification)
    centers <- polyCent$centers
    warning('Only ploidy="tetra" is implemented. A non-ranked list of starting points
is returned for other ploidy, however accurate clustering is NOT expected')
  #}else if (gO$ploidy%in%'di'){
  #  i1 <- polyCent$size==1
  #  i3 <- polyCent$size==3
  #  if (polyCl$nCl==1){
  #    
  #  }else if (polyCl$nCl==2){
  #    
  #  }else if (polyCl$nCl==3){
  #    
  #  }else{
  #    
  #  } #end if polyCl$nCl==1
  } #end if polyCl$nCl==0
  sConf <- list(ix=ix,centers=centers)
  sConf
}



getSpecificCenters <- function(theta,classification,gO=setGenoOptions(),breaks=seq(-.25,1.25,gO$binWidth),polyCent=generatePolyCenters(ploidy=gO$ploidy)){
##### for testing:
  #breaks <- seq(-.25,1.25,.1)
  #theta <- assayData(BSRed)$theta[8,]; theta <- theta[theta>breaks[1]&theta<max(breaks)]
  #hX <- hist(theta,breaks=breaks,plot=TRUE,ylim=c(0,50));abline(h=gO$minBin-1)
#####
  
  theta <- theta[theta>breaks[1] & theta<max(breaks)]
  nSamp <- length(theta)
  polyCl <- findClusters(theta=theta,breaks=breaks,minBin=gO$minBin)
  iConf <- which(polyCent$classification %in% classification)
  if (polyCl$nCl==polyCent$size[iConf]){
    return(polyCl$clPeaks)
  }else if (polyCl$nCl<polyCent$size[iConf]){
    breaks1 <- seq(breaks[1],max(breaks),gO$binWidth/2)
    polyCl1 <- findClusters(theta=theta,breaks=breaks1,minBin=gO$minBin)
    if (polyCl1$nCl==polyCent$size[iConf]){
      return(polyCl1$clPeaks)
    }
  }else{
    breaks1 <- seq(breaks[1],max(breaks),gO$binWidth*2)
    polyCl1 <- findClusters(theta=theta,breaks=breaks1,minBin=gO$minBin)
    if (polyCl1$nCl==polyCent$size[iConf]){
      return(polyCl1$clPeaks)
    }
  }
  return(polyCent$centers[[iConf]])
}


  
#Suggests clusters based on histograms
findClusters <- function(theta,breaks=seq(-.25,1.25,.05),minBin=2,plot=FALSE){
  hX <- hist(theta,breaks=breaks,plot=plot)
  spec <- hX$counts
  spec[spec<minBin] <- 0
  spec <- c(0,spec,0)
  indPeaks <- which(apply(embed(spec,3),1,function(x) rank(x)[2]>2))
  clPeaks <- hX$mids[indPeaks]
  iSame <- c(diff(indPeaks)==1,FALSE)
  if (any(iSame)){
    clPeaksRed <- clPeaks[!iSame]
    indx <- which(iSame)
    for (i in 1:sum(iSame)){
      clPeaksRed[indx[1]+1-i] <- mean(clPeaks[indx[1]:(indx[1]+1)])
      indx <- indx[-1]
    }
    clPeaks <- clPeaksRed
  }
  nCl <- length(clPeaks)
  clSizes <- apply(matrix(embed(c(0,hX$counts,0),3)[indPeaks[!iSame],],nrow=nCl),1,sum)
  polyCl <- list(clPeaks=clPeaks,clSizes=clSizes,nCl=nCl)
  polyCl
}



findPolyploidClusters <- function(X,indSE=rep(TRUE,nrow(X)),centers,plot=FALSE,wss.update=TRUE,...){
  clObj <- NULL
  nCl <- length(centers)
  if (sum(indSE)>=nCl+1){   #Requires more samples than groups
    breaks <- c(-2,-.25,centers[-nCl]+diff(centers)/2,1.25,2)
    hX <- hist(X[indSE,1],breaks=breaks,plot=plot,...)
    indCall <- hX$counts[-c(1,nCl+2)]>0
  }else{
    stop('k-means algorithm requires more samples than cluster-centers')
  }
  centers <- cbind(centers,rep(NA,nCl))
  if (sum(indCall)>1){
    for (i in which(indCall)){
      indi <- X[indSE,1]>hX$breaks[i+1] & X[indSE,1]<=hX$breaks[i+2]+1e-7  #Accounts for numerical tolerance of 'hist'
      centers[i,2] <- median(X[indSE,2][indi],na.rm=TRUE)
    }
    clObj <- try(kmeans(X[indSE,1:2],centers[indCall,]),silent=TRUE)
    if (inherits(clObj, "try-error")){
      message(clObj)
      warning('Hartigan-Wong clustering failed for one marker (kmeans).\n  Will attempt MacQueen-algorithm instead.')
      clObj <- kmeans(X[indSE,1:2],centers[indCall,],algorithm='MacQueen')
    }
  }else if (sum(indCall)==1){
    clObj$centers <- matrix(apply(X[indSE,1:2],2,median,na.rm=TRUE),nrow=1)
    clObj$cluster <- rep(1,sum(indSE))
  }else{
    stop('No clusters found')
  }

  ## Assign left-outs to clusters
  clCluster <- rep(NA,nrow(X))
  clCluster[indSE] <- which(indCall)[clObj$cluster]
  centers[indCall,] <- clObj$centers
  centers[!indCall,2] <- mean(X[!indSE,2],na.rm=TRUE)
  if (!all(indSE)){
    #dAll <- (matrix(clObj$centers[,1],nrow=sum(!indSE),ncol=sum(indCall),byrow=TRUE) - matrix(X[!indSE,1],nrow=sum(!indSE),ncol=sum(indCall)))^2 + (matrix(clObj$centers[,2],nrow=sum(!indSE),ncol=sum(indCall),byrow=TRUE) - matrix(X[!indSE,2],nrow=sum(!indSE),ncol=sum(indCall)))^2
    dAll <- (matrix(centers[,1],nrow=sum(!indSE),ncol=nCl,byrow=TRUE) - matrix(X[!indSE,1],nrow=sum(!indSE),ncol=nCl))^2 + (matrix(centers[,2],nrow=sum(!indSE),ncol=nCl,byrow=TRUE) - matrix(X[!indSE,2],nrow=sum(!indSE),ncol=nCl))^2
    clCluster[!indSE] <- apply(dAll,1,function(x) which(rank(x,ties.method='random')==1))
  }
  #for (a in 1:nCl){
  #  if (!indCall[a]){
  #    clCluster[clCluster>=a] <- clCluster[clCluster>=a] + 1
  #  }
  #}
  clObj$cluster <- clCluster

  ## Update cluster-centers etc.
  clObj$centers <- matrix(NA,nrow=nCl,ncol=2)
  newCent1 <- tapply(X[,1],clCluster,median,na.rm=TRUE) #NB! median
  newCent2 <- tapply(X[,2],clCluster,median,na.rm=TRUE) #NB! median
  clObj$centers[as.numeric(names(newCent1)),] <- cbind(newCent1,newCent2)
  
  clObj$size <- rep(0,nCl)
  newSizes <- tapply(clCluster,clCluster,length)
  clObj$size[as.numeric(names(newSizes))] <- newSizes

  if (sum(!indSE)>0){
    clObj$withinss <- rep(NA,nCl)
    if (wss.update){   #Save time by not calculating withinss
      ws1 <- tapply(X[,1],clCluster,function(x) (x-rep(mean(x,na.rm=TRUE),length(x)))^2)
      ws2 <- tapply(X[,2],clCluster,function(x) (x-rep(mean(x,na.rm=TRUE),length(x)))^2)
      newWss <- rep(NA,nCl)
      for (i in 1:nCl){
        newWss[i] <- sum(ws1[[i]] + ws2[[i]])
      }
      clObj$withinss[as.numeric(names(ws1))] <- newWss
    }
  }
  class(clObj) <- 'kmeans'
  clObj
}



#Null hypothesis: The population is in Hardy-Weinberg frequencies.
#Different allele-frequencies for testing tetraploid loci given in afList.
#A value of afList=.5 means that the AF is the same at both loci.
testHardyWeinberg <- function(sizes,bestConf,polyCent=generatePolyCenters(ploidy='di'),afList=seq(0,.5,.05)){
  if (is.character(bestConf))
    bestConf <- which(polyCent$classification==bestConf)
  nPheno <- polyCent$size[bestConf]
  nTot <- sum(sizes)
  af <- rep(NA,2)
  if (nPheno==1){
    chiSq <- 0
    df <- 1   #dummy-variable
    if (bestConf==1){
      af[1] <- 1
    }else if (bestConf==2){
      af <- 0:1
    }else if (bestConf==3){
      af[1] <- 0
    }
  }else if (nPheno==3){
    p <- (2*sizes[1]+sizes[2])/(2*nTot)
    if (any(p==c(0,1))){
      chiSq <- 0
    }else{
      q <- 1-p
      expCounts <- rep(NA,nPheno)
      expCounts[1] <- p^2*nTot          #AA
      expCounts[2] <- 2*p*q*nTot        #AB
      expCounts[3] <- q^2*nTot          #BB
      chiSq <- sum((sizes-expCounts)^2/expCounts)
    }
    df <- 1
    af[1] <- p
    if (bestConf==5){
      af[2] <- 1
    }else if (bestConf==6){
      af[2] <- 0
    }
  }else if (nPheno==5){
    p <- (4*sizes[1]+3*sizes[2]+2*sizes[3]+sizes[4])/(4*nTot)  #'sizes' contains the 5 cluster-sizes
    if (any(p==c(0,1))){
      chiSq <- 0
      af <- rep(p,2)
    }else{
      iLim <- min(1-1/(2*p),1/(2*p))
      if (any(afList<iLim))  #Negative af's not allowed
        afList <- c(iLim,afList[afList>iLim])
      nTry <- length(afList)   #Try different allowed weights between p1 and p2
      expCounts <- matrix(nrow=nTry,ncol=nPheno)
      chiSqs <- rep(Inf,nTry)
      afs <- matrix(nrow=nTry,ncol=2)
      for (i in 1:nTry){
        p1 <- 2*p*afList[i]
        p2 <- 2*p*(1-afList[i])   #mean(AFs)=p (AFs sum to 2*p)
        if (all(c(p1,p2)<=1)){
          q1 <- 1-p1
          q2 <- 1-p2
          expCounts[i,1] <- p1^2*p2^2*nTot              #AAAA
          expCounts[i,2] <- 2*p1*p2*(p1*q2+p2*q1)*nTot      #AAAB
          expCounts[i,3] <- (p1^2*q2^2+p2^2*q1^2+4*p1*q1*p2*q2)*nTot #AABB
          expCounts[i,4] <- 2*q1*q2*(p1*q2+p2*q1)*nTot      #ABBB
          expCounts[i,5] <- q1^2*q2^2*nTot          #BBBB
          indi <- expCounts[i,]>0
          chiSqs[i] <- sum((sizes[indi]-expCounts[i,indi])^2/expCounts[i,indi])
          afs[i,] <- c(p1,p2)
        }
      }
      sChi <- sort(chiSqs,index.return=TRUE)
      chiSq <- chiSqs[sChi$ix[1]]
      af <- afs[sChi$ix[1],]
    }
    df <- 3    #Number of genotypes/clusters - number of alleles
  }else{
    stop(paste('HW-testing not implemented for',nPheno,'cluster(s)'))
  }
  pVal <- pchisq(chiSq,df=df,lower.tail=FALSE)
  HWstats <- data.frame(chiSq=chiSq,df=df,pVal=pVal,af1=af[1],af2=af[2])
  HWstats
}



#Validates BSRed, returns updated BSRed
#PedigreeID given as '<p><mmm><fff><oo>', where 'p', 'mmm', 'fff' and 'oo' are unique
#identifiers for population, mother, father, and individual within full-sib group,
#respectively. '000' means founding parent, '999' means unknown.
validateCallsPedigree <- function(BSRed){
  assayData(BSRed)$ped.check <- matrix(nrow=nrow(BSRed),ncol=ncol(BSRed),dimnames=list(featureNames(BSRed),sampleNames(BSRed)))
  if (!'PedigreeID' %in% colnames(pData(BSRed))){
    stop('PedigreeID not provided!')
  }
  sampIDs <- pData(BSRed)$PedigreeID
  indFamily <- unique(substr(sampIDs,1,7)[!sampIDs%in%''])
  indFounding <- substr(indFamily,2,4)=='000' | substr(indFamily,5,7)=='000'
  indFamily <- indFamily[!indFounding]
  nFull <- length(indFamily)
  
  #Establish correct classification and adjust MSV3-calls to diploid
  callType <- fData(BSRed)$Classification
  if (!is.null(fData(BSRed)$Manual.Calls.R)){
    iManual <- !fData(BSRed)$Manual.Calls.R%in%''
    callType[iManual] <- fData(BSRed)$Manual.Calls.R[iManual]
  }
  callMatr <- assayData(BSRed)$call
  callMatr[callType%in%'MSV-a',] <- callMatr[callType%in%'MSV-a',]*2
  callMatr[callType%in%'MSV-b',] <- callMatr[callType%in%'MSV-b',]*2-1

  for (j in 1:nFull){
    jf <- sampIDs %in% paste(substr(indFamily[j],1,1),'000',substr(indFamily[j],5,7),'00',sep='')
    jm <- sampIDs %in% paste(substr(indFamily[j],1,4),'00000',sep='')
    j.fam <- substr(sampIDs,1,7) %in% indFamily[j] &! (jf|jm)
    pCalls <- matrix(callMatr[,jf|jm],nrow=nrow(BSRed))
    oCalls <- matrix(callMatr[,j.fam],nrow=nrow(BSRed))
    if (ncol(pCalls)==1){
      assayData(BSRed)$ped.check[,j.fam] <- oCalls>=rep(pCalls,sum(j.fam))-.5 & oCalls<=rep(pCalls,sum(j.fam))+.5
    }else if (ncol(pCalls)==2){
      indNA <- apply(pCalls,1,function(x) any(is.na(x)))
      indLH <- pCalls[,1]<=.5 & pCalls[,2]>=.5 & !indNA
      assayData(BSRed)$ped.check[indLH,j.fam] <- oCalls[indLH,]<=rep(pCalls[indLH,1],sum(j.fam))+.5 & oCalls[indLH,]>=rep(pCalls[indLH,2],sum(j.fam))-.5
      indHL <- pCalls[,2]<=.5 & pCalls[,1]>=.5 & !indLH & !indNA
      assayData(BSRed)$ped.check[indHL,j.fam] <- oCalls[indHL,]<=rep(pCalls[indHL,2],sum(j.fam))+.5 & oCalls[indHL,]>=rep(pCalls[indHL,1],sum(j.fam))-.5
      indL <- which(apply(pCalls,1,function(x) all(x<.5)))
      assayData(BSRed)$ped.check[indL,j.fam] <- oCalls[indL,]<=apply(matrix(pCalls[indL,],nrow=length(indL),ncol=2),1,sum)
      indH <- which(apply(pCalls,1,function(x) all(x>.5)))
      assayData(BSRed)$ped.check[indH,j.fam] <- oCalls[indH,]>=apply(matrix(pCalls[indH,],nrow=length(indH),ncol=2),1,sum)-1
      indNA1 <- apply(pCalls,1,function(x) sum(is.na(x))==1)
      pRed <- t(pCalls[indNA1,])[!is.na(t(pCalls[indNA1,]))]
      assayData(BSRed)$ped.check[indNA1,j.fam] <- oCalls[indNA1,]>=rep(pRed,sum(j.fam))-.5 & oCalls[indNA1,]<=rep(pRed,sum(j.fam))+.5
    }else if (ncol(pCalls)>2){
      stop('More than two parents. Genotyped repeatedly, or ped-error?')
    }
  }
  nErrors <- apply(assayData(BSRed)$ped.check,1,function(x) sum(!x,na.rm=TRUE))
  fData <- data.frame(Ped.Errors=nErrors)
  fMetadata <- data.frame(labelDescription='Number of pedigree-errors, summed over offspring')
  featureData <- new('AnnotatedDataFrame', data=fData, varMetadata=fMetadata, dimLabels=c('featureNames','featureColumns'))
  indx <- !colnames(fData(BSRed))%in%'Ped.Errors'
  featureData(BSRed) <- combine(featureData(BSRed)[,indx], featureData)
  validObject(BSRed)
  BSRed
}



#Calculates the ratio of non-FAILed markers which are called for each array.
#If ped-error flag is TRUE, calls violating pedigree count as missing. 
countFailedSNP <- function(BSRed,inclPedErrors=TRUE){
  classMatr <- assayData(BSRed)$call
  if (inclPedErrors & ! is.null(assayData(BSRed)$ped.check)){
    classMatr[!assayData(BSRed)$ped.check] <- NA
    lStr <- 'Ratio of passed markers (ped-errors count as failed)'
  }else{
    lStr <- 'Ratio of passed markers (ped-errors not considered)'
  }
  iPass <- !fData(BSRed)$Classification %in% c('FAIL','MONO-filt')
  nErrors <- apply(classMatr[iPass,],2,function(x) sum(is.na(x)))
  pData <- data.frame(passRatio = 1-nErrors/sum(iPass))
  pMetadata <- data.frame(labelDescription = lStr)
  phenoData <- new('AnnotatedDataFrame', data=pData, varMetadata=pMetadata, dimLabels=c('sampleNames','sampleColumns'))
  phenoData(BSRed) <- combine(phenoData(BSRed),phenoData)
  validObject(BSRed)
  BSRed
}



#Plotting function
plotGenotypes <- function(BSRed,markers=1:min(nrow(BSRed),64),indHighlight=NULL,ploidy='tetra',indicate.SE=FALSE,retFrames=FALSE,nC=NULL,mai=NULL,mNoise=NULL,main=NULL){
  nSNPs <- length(markers)
  pp <- par()
  if (is.null(nC)){
    nC <- ceiling(sqrt(nSNPs))
  }
  nR <- ceiling(nSNPs/nC)
  if (is.null(mai)){
    mai <- c(1,1,1,.5)/sqrt(min(nC,nR))
  }
  par(mfrow=c(nR,nC),mai=mai,bg='white')
  
  Theta <- matrix(assayData(BSRed)$theta[markers,],nrow=nSNPs)
  R <- matrix(assayData(BSRed)$intensity[markers,],nrow=nSNPs)
  Call <- matrix(assayData(BSRed)$call[markers,],nrow=nSNPs)
  if (indicate.SE)
    callProb <- 1/matrix(assayData(BSRed)$SE[markers,],nrow=nSNPs)
  pedCheckP <- assayData(BSRed)$ped.check.parents[markers,]
  if (!is.null(pedCheckP)){  #Use 'ped.check.parents' only if valid for all markers
    classR <- fData(BSRed)$Manual.Calls.R[markers]
    if (is.null(classR)){
      ii <- rep(TRUE,nSNPs)
    }else{
      ii <- !classR%in%'FAIL'
    }
    useP <- !any(apply(is.na(pedCheckP[ii,]),1,all))
  }else{
    useP <- FALSE
  }
  if (!useP){ #Use 'ped.check' rather than 'ped.check.parents'
    pedCheck <- assayData(BSRed)$ped.check[markers,]
    if (!is.null(pedCheck)){
      pedCheck <- matrix(pedCheck,nrow=nSNPs)
      pedCheck[is.na(pedCheck)] <- TRUE
      pedParents <- matrix(0,nrow=nSNPs,ncol=ncol(BSRed))
      pedParents[!pedCheck] <- 1
    }else{
      pedParents <- NULL
    }
  }else{
    pedParents <- NULL
  }
  if (is.null(mNoise)){
    mNoise <- rep(mean(pData(BSRed)$noiseIntensity,na.rm=TRUE),nSNPs)
  }
  hetLine <- c(.5,.5,0,1e6)
  cex <- 50/max(50,min(ncol(BSRed),200))   #No larger than one, no smaller than 1/4
  if (!is.null(indHighlight))
    cexhl <- 50/max(50,min(length(indHighlight),100))
  
  switch(ploidy,
         di={
           stop('ploidy=\'di\' not implemented')
         },
         tetra={
           palette(c('red','magenta','blue','cyan','seagreen','black','purple','lightblue'))
           indTriplets <- Call==1/3 | Call==2/3 | is.na(Call)
           Call[!indTriplets] <- Call[!indTriplets]*4 + 1
           Call[is.na(Call)] <- 6
           Call[Call==1/3] <- 7
           Call[Call==2/3] <- 8
           for (i in 1:nSNPs){
             xylim <- c(-0.6,1.6,0,max(R[i,],na.rm=TRUE))
             plot.new()
             plot.window(xylim[1:2],xylim[3:4]); axis(1); axis(2)
             if (indicate.SE){
               points(Theta[i,],R[i,],pch=16,col=Call[i,],cex=callProb[i,]*2/max(callProb[1,],na.rm=TRUE))
             }else if (!is.null(indHighlight)){
               points(Theta[i,-indHighlight],R[i,-indHighlight],pch=16,cex=cex,col=Call[i,-indHighlight])
               points(Theta[i,indHighlight],R[i,indHighlight],pch=16,cex=cexhl,col='yellow2')
             }else{
               points(Theta[i,],R[i,],pch=16,cex=cex,col=Call[i,])
               if (useP){
                 points(Theta[i,pedCheckP[i,]%in%1],R[i,pedCheckP[i,]%in%1],pch=4,lwd=3,col='yellow2')
                 points(Theta[i,pedCheckP[i,]%in%2],R[i,pedCheckP[i,]%in%2],pch=1,lwd=3,col='yellow2')
               }else if (!is.null(pedCheck)){
                 points(Theta[i,!pedCheck[i,]],R[i,!pedCheck[i,]],pch=4,lwd=3,col='yellow2')
                 famInd <- unique(substr(pData(BSRed)$PedigreeID[!pedCheck[i,]],1,7))
                 for (j in 1:length(famInd)){
                   jf <- pData(BSRed)$PedigreeID %in% paste(substr(famInd[j],1,1),'000',substr(famInd[j],5,7),'00',sep='')
                   jm <- pData(BSRed)$PedigreeID %in% paste(substr(famInd[j],1,4),'00000',sep='')
                   points(Theta[i,jf|jm],R[i,jf|jm],pch=1,lwd=3,col='yellow2')
                   pedParents[i,jf|jm] <- 2
                 }
               }
             }
             if (is.null(main)){
               title(main=featureNames(BSRed)[markers[i]],xlab='theta',ylab='intensity')
             }else{
               title(main=main[i],xlab='theta',ylab='intensity')
             }
             lines(c(0,0),xylim[3:4],col='gray')
             lines(c(1,1),xylim[3:4],col='gray')
             lines(hetLine[1:2],hetLine[3:4],col='gray')
             lines(xylim[1:2],rep(mNoise[i],2),col='red')
           }},
         pooled={
           stop('ploidy=\'pooled\' not implemented')
         }
         ) #end switch
  par(pp[c('mfrow','mai')])
  if (retFrames){
    xy <- list()
    for (i in 1:nSNPs){
      xy[[i]] <- data.frame(Theta=Theta[i,],R=R[i,],Call=Call[i,],pedErrors=pedParents[i,])
    }
    xy
  }
}



#Interactively define clusters for a given marker.
#Optional parameter gg is an instance of a GObject
manualCall <- function(marker,cntIdeal,classification,gg=NULL,close.gg=TRUE){
  nCl <- length(cntIdeal)
  if (is.null(marker$Call))
    marker$Call <- rep(NA,nrow(marker))
  if (is.null(marker$PedCheck) & !is.null(marker$PedigreeID))
    marker$PedCheck <- validateSingleCall(marker,classification)
  if (is.null(gg)){
    gg <- ggobi(marker[,c('Theta','R')],name='marker')
  }else{
    if (!'GGobi'%in%class(gg) | nrow(gg$marker)!=nrow(marker))
      stop('Arguments \'marker\' and \'gg\' not compatible')
  }
  dd <- displays(gg)[[1]]
  imode(dd) <- 'Brush'
  str.3 <- sprintf('Modify which cluster? (in [1,%d], or use "8" to assign NAs) ',nCl)
  str.8 <- sprintf("Brush calls to erase, OR
press <1> to return: ")
  str.end <- 'Satisfied with assignments? (y/n) '
  glc <- unlist(as.numeric(marker$Call)*4+1)
  glc[is.na(glc)] <- rep(8,sum(is.na(glc)))
  glyph_colour(gg$marker) <- glc
  glt <- glyph_type(gg$marker) <- rep(6,nrow(marker))
  if (!is.null(marker$PedigreeID)){
    iErr <- marker$PedCheck>0
    glyph_colour(gg$marker)[iErr] <- 5 + unlist(marker$PedCheck[iErr])
    glyph_type(gg$marker)[iErr] <- 7
  }
  OK <- FALSE
  while (!OK){
    cl.i <- 1
    shadowed(gg$marker) <- rep(FALSE,nrow(marker))
    while (cl.i<=nCl){
      str.i <- sprintf("Brush points belonging to cluster %d (of %d), OR
press <1> to skip cluster, OR
press <2> to skip remaining clusters, OR
press <3> to modify a specific cluster or set points to NA, OR
press <4> to start from scratch: ",cl.i,nCl)
      actn <- readline(prompt=str.i)
      if (actn=='1'){
        cl.i <- cl.i + 1
      }else if (actn=='2'){
        cl.i <- nCl + 1
      }else if (actn=='3'){
        OK.3 <- FALSE
        while (!OK.3){
          cl.i <- as.numeric(readline(str.3))
          OK.3 <- cl.i %in% c(1:nCl,8)
        } #end while OK.3
      }else if (actn=='4'){
        cl.i <- 1
        glc <- glyph_colour(gg$marker) <- rep(8,nrow(marker))
        marker$Call <- rep(NA,nrow(marker))
        if (!is.null(marker$PedCheck))
          marker$PedCheck <- rep(0,nrow(marker))
      }else{
        iCluster <- selected(gg$marker)
        marker$Call[iCluster] <- rep(cntIdeal[cl.i],sum(iCluster)) #(glc[iCluster]-1)/4
        glc[iCluster] <- rep(cntIdeal[cl.i]*4+1,sum(iCluster))
        assignedCl <- unique(glc[glc<=5])
        glyph_colour(gg$marker) <- glc
        if (!is.null(marker$PedigreeID) & length(assignedCl)>1){
          selM <- unique(substr(marker$PedigreeID[iCluster],1,4))
          selM <- selM[!substr(selM,2,4)%in%c('000','999')]
          redF <- paste(substr(marker$PedigreeID,1,1),substr(marker$PedigreeID,5,7),sep='')
          selF <- unique(redF[iCluster])
          selF <- selF[!substr(selF,2,4)%in%c('000','999')]
          indSel <- substr(marker$PedigreeID,1,4)%in%selM | redF%in%selF | marker$PedCheck>0
          marker$PedCheck[indSel] <- validateSingleCall(marker[indSel,],classification)
          iErr <- marker$PedCheck>0
          glyph_colour(gg$marker)[iErr] <- 5 + unlist(marker$PedCheck[iErr])
          glyph_type(gg$marker) <- glt
          glyph_type(gg$marker)[iErr] <- 7
        } #end if is.null
        cl.i <- cl.i + 1
      } #end if actn
    } #end while cl.i
    
    if (cl.i==8){
      actn <- readline(prompt=str.8)
      if (actn=='1'){
      }else{
        iErase <- selected(gg$marker)
        glc[iErase] <- rep(8,sum(iErase))
        marker$Call[iErase] <- rep(NA,sum(iErase))
        glyph_colour(gg$marker) <- glc
        if (!is.null(marker$PedCheck)){
          indSel <- marker$PedCheck>0
          marker$PedCheck[indSel] <- validateSingleCall(marker[indSel,],classification)
          iErr <- marker$PedCheck>0
          glyph_colour(gg$marker)[iErr] <- 5 + unlist(marker$PedCheck[iErr])
          glyph_type(gg$marker) <- glt
          glyph_type(gg$marker)[iErr] <- 7
        } #end if is.null
      } #end if actn
    }else{
      nErr <- sum(marker$PedCheck>0)
      if (nErr>0){
        shadowed(gg$marker) <- glyph_colour(gg$marker)%in%c(1:5,8)
        str.ped <- sprintf('Warning: %d ped-errors!
Brush points to un-assign, OR
press <1> to un-assign all erroneous offspring, OR
press <2> to modify clusters, OR
press <3> to disregard errors and finish: ',nErr)
        actn <- readline(prompt=str.ped)
        if (actn=='1'){
          iErase <- marker$PedCheck==1
          glc[iErase] <- rep(8,sum(iErase))
          marker$Call[iErase] <- rep(NA,sum(iErase))
          marker$PedCheck[marker$PedCheck>0] <- 0
          glyph_colour(gg$marker) <- glc
          glyph_type(gg$marker) <- glt
          #OK <- TRUE
        }else if (actn=='2'){
        }else if (actn=='3'){
          OK <- TRUE
        }else{
          iErase <- selected(gg$marker)
          glc[iErase] <- rep(8,sum(iErase))
          marker$Call[iErase] <- rep(NA,sum(iErase))
          glyph_colour(gg$marker) <- glc
          if (!is.null(marker$PedCheck)){
            indSel <- marker$PedCheck>0
            marker$PedCheck[indSel] <- validateSingleCall(marker[indSel,],classification)
            iErr <- marker$PedCheck>0
            glyph_colour(gg$marker)[iErr] <- 5 + unlist(marker$PedCheck[iErr])
            glyph_type(gg$marker) <- glt
            glyph_type(gg$marker)[iErr] <- 7
          } #end if is.null
        } #end if actn
      }else{
        finished <- readline(prompt=str.end)
        OK <- finished %in% c('y','Y','yes','Yes','YES','1','T','TRUE')
      } #end if nErr
    } #end if cl.i
  } #while !OK
  if (close.gg){
    close(gg)
    return(marker)
  }else{
    shadowed(gg$marker) <- rep(FALSE,nrow(marker))
    gg$marker$Call <- as.numeric(marker$Call)
    gg$marker$PedCheck <- marker$PedCheck
    return(gg)
  } #end if close.gg
}



#Validates a single SNP/marker, returns vector of pedErrors (o=1,p=2) 
validateSingleCall <- function(marker,classification){
  #indFamily <- unique(substr(marker$PedigreeID,2,7))
  #indFounding <- substr(indFamily,1,3)=='000' | substr(indFamily,4,6)=='000'
  indFamily <- unique(substr(marker$PedigreeID,1,7))
  indFamily <- indFamily[!indFamily%in%'']
  indFounding <- substr(indFamily,2,4)=='000' | substr(indFamily,5,7)=='000'
  indFamily <- indFamily[!indFounding]
  nFull <- length(indFamily)
  pedErrors <- rep(0,nrow(marker))
  if (classification%in%'MSV-a'){
    callVec <- marker$Call*2
  }else if (classification%in%'MSV-b'){
    callVec <- marker$Call*2-1
  }else{
    callVec <- marker$Call
  }
  for (j in 1:nFull){
    #jf <- marker$PedigreeID==paste('1000',substr(indFamily[j],4,6),'00',sep='')
    #jm <- marker$PedigreeID==paste('1',substr(indFamily[j],1,3),'00000',sep='')
    #j.fam <- substr(marker$PedigreeID,2,7)==indFamily[j] &! (jf|jm)
    jf <- marker$PedigreeID %in% paste(substr(indFamily[j],1,1),'000',substr(indFamily[j],5,7),'00',sep='')
    jm <- marker$PedigreeID %in% paste(substr(indFamily[j],1,4),'00000',sep='')
    j.fam <- substr(marker$PedigreeID,1,7) %in% indFamily[j] &! (jf|jm)
    bothP <- c(which(jf),which(jm))
    pCalls <- sort(callVec[bothP][!is.na(callVec[bothP])],index.return=TRUE)
    oCalls <- callVec[j.fam]
    if (length(pCalls$x)==1){
      iErr <- (oCalls<pCalls$x-.5 | oCalls>pCalls$x+.5) &! is.na(oCalls)
      if (any(iErr)){
        pedErrors[j.fam][iErr] <- 1
        pedErrors[bothP] <- 2
      } #if any
    }else if (length(pCalls$x)==2){
      if (pCalls$x[1]<=.5 & pCalls$x[2]>=.5){    #Alt. 1
        iHi <- oCalls>pCalls$x[1]+.5 &! is.na(oCalls)
        if (any(iHi)){
          pedErrors[j.fam][iHi] <- 1
          pedErrors[bothP[pCalls$ix[1]]] <- 2
        } #if any
        iLo <- oCalls<pCalls$x[2]-.5 &! is.na(oCalls)
        if (any(iLo)){
          pedErrors[j.fam][iLo] <- 1
          pedErrors[bothP[pCalls$ix[2]]] <- 2
        } #if any
      }else if(all(pCalls$x<.5)){    #Alt. 2
        iHi <- oCalls>sum(pCalls$x) &! is.na(oCalls)
        if (any(iHi)){
          pedErrors[j.fam][iHi] <- 1
          pedErrors[bothP] <- 2
        } #if any
      }else if(all(pCalls$x>.5)){      #Alt. 3
        iLo <- oCalls<sum(pCalls$x)-1 &! is.na(oCalls)
        if (any(iLo)){
          pedErrors[j.fam][iLo] <- 1
          pedErrors[bothP] <- 2
        } #if any
      } #if alternatives
    }else if (length(pCalls$x)>2){
      stop('More than two parents. Genotyped repeatedly, or ped-error?')
    } #end if length(pCalls)
  } #end for j
  pedErrors
}


  
#Set values for transformation/normalisation and return in a list. Suggests default values.
setNormOptions <- function(shearInf1=TRUE,transf='root',method='medianAF',minSize=suggestSh(shearInf1)$minSize,prob=suggestSh(shearInf1)$prob,nBins=suggestSh(shearInf1)$nBins,dist=suggestTr(transf)$dist,pNorm=suggestTr(transf)$pNorm,nthRoot=suggestTr(transf)$nthRoot,offset=suggestTr(transf)$offset,scale=suggestNo(method)$scale,nSD=3,breaks=200){
  suggestSh <- function(shearInf1){
    opt <- NULL
    if (shearInf1){
      opt$minSize <- 20; opt$prob <- .1; opt$nBins <- 50
    }else{
      opt$minSize <- 40; opt$prob <- .05; opt$nBins <- 50
    }
    opt
  }
  suggestTr <- function(transf){
    opt <- NULL
    if (transf=='none'){
      opt$dist <- 'manhattan'; opt$pNorm <- NULL; opt$nthRoot <- NULL; opt$offset <- NULL
    }else if (transf=='log'){
      opt$dist <- 'minkowski'; opt$pNorm <- 5; opt$nthRoot <- NULL; opt$offset <- 200
    }else if (transf=='root'){
      opt$dist <- 'minkowski'; opt$pNorm <- 5; opt$nthRoot <- 4; opt$offset <- 200
    }
    opt
  }
  suggestNo <- function(method){
    opt <- NULL
    if (method=='medianAF'){
      opt$scale <- NULL
    }else if (method=='linPeak'){
      opt$scale <- 0.75
    }
    opt
  }
  normOpts <- list(shearInf1=shearInf1,minSize=minSize,prob=prob,nBins=nBins,transf=transf,nthRoot=nthRoot,offset=offset,dist=dist,pNorm=pNorm,method=method,scale=scale,nSD=nSD,breaks=breaks)
  normOpts
}



#Set parameters for genotype-calling and return in a list. Suggests default values.
#Default values for largeSample==TRUE calibrated for more than 3000 objects,
#largeSample==FALSE calibrated for 100 objects.
#Calibrated for transf='root'; if transf='none', increase devCentLim and wSpreadLim
#Earlier versions: mseLim=1, probsQRange=c(.01,.5,.99)
setGenoOptions <- function(largeSample=FALSE,snpPerArrayLim=.8,arrayPerSnpLim=0,ploidy='tetra',filterLim=0,detectLim=.8,wSpreadLim=suggestGeno(largeSample)$wSpreadLim,devCentLim=.35,hwAlpha=suggestGeno(largeSample)$hwAlpha,probsIndSE=suggestGeno(largeSample)$probsIndSE,afList=seq(0,.5,.05),clAlpha=suggestGeno(largeSample)$clAlpha,rPenalty=2,rotationLim=suggestGeno(largeSample)$rotationLim,minClLim=5,nSdOverlap=2,minBin=suggestGeno(largeSample)$minBin,binWidth=suggestGeno(largeSample)$binWidth){
  suggestGeno <- function(largeSample){
    opt <- NULL
    if (largeSample){
      opt$probsIndSE=.9; opt$clAlpha=.01; opt$rotationLim=2; opt$wSpreadLim=.07; opt$hwAlpha=1e-10; opt$minBin=2; opt$binWidth=.05
    }else{
      opt$probsIndSE=NULL; opt$clAlpha=.05; opt$rotationLim=5; opt$wSpreadLim=.1; opt$hwAlpha=1e-4; opt$minBin=2; opt$binWidth=.1
    }
    opt
  }
  genoOpts <- list(snpPerArrayLim=snpPerArrayLim,arrayPerSnpLim=arrayPerSnpLim,ploidy=ploidy,filterLim=filterLim,detectLim=detectLim,wSpreadLim=wSpreadLim,devCentLim=devCentLim,hwAlpha=hwAlpha,probsIndSE=probsIndSE,afList=afList,clAlpha=clAlpha,rPenalty=rPenalty,rotationLim=rotationLim,minClLim=minClLim,nSdOverlap=nSdOverlap,minBin=minBin,binWidth=binWidth)
  genoOpts
}


#Set parameters for merging of genotypes and return in a list. Suggests default values.
setMergeOptions <- function(minC=NULL,noiseQuantile=.75,offspringLim=7,ratioLim=.9,rngLD=5){
  if (!is.null(minC))
    noiseQuantile <- NULL
  mergeOpts <- list(minC=minC,noiseQuantile=noiseQuantile,offspringLim=offspringLim,ratioLim=ratioLim,rngLD=rngLD)
  mergeOpts
}



#Takes MSV5-calls as input and partially resolves the true calls for each paralog.
#'Partially' because the actual chromosome number is not known, and the duplicate
#calls therefore have to be assigned assuming known location of paralogs of only
#one parent at the time. Where informative meiosis is found, each offspring is
#called as {Para1,Para2} <- {0,0}, {0,1}, or {1,1}. Can theoretically be used for
#LA by assuming only one known parent.
unmixParalogues <- function(BSRed,singleCalls=getSingleCalls(BSRed)){
  unmixTheta <- function(x){
    X <- rep(NA,2)
    if (x!=.5){
      X[1] <- round(x)
      X[2] <- 2*x - X[1]
      X <- sort(X)
    }
    X
  }
  lookupOffspring <- function(x,LookUp){
    X <- cbind(LookUp[x[-c(1,2)],x[1],x[2]],LookUp[x[-c(1,2)],x[1],x[2]])
    iHalf <- X[,1]%in%.5
    X[iHalf,1] <- 0
    X[iHalf,2] <- 1
    X <- as.vector(X)
    X
  }

  #Input control
  callType <- fData(BSRed)$Classification
  if (!is.null(fData(BSRed)$Manual.Calls.R)){
    iManual <- !fData(BSRed)$Manual.Calls.R %in% ''
    callType[iManual] <- fData(BSRed)$Manual.Calls.R[iManual]
  }
  if (!all(callType %in% 'MSV-5')){
    stop('This function accepts nothing but MSV-5 markers as input')
  }
  
  #Expand call-tables for father and mother to include both paralogues
  str1 <- paste(rownames(singleCalls),'Para1',sep='_')
  str2 <- paste(rownames(singleCalls),'Para2',sep='_')
  rNames <- as.vector(t(cbind(str1,str2)))
  nRows <- 2*nrow(BSRed)
  paraCalls <- NULL
  paraCalls$father <- paraCalls$mother <- matrix(nrow=nRows,ncol=ncol(BSRed),dimnames=list(rNames,sampleNames(BSRed)))
  paraCalls$mother[seq(1,nRows,2),] <- paraCalls$mother[seq(2,nRows,2),] <- paraCalls$father[seq(1,nRows,2),] <- paraCalls$father[seq(2,nRows,2),] <- singleCalls

  #Create look-up table for valid meioses (dim: 5(offspring) x 5(parent1) x 5(parent2))
  LookUp <- array(dim=rep(5,3))
  LookUp[c(1,6,11,26,27,31,36,51:53,56,61,77,78,82,87,103,108,113)] <- 0
  LookUp[c(7,17,33,42,59,67,84,93,109,119)] <- .5
  LookUp[c(13,18,23,39,44,48,49,65,70,73:75,90,95,99,100,115,120,125)] <- 1

  #Index including only <p><fff> (population + father-ID)
  pp <- paste(substr(pData(BSRed)$'PedigreeID',1,1),substr(pData(BSRed)$'PedigreeID',5,7),sep='')
  
  #Inheritance from mothers, fathers not fixed to paralogues yet
  indMothers <- unique(substr(pData(BSRed)$PedigreeID,1,4))
  indMothers <- indMothers[!substr(indMothers,2,4) %in% c('000','999')]
  for (j in 1:length(indMothers)){
    jm <- pData(BSRed)$PedigreeID %in% paste(indMothers[j],'00000',sep='')
    mCalls <- assayData(BSRed)$call[,jm]
    iNA <- is.na(singleCalls[,jm]) &! is.na(mCalls)
    if (sum(iNA)>0){
      paraCalls$mother[sort(c(which(iNA)*2,which(iNA)*2-1)),jm] <- as.vector(sapply(mCalls[iNA],unmixTheta))
      #Don't bother with homozygotes and unknown 0.5's:
      iCand.m <- !is.na(paraCalls$mother[seq(1,nRows,2),jm]) &! mCalls%in%c(0,1)
      j.fam <- substr(pData(BSRed)$PedigreeID,1,4) %in% indMothers[j] &! jm
      if (sum(j.fam)>2){   #Which minimum family-size to allow?
        indMates <- unique(pp[j.fam])
        indMates <- indMates[!substr(indMates,2,4) %in% '999']
        for (i in 1:length(indMates)){
          jf <- pData(BSRed)$PedigreeID %in% paste(substr(indMates[i],1,1),'000',substr(indMates[i],2,4),'00',sep='')
          fCalls <- assayData(BSRed)$call[,jf]
          i.fam <- pp[j.fam] %in% indMates[i]
          oCalls <- assayData(BSRed)$call[,j.fam][,i.fam]
          iCand.f <- fCalls!=.5 | !is.na(paraCalls$mother[seq(1,nRows,2),jf])
          iCand <- which(iCand.m & iCand.f)
          if (length(iCand)>0){
            indLookup <- matrix(cbind(mCalls,fCalls,oCalls)[iCand,]*4+1,nrow=length(iCand))
            cc <- apply(indLookup,1,lookupOffspring,LookUp)
            dim(cc) <- c(ncol(oCalls),length(iCand)*2)
            paraCalls$mother[sort(c(iCand*2-1,iCand*2)),j.fam][,i.fam] <- t(cc)
          }
        } #for i
      } #if sum(j.fam)
    } #if sum(iNA)
  } #for j
  
  #Inheritance from fathers, mothers not fixed to paralogues yet
  indFathers <- unique(pp)
  indFathers <- indFathers[!substr(indFathers,2,4) %in% c('000','999')]
  for (j in 1:length(indFathers)){
    jf <- pData(BSRed)$PedigreeID %in% paste(substr(indFathers[j],1,1),'000',substr(indFathers[j],2,4),'00',sep='')
    fCalls <- assayData(BSRed)$call[,jf]
    iNA <- is.na(singleCalls[,jf]) &! is.na(fCalls)
    if (sum(iNA)>0){
      paraCalls$father[sort(c(which(iNA)*2,which(iNA)*2-1)),jf] <- as.vector(sapply(fCalls[iNA],unmixTheta))
      #Don't bother with homozygotes and unknown 0.5's:
      iCand.f <- !is.na(paraCalls$father[seq(1,nRows,2),jf]) &! fCalls%in%c(0,1)
      j.fam <- pp %in% indFathers[j] &! jf
      if (sum(j.fam)>2){   #Which minimum family-size to allow?
        indMates <- unique(substr(pData(BSRed)$PedigreeID[j.fam],1,4))
        indMates <- indMates[!substr(indMates,2,4) %in% '999']
        for (i in 1:length(indMates)){
          jm <- pData(BSRed)$PedigreeID %in% paste(indMates[i],'00000',sep='')
          mCalls <- assayData(BSRed)$call[,jm]
          i.fam <- substr(pData(BSRed)$PedigreeID[j.fam],1,4) %in% indMates[i]
          oCalls <- assayData(BSRed)$call[,j.fam][,i.fam]
          iCand.m <- mCalls!=.5 | !is.na(paraCalls$father[seq(1,nRows,2),jm])
          iCand <- which(iCand.m & iCand.f)
          if (length(iCand)>0){
            indLookup <- matrix(cbind(fCalls,mCalls,oCalls)[iCand,]*4+1,nrow=length(iCand))
            cc <- apply(indLookup,1,lookupOffspring,LookUp)
            dim(cc) <- c(ncol(oCalls),length(iCand)*2)
            paraCalls$father[sort(c(iCand*2-1,iCand*2)),j.fam][,i.fam] <- t(cc)
          }
        } #for i
      } #if sum(j.fam)
    } #if sum(iNA)
  } #for j
  j.allFathers <- substr(pData(BSRed)$PedigreeID,2,4)%in%'000' & substr(pData(BSRed)$PedigreeID,8,9)%in%'00'
  paraCalls$mother[,j.allFathers] <- NA
  j.allMothers <- substr(pData(BSRed)$PedigreeID,5,9)%in%'00000'
  paraCalls$father[,j.allMothers] <- NA
  paraCalls
}



#Identify MSV-5's for which both paralogues are either AA, BB, or AB
getSingleCalls <- function(BSRed){
  compTheta <- function(x){
    rng <- c(max(0,x[1]-.5),min(1,x[1]+.5))
    indAB <- any(rng %in% x[-1]) &! is.na(x[1])
    indAB
  }
  singleCalls <- matrix(nrow=nrow(BSRed),ncol=ncol(BSRed),dimnames=list(featureNames(BSRed),sampleNames(BSRed)))
  singleCalls[assayData(BSRed)$call==0] <- 0
  singleCalls[assayData(BSRed)$call==1] <- 1
  indFamily <- unique(substr(pData(BSRed)$PedigreeID,1,7))
  indFounding <- substr(indFamily,2,4)=='000' | substr(indFamily,5,7)=='000'
  indFamily <- indFamily[!indFounding]
  nFull <- length(indFamily)
  for (j in 1:nFull){ #Parents with [AB,AB]-markers
    jf <- pData(BSRed)$PedigreeID %in% paste(substr(indFamily[j],1,1),'000',substr(indFamily[j],5,7),'00',sep='')
    jm <- pData(BSRed)$PedigreeID %in% paste(substr(indFamily[j],1,4),'00000',sep='')
    j.fam <- substr(pData(BSRed)$PedigreeID,1,7) %in% indFamily[j] &! (jf|jm)
    if (sum(jf)==1 & sum(jm)==1 & sum(j.fam)>=1){
      fCalls <- assayData(BSRed)$call[,jf]
      mCalls <- assayData(BSRed)$call[,jm]
      oCalls <- assayData(BSRed)$call[,j.fam]
      indF <- fCalls==.5 & !is.na(fCalls)
      if (any(indF)){
        indF[indF] <- apply(matrix(cbind(mCalls,oCalls)[indF,],nrow=sum(indF)),1,compTheta)
        singleCalls[indF,jf] <- .5
      }
      indM <- mCalls==.5 & !is.na(mCalls)
      if (any(indM)){
        indM[indM] <- apply(matrix(cbind(fCalls,oCalls)[indM,],nrow=sum(indM)),1,compTheta)
        singleCalls[indM,jm] <- .5
      }
    }
  }
  singleCalls
}



#Enables assignments to markers of 'AlleleSetIllumina', given phenoData are identical
assignToAlleleSet <- function(BSRed,BSAdd){
  if (!identical(phenoData(BSRed),phenoData(BSAdd)))
    stop('The phenoData of the AlleleSets differ. Will not merge!')
  indReplace <- which(featureNames(BSRed) %in% featureNames(BSAdd))
  orderRed <- order(featureNames(BSRed)[indReplace])
  indReplace <- indReplace[orderRed]
  orderAdd <- order(fData(BSAdd)$Name)
  if (!identical(featureData(BSRed),featureData(BSAdd))){
    varRed <- colnames(fData(BSRed))
    varAdd <- colnames(fData(BSAdd))
    if (!all(varAdd %in% varRed)){
      fData <- data.frame(stringsAsFactors=FALSE, row.names=featureNames(BSRed))
      varExtra <- varAdd[!varAdd %in% varRed]
      for (i in 1:length(varExtra)){
        fData[[varExtra[i]]] <- as('',class(fData(BSAdd)[[varExtra[i]]]))
      }
      fMetadata <- data.frame(labelDescription=featureData(BSAdd)@varMetadata[varExtra,], row.names=varExtra)
      featureData <- new('AnnotatedDataFrame', data=fData, varMetadata=fMetadata, dimLabels=c('featureNames','featureColumns'))
      featureData(BSRed) <- combine(featureData(BSRed), featureData)
      message(sprintf('%d new column(s) of featureData added to target',i))
      varRed <- colnames(fData(BSRed))
   }
    featureData(BSRed)@data[indReplace, varRed %in% varAdd] <- fData(BSAdd)[orderAdd, varRed[varRed %in% varAdd]]
    message('New featureData assigned to target')
  } #end if
  assayRed <- names(assayData(BSRed))
  assayAdd <- names(assayData(BSAdd))
  if (!all(assayAdd%in%assayRed)){
    assayExtra <- assayAdd[!assayAdd%in%assayRed]
    for (i in 1:length(assayExtra)){
      assayData(BSRed)[[assayExtra[i]]] <- matrix(nrow=nrow(BSRed), ncol=ncol(BSRed), dimnames=list(featureNames(BSRed),sampleNames(BSRed)))
      message(sprintf('Added assayData \'%s\' to target',assayExtra[i]))
    } #end for
  } #end if
  if (!all(assayRed%in%assayAdd))
    warning('One or more assayData-arrays in target not found in added AlleleSetIllumina-object')
  assayNames <- intersect(assayAdd, names(assayData(BSRed)))
  nAssays <- length(assayNames)
  for (i in 1:nAssays){
    Xred <- assayData(BSRed)[[assayNames[i]]][indReplace,]
    Xadd <- assayData(BSAdd)[[assayNames[i]]][orderAdd,]
    if (!identical(Xred,Xadd)){
      assayData(BSRed)[[assayNames[i]]][indReplace,] <- Xadd
      message(sprintf('New %s-data assigned to target',assayNames[i]))
    } #end if
  } #end for
  validObject(BSRed)
  BSRed
}



#Translates calls (in {0,1/2,1}) into true genotypes (A,T,C,G). 'type' is either
#'regular', 'single' (markernames w. 'Para'), or 'merged' (markernames w. 'Chrom')
#Output is a matrix of the same size as 'calls', each element given as "x y",
#where x and y are the two alleles. Missing: "- -"
#If type is 'regular', and non-diploid markers are found, these are attempted
#translated into a diploid representation where meaningful.
translateTheta <- function(calls,resInfo,type='regular'){
  subTrans <- function(x,alleles){
    tG <- character(length(x))
    tG[x%in%0] <- paste(alleles[,1],alleles[,1])[x%in%0]
    tG[x%in%0.5] <- paste(alleles[,1],alleles[,2])[x%in%0.5]
    tG[x%in%1] <- paste(alleles[,2],alleles[,2])[x%in%1]
    tG
  }
  if (is.vector(calls))
    calls <- as.matrix(calls)
  if (type%in%'regular'){
    markerNames <- rownames(calls)
    fData <- data.frame(Classification=resInfo[markerNames,'Classification'], row.names=markerNames)
    #if (any(c(.25,.75) %in% unique(as.vector(calls)))){
    if (any(!unique(fData$Classification) %in% c('MONO-a','MONO-b','SNP','FAIL','MONO-filt'))){
      message('Non-diploid calls found. Will be converted to a diploid representation where possible')
      calls <- makeDiploidCalls(calls,fData)
    }
  }else if (type%in%'single'){
    markerNames <- sub('_Para[12]','',rownames(calls))
  }else if (type%in%'merged'){
    markerNames <- sub('_Chrom[0-9XY]+','',rownames(calls))
  }
  cc <- data.frame(Base=c('A','T','C','G'),Compl=c('t','a','g','c'))
  snpInfo <- resInfo[markerNames,c('SNP','ILMN.Strand')]
  topSNP <- snpInfo[,1]
  for (i in 1:4)
    topSNP[snpInfo[,2]%in%'BOT'] <- sub(cc[i,1],cc[i,2],topSNP[snpInfo[,2]%in%'BOT'])
  alleles <- toupper(do.call('rbind',strsplit(substr(topSNP,2,4),'/')))
  genotype <- apply(calls,2,subTrans,alleles)
  dimnames(genotype) <- dimnames(calls)
  genotype[is.na(calls)] <- '- -'
  genotype
}



#Will convert non-diploid calls to a diploid representation where possible (i.e. for MSV-a and -b's)
makeDiploidCalls <- function(calls,fData){
  diploidCalls <- calls
  callType <- fData$Classification
  if (!is.null(fData$Manual.Calls.R)){
    iManual <- !fData$Manual.Calls.R %in% ''
    callType[iManual] <- fData$Manual.Calls.R[iManual]
  }
  diploidCalls[callType%in%'MSV-a'] <- diploidCalls[callType%in%'MSV-a']*2
  diploidCalls[callType%in%'MSV-b'] <- diploidCalls[callType%in%'MSV-b']*2-1
  diploidCalls[callType%in%'MSV-5'] <- NA
  diploidCalls[callType%in%'PSV'] <- NA
  diploidCalls
}



#Translates calls (in {0, 1/4, 1/2, 3/4, 1}) into true genotypes (A,T,C,G). Genotypes
#for duplicated markers given as 4 letters, 2 letters for diploid loci. Resolved
#paralogs separated by comma.
translateThetaCombined <- function(BSRed,mergedCalls=NULL){
  subTransDiploid <- function(x,alleles){
    tG <- character(length(x))
    tG[x%in%0] <- paste(alleles[,1],alleles[,1],sep='')[x%in%0]
    tG[x%in%0.5] <- paste(alleles[,1],alleles[,2],sep='')[x%in%0.5]
    tG[x%in%1] <- paste(alleles[,2],alleles[,2],sep='')[x%in%1]
    tG[tG %in% ""] <- '--'
    tG
  }
  subTransMSVa <- function(x,alleles){
    tG <- character(length(x))
    tG[x%in%0] <- paste(alleles[,1],alleles[,1],',',alleles[,1],alleles[,1],sep='')[x%in%0]
    tG[x%in%0.25] <- paste(alleles[,1],alleles[,2],',',alleles[,1],alleles[,1],sep='')[x%in%0.25]
    tG[x%in%.5] <- paste(alleles[,2],alleles[,2],',',alleles[,1],alleles[,1],sep='')[x%in%.5]
    tG[tG %in% ""] <- '--,--'
    tG
  }
  subTransMSVb <- function(x,alleles){
    tG <- character(length(x))
    tG[x%in%0.5] <- paste(alleles[,1],alleles[,1],',',alleles[,2],alleles[,2],sep='')[x%in%0.5]
    tG[x%in%0.75] <- paste(alleles[,1],alleles[,2],',',alleles[,2],alleles[,2],sep='')[x%in%0.75]
    tG[x%in%1] <- paste(alleles[,2],alleles[,2],',',alleles[,2],alleles[,2],sep='')[x%in%1]
    tG[tG %in% ""] <- '--,--'
    tG
  }
  subTransMSV5 <- function(x,alleles){
    tG <- character(length(x))
    tG[x%in%0] <- paste(alleles[,1],alleles[,1],',',alleles[,1],alleles[,1],sep='')[x%in%0]
    tG[x%in%0.25] <- paste(alleles[,1],alleles[,1],alleles[,1],alleles[,2],sep='')[x%in%0.25]
    tG[x%in%0.5] <- paste(alleles[,1],alleles[,1],alleles[,2],alleles[,2],sep='')[x%in%0.5]
    tG[x%in%0.75] =paste(alleles[,1],alleles[,2],alleles[,2],alleles[,2],sep='')[x%in%0.75]
    tG[x%in%1] <- paste(alleles[,2],alleles[,2],',',alleles[,2],alleles[,2],sep='')[x%in%1]
    tG[tG %in% ""] <- '--,--'
    tG
  }
  if (!is.null(mergedCalls)){
    if (!identical(colnames(mergedCalls),sampleNames(BSRed))){
      stop('"mergedCalls" must be either NULL or a matrix with sample-info in the columns')
    }
    mergedNames <- sub('_Chrom[0-9XY]+','',rownames(mergedCalls))
  }
  
  call <- assayData(BSRed)$call
  callType <- fData(BSRed)$Classification
  if (!is.null(fData(BSRed)$Manual.Calls.R)){
    iManual <- !fData(BSRed)$Manual.Calls.R%in%''
    callType[iManual] <- fData(BSRed)$Manual.Calls.R[iManual]
  }
  cc <- data.frame(Base=c('A','T','C','G'),Compl=c('t','a','g','c'))
  topSNP <- fData(BSRed)$SNP
  for (i in 1:4)
    topSNP[fData(BSRed)$ILMN.Strand%in%'BOT'] <- sub(cc[i,1],cc[i,2],topSNP[fData(BSRed)$ILMN.Strand%in%'BOT'])
  alleles <- toupper(do.call('rbind',strsplit(substr(topSNP,2,4),'/')))
  rownames(alleles) <- featureNames(BSRed)
  genotype <- matrix('',nrow=nrow(call),ncol=ncol(call),dimnames=dimnames(call))
  iDiploid <- callType %in% c('MONO-a','MONO-b','SNP')
  call.i <- matrix(call[iDiploid,],nrow=sum(iDiploid),ncol=ncol(call))
  alleles.i <- matrix(alleles[iDiploid,],nrow=sum(iDiploid),ncol=2)
  genotype[iDiploid,] <- apply(call.i,2,subTransDiploid,alleles.i)
  iPSV <- callType %in% 'PSV'
  genotype[iPSV,] <- paste(alleles[iPSV,1],alleles[iPSV,1],',',alleles[iPSV,2],alleles[iPSV,2],sep='')
  iMSVa <- callType %in% 'MSV-a'
  call.i <- matrix(call[iMSVa,],nrow=sum(iMSVa),ncol=ncol(call))
  alleles.i <- matrix(alleles[iMSVa,],nrow=sum(iMSVa),ncol=2)
  genotype[iMSVa,] <- apply(call.i,2,subTransMSVa,alleles.i)
  iMSVb <- callType %in% 'MSV-b'
  call.i <- matrix(call[iMSVb,],nrow=sum(iMSVb),ncol=ncol(call))
  alleles.i <- matrix(alleles[iMSVb,],nrow=sum(iMSVb),ncol=2)
  genotype[iMSVb,] <- apply(call.i,2,subTransMSVb,alleles.i)
  iMSV5 <- callType %in% 'MSV-5'
  call.i <- matrix(call[iMSV5,],nrow=sum(iMSV5),ncol=ncol(call))
  alleles.i <- matrix(alleles[iMSV5,],nrow=sum(iMSV5),ncol=2)
  genotype[iMSV5,] <- apply(call.i,2,subTransMSV5,alleles.i)
  iFAIL <- callType %in% 'FAIL'
  genotype[iFAIL,] <- '--'
  if (!is.null(mergedCalls)){
    alleles.i <- matrix(alleles[mergedNames,],nrow=nrow(mergedCalls),ncol=2)
    trueSplit <- apply(mergedCalls,2,subTransDiploid,alleles.i)
    rownames(trueSplit) <- mergedNames
    trueMerged <- genotype[iMSV5,]
    for (i in 1:sum(iMSV5)){
      #i=i+1
      marker.i <- rownames(trueMerged)[i]
      ii <- which(rownames(trueSplit)%in%marker.i)
      split.i <- matrix(trueSplit[ii,],nrow=length(ii),ncol=ncol(trueSplit))
      #tt <- rbind(mergedCalls=mergedCalls[ii,],trueSplit=trueSplit[ii,],call=call[marker.i,],genotype=trueMerged[marker.i,],oldGeno=oldGeno[marker.i,])
      #print(tt[,101:110])
      if (nrow(split.i)==2){
        iRes <- apply(split.i,2,function(x) !all(x %in% '--'))
        trueMerged[i,iRes] <- apply(split.i[,iRes],2,function(x) paste(x[1],x[2],sep=','))
      }else if (nrow(split.i)==1){
        iRes <- !(split.i%in%'--' | call[marker.i,]%in%c(0,1))
        split.i2 <- sub(split.i[iRes],'',genotype[marker.i,iRes])
        trueMerged[i,iRes] <- paste(split.i[iRes],split.i2,sep=',')
      }
    }
    genotype[iMSV5,] <- trueMerged
  }
  
  assayData(BSRed)$genotype <- genotype
  validObject(BSRed)
  BSRed
}



#Reads call-data from file and writes genotype to a file of similar name
translateThetaFromFiles <- function(dataFiles,mergedCalls=NULL,markerStep=1000,sep='\t',quote=''){
  message("Loading featureData...")
  beadInfo <- read.table(dataFiles[['resFile']],header=TRUE,sep=sep,quote=quote,as.is=TRUE)
  callHeader <- c('',unlist(read.table(dataFiles[['callFile']],nrows=1,as.is=TRUE,sep=sep)))
  callCols <- rep('numeric',length(callHeader))
  callCols[1] <- 'character'
  file <- sub('Call','Alleles',dataFiles[['callFile']])
  for (i in seq(1,nrow(beadInfo),markerStep)){
    minM <- min(nrow(beadInfo),i+markerStep-1)
    message(paste('Translating markers ',i,':',minM,' (of ',nrow(beadInfo),') ...',sep=''))
    call <- read.table(dataFiles[['callFile']],header=FALSE,sep=sep,quote=quote,nrows=minM+1-i,colClasses=callCols,skip=i,row.names=1,col.names=callHeader)
    BSRed <- new('MultiSet',call=as.matrix(call),featureData=as(beadInfo[i:minM,],'AnnotatedDataFrame'),storage.mode='list')
    BSRed <- translateThetaCombined(BSRed,mergedCalls)
    write.table(assayData(BSRed)$genotype,file=file,append=i>1,sep=sep,col.names=i==1)
  }
}



#Resolves the within half-sib family inheritance from heterozygous parents. NB! No MSV's accepted
resolveInheritanceSNP <- function(BSSnp){
  uCalls <- unique(as.vector(assayData(BSSnp)$call))
  if (any(!uCalls %in% c(0,.5,1,NA))){
    stop('Only SNPs, no MSVs, are allowed. Calls found with illegal values')
  }
  inheritance <- NULL
  inheritance$mother <- inheritance$father <- matrix(nrow=nrow(BSSnp),ncol=ncol(BSSnp),dimnames=list(featureNames(BSSnp),sampleNames(BSSnp)))
  i0 <- assayData(BSSnp)$call==0
  inheritance$mother[i0] <- inheritance$father[i0] <- 0
  i1 <- assayData(BSSnp)$call==1
  inheritance$mother[i1] <- inheritance$father[i1] <- 1
  #LookUp <- matrix(c(0,0,NA,1,NA,0,NA,1,1),nrow=3,byrow=T) # Offspring x P2 (P1 is AB)

  #Index including only <p><fff> (population + father-ID)
  pp <- paste(substr(pData(BSSnp)$'PedigreeID',1,1),substr(pData(BSSnp)$'PedigreeID',5,7),sep='')
  
  #Inheritance from mothers
  indMothers <- unique(substr(pData(BSSnp)$PedigreeID,1,4))
  indMothers <- indMothers[!substr(indMothers,2,4) %in% c('000','999')]
  jm <- pData(BSSnp)$PedigreeID %in% paste(indMothers,'00000',sep='')
  inheritance$mother[,jm] <- inheritance$father[,jm] <- assayData(BSSnp)$call[,jm]
  for (j in 1:length(indMothers)){
    jm <- pData(BSSnp)$PedigreeID %in% paste(indMothers[j],'00000',sep='')
    mCalls <- assayData(BSSnp)$call[,jm]
    iCand.m <- mCalls %in% .5
    j.fam <- substr(pData(BSSnp)$PedigreeID,1,4) %in% indMothers[j] &! jm
    if (sum(j.fam)>2){   #Which minimum family-size to allow?
      indMates <- unique(pp[j.fam])
      indMates <- indMates[!substr(indMates,2,4) %in% '999']
      for (i in 1:length(indMates)){
        jf <- pData(BSSnp)$PedigreeID %in% paste(substr(indMates[i],1,1),'000',substr(indMates[i],2,4),'00',sep='')
        fCalls <- assayData(BSSnp)$call[iCand.m,jf]
        iCand.f <- fCalls!=.5  & !is.na(fCalls)
        i.fam <- pp[j.fam] %in% indMates[i]
        oCalls <- assayData(BSSnp)$call[iCand.m,j.fam][iCand.f,i.fam]
        iCand.o <- oCalls %in% .5
        cc <- 1 - fCalls[iCand.f]%*%t(rep(1,sum(i.fam)))
        inheritance$mother[iCand.m,j.fam][iCand.f,i.fam][iCand.o] <- cc[iCand.o]
        } #for i
    } #if sum
  } #for j
  
  #Inheritance from fathers
  indFathers <- unique(pp)
  indFathers <- indFathers[!substr(indFathers,2,4) %in% c('000','999')]
  jf <- pData(BSSnp)$PedigreeID %in% paste(substr(indFathers,1,1),'000',substr(indFathers,2,4),'00',sep='')
  inheritance$mother[,jf] <- inheritance$father[,jf] <- assayData(BSSnp)$call[,jf]
  for (j in 1:length(indFathers)){
    jf <- pData(BSSnp)$PedigreeID %in% paste(substr(indFathers[j],1,1),'000',substr(indFathers[j],2,4),'00',sep='')
    fCalls <- assayData(BSSnp)$call[,jf]
    iCand.f <- fCalls %in% .5
    j.fam <- pp %in% indFathers[j] &! jf
    if (sum(j.fam)>2){   #Which minimum family-size to allow?
      indMates <- unique(substr(pData(BSSnp)$PedigreeID[j.fam],1,4))
      indMates <- indMates[!substr(indMates,2,4) %in% '999']
      for (i in 1:length(indMates)){
        jm <- pData(BSSnp)$PedigreeID %in% paste(indMates[i],'00000',sep='')
        mCalls <- assayData(BSSnp)$call[iCand.f,jm]
        iCand.m <- mCalls!=.5  & !is.na(mCalls)
        i.fam <- substr(pData(BSSnp)$PedigreeID[j.fam],1,4) %in% indMates[i]
        oCalls <- assayData(BSSnp)$call[iCand.f,j.fam][iCand.m,i.fam]
        iCand.o <- oCalls %in% .5
        cc <- 1 - mCalls[iCand.m]%*%t(rep(1,sum(i.fam)))
        inheritance$father[iCand.f,j.fam][iCand.m,i.fam][iCand.o] <- cc[iCand.o]
       } #for i
    } #if sum
  } #for j
  inheritance
}



#Compares each split MSV-5 with known markers within each father and mother half-sib family
#and counts the IBD-occurences. Returns a matrix of dim(#MSV-5, #Chromosomes, 2 parents)
#from which the peaks can be used to assign each homeolog to its respective chromosome.
locateParalogues <- function(BSSnp,paraCalls,inheritP,offspringLim=7,ratioLim=.9){
  indMothers <- unique(substr(pData(BSSnp)$PedigreeID,1,4))
  indMothers <- indMothers[!substr(indMothers,2,4) %in% c('000','999')]
  nMothers <- length(indMothers)
  pp <- paste(substr(pData(BSSnp)$'PedigreeID',1,1),substr(pData(BSSnp)$'PedigreeID',5,7),sep='')
  indFathers <- unique(pp)
  indFathers <- indFathers[!substr(indFathers,2,4) %in% c('000','999')]
  nFathers <- length(indFathers)
  nPara <- nrow(paraCalls$m)
  countsGood <- NULL
  countsGood$M <- countsGood$F <- matrix(0,nrow=nPara/2,ncol=nrow(BSSnp),dimnames=list(unique(sub('_Para[12]','',rownames(paraCalls$m))),fData(BSSnp)$Name))
  for (i in 1:nPara){    #One paralogue at the time
    if (i %in% seq(1,nPara,2))
      message(sprintf('Mapping paralogue %d/%d: %s',i,nPara,rownames(paraCalls$m)[i]))
    nDen <- nNum <- matrix(nrow=nMothers,ncol=nrow(BSSnp))
    for (j in 1:nMothers){    #Inheritance from mothers
      jm <- pData(BSSnp)$PedigreeID %in% paste(indMothers[j],'00000',sep='')
      if (paraCalls$m[i,jm] %in% .5){    #Mother j heterozygous in paralogue i
        iHetSnp <- assayData(BSSnp)$call[,jm] %in% .5   #Polymorphic markers
        j.fam <- substr(pData(BSSnp)$PedigreeID,1,4) %in% indMothers[j] &! jm
        iSamp <- !is.na(paraCalls$m[i,j.fam])
        if (sum(iSamp)>1){   #Minimum 2 offspring w. known genotype
          snps <- as.matrix(inheritP$m[iHetSnp,j.fam][,iSamp],nrow=sum(iHetSnp))
          nDen[j,iHetSnp] <- apply(snps,1,function(x) sum(!is.na(x)))
          nNum[j,iHetSnp] <- apply(snps,1,function(x,p) max(sum(x==p,na.rm=TRUE),sum(x==1-p,na.rm=TRUE)),paraCalls$m[i,j.fam][iSamp])
        } #if sum
      } #if polymorphic
    } #for j
    indGood <- matrix(nNum>=offspringLim & nNum/nDen>=ratioLim,nrow=nrow(nNum))
    cGood <- apply(indGood,2,function(x) sum(x[!is.na(x)]))
    ii <- sub('_Para[12]','',rownames(paraCalls$m)[i])
    countsGood$M[ii,] <- countsGood$M[ii,] + cGood
      
    nDen <- nNum <- matrix(nrow=nFathers,ncol=nrow(BSSnp))
    for (j in 1:nFathers){    #Inheritance from fathers
      jf <- pData(BSSnp)$PedigreeID %in% paste(substr(indFathers[j],1,1),'000',substr(indFathers[j],2,4),'00',sep='')
      if (paraCalls$f[i,jf] %in% .5){    #Father j heterozygous in paralogue i
        iHetSnp <- assayData(BSSnp)$call[,jf] %in% .5   #Polymorphic markers
        j.fam <- pp %in% indFathers[j] &! jf
        iSamp <- !is.na(paraCalls$f[i,j.fam])
        if (sum(iSamp)>1){   #Minimum 2 offspring w. known genotype
          snps <- as.matrix(inheritP$f[iHetSnp,j.fam][,iSamp],nrow=sum(iHetSnp))
          nDen[j,iHetSnp] <- apply(snps,1,function(x) sum(!is.na(x)))
          nNum[j,iHetSnp] <- apply(snps,1,function(x,p) max(sum(x==p,na.rm=TRUE),sum(x==1-p,na.rm=TRUE)),paraCalls$f[i,j.fam][iSamp])
        } #if sum
      } #if polymorphic
    } #for j
    indGood <- matrix(nNum>=offspringLim & nNum/nDen>=ratioLim,nrow=nrow(nNum))
    cGood <- apply(indGood,2,function(x) sum(x[!is.na(x)]))
    ii <- sub('_Para[12]','',rownames(paraCalls$f)[i])
    countsGood$F[ii,] <- countsGood$F[ii,] + cGood
  } #for i
  #chromNames <- sort(unique(lMap$Ssa))
  chromNames <- sort(unique(fData(BSSnp)$Chr.Name))
  nChroms <- length(chromNames)
  chromHits <- array(0,dim=c(nPara/2,nChroms,2),dimnames=list(rownames(countsGood$F),chromNames,c('mother','father')))
  for (i in 1:(nPara/2)){
    #rr <- data.frame(Chr.Name=fData(BSSnp)$Chr.Name, mCounts=countsGood$M[i,], fCounts=countsGood$F[i,])
    #chromHits[i,,1] <- tapply(rr$mCounts,rr$Chr.Name,sum)
    #chromHits[i,,2] <- tapply(rr$fCounts,rr$Chr.Name,sum)
    #rr <- data.frame(mCounts=countsGood$M[i,], fCounts=countsGood$F[i,])
    chromHits[i,,1] <- tapply(countsGood$M[i,],fData(BSSnp)$Chr.Name,sum)
    chromHits[i,,2] <- tapply(countsGood$F[i,],fData(BSSnp)$Chr.Name,sum)
  } #for
  lenChrom <- tapply(rep(1,nrow(BSSnp)),fData(BSSnp)$Chr.Name,sum)
  cHitsPerMarker <- chromHits/array(rep(rep(1,nPara/2)%*%t(lenChrom),2),dim=dim(chromHits))
  nCountsTot <- cbind(mother=apply(countsGood$M,1,sum), father=apply(countsGood$F,1,sum))
  
  testP <- FALSE
  if (testP){
    plotCountsChrom(cHitsPerMarker,1:min(16,(nPara/2)))
    countG <- countsGood$F[1,]
    rr <- cbind(fData(BSSnp)[,c('Chr.Name','Chr.Index')],countG)
    rr <- rr[order(rr$Chr.Index),]
    rr <- rr[order(rr$Chr.Name),]
    print(rr[order(rr$countG),][-(1:(nrow(BSSnp)-30)),])
  }
  chromHits <- NULL
  chromHits$cPerMarker <- cHitsPerMarker
  chromHits$nCountsTot <- nCountsTot
  chromHits
}



#For each marker, plot the number of hits per chromosome for each chromosome
plotCountsChrom <- function(chromHits,markers=1:16,...){
  nSnp <- length(markers)
  pp <- par()
  nC <- ceiling(sqrt(nSnp))
  nR <- ceiling(nSnp/nC)
  par(mfrow=c(nR,nC),mai=c(1,1,1,.5)/sqrt(nC))
  xlim <- c(1,dim(chromHits)[2])
  for (i in 1:nSnp){
    plot.new()
    ylim <- c(0,max(chromHits[markers[i],,]))
    plot.window(xlim,ylim); axis(1,...); axis(2)
    lines(xlim[1]:xlim[2],chromHits[markers[i],,1],col='black')
    lines(xlim[1]:xlim[2],chromHits[markers[i],,2],col='red')
    title(main=dimnames(chromHits)[[1]][markers[i]],xlab='Chrom. #',ylab='Counts')
  }
  par(pp[c('mfrow','mai')])
}



#Assign paralogs to specific chromosomes
assignParalogues <- function(BSSnp,BSRed,paraCalls=unmixParalogues(BSRed,singleCalls),inheritP=resolveInheritanceSNP(BSSnp),singleCalls=getSingleCalls(BSRed),cHits=locateParalogues(BSSnp,paraCalls,inheritP,mO$offspringLim,mO$ratioLim)$cPerMarker,mO=setMergeOptions()){
  indMothers <- unique(substr(pData(BSSnp)$PedigreeID,1,4))
  indMothers <- indMothers[!substr(indMothers,2,4) %in% c('000','999')]
  nMothers <- length(indMothers)
  pp <- paste(substr(pData(BSSnp)$'PedigreeID',1,1),substr(pData(BSSnp)$'PedigreeID',5,7),sep='')
  indFathers <- unique(pp)
  indFathers <- indFathers[!substr(indFathers,2,4) %in% c('000','999')]
  nFathers <- length(indFathers)

  #Select markers that can be assigned to chromosomes and name them accordingly
  dimCH =dim(cHits)
  bestCH <- t(apply(cHits,1,function(x) apply(x,1,max)))    #Use largest value betw. parents
  orderC <- apply(bestCH,1,function(x) order(x,decreasing=TRUE))[1:3,]
  if (is.null(mO$minC)){
    noiseLevels <- apply(cbind(bestCH,orderC[3,]),1,function(x,n) x[x[n]],n=dimCH[2]+1)
    mO$minC <- quantile(noiseLevels, probs=mO$noiseQuantile)   #Quantile of 3rd largest peak across markers
    str <- 'The %.0fth percentile of the 3rd largest chromosome peaks gives minC = %.3f'
    message(sprintf(str, mO$noiseQuantile*100, mO$minC))
  }
  maxCounts <- apply(cHits,1,max)
  iOK <- which(maxCounts>mO$minC)
  iRunning <- sort(c(2*iOK-1,2*iOK))
  rNames <- NULL
  step <- 1
  altC <- matrix(nrow=dimCH[1],ncol=2,dimnames=list(names(maxCounts),c('P1','P2')))
  cNumbers <- sub('[a-z0]+','',colnames(cHits))
  for (i in iOK){
    altInd <- orderC[,i]
    ii <- bestCH[i,altInd]>mO$minC & c(TRUE,TRUE,FALSE)   #More than two matches not allowed
    altCH <- bestCH[i,altInd]
    if (ii[2]){
      ii[2] <- altCH[2]>2*altCH[3]   #Requires twice as many matches as max(noise)
    }
    altCH <- altCH[ii]
    #print(i);print(altInd);print(altCH)
    if (length(altCH)==2){
      sAltInd <- sort(altInd[1:2])
      rNames[step:(step+1)] <- paste(colnames(orderC)[i],'_Chrom',cNumbers[sAltInd],sep='')
      altC[i,] <- cNumbers[sAltInd]
      step <- step + 2
    }else{
      rNames[step] <- paste(colnames(orderC)[i],'_Chrom',cNumbers[altInd[1]],sep='')
      altC[i,1] <- cNumbers[altInd[1]]
      iAlt <- which(names(iRunning)%in%colnames(orderC)[i])
      iRunning <- iRunning[-iAlt[2]]
      step <- step + 1
    } #if
  } #for
  nNew <- length(rNames)
  nPara <- nrow(paraCalls$m)
  markerNames <- sub('_Chrom[0-9XY]+','',rNames)
  mergedCalls <- singleCalls[markerNames,]
  filledIn <- !is.na(singleCalls)
  rownames(mergedCalls) <- rNames
  indSCalls <- NULL    #each element is the corresponding element-index in paraCalls
  colOffset <- rep(1,nNew) %*% t(seq(0,ncol(mergedCalls)*nPara-1,nPara))
  indSCalls$running <- indSCalls$m <- indSCalls$f <- matrix(iRunning,nrow=nNew,ncol=ncol(paraCalls$m),dimnames=dimnames(mergedCalls)) + colOffset
  indSCalls$m[!filledIn[markerNames,]] <- indSCalls$f[!filledIn[markerNames,]] <- NA
  #partnerCandidates <- list()   #Keeps ranking of chromosomes where both paralogues not known
  positionFemale <- matrix(nrow=nNew,ncol=nMothers,dimnames=list(rNames,indMothers))
  positionMale <- matrix(nrow=nNew,ncol=nFathers,dimnames=list(rNames,indFathers))
  
  #Fill in first for mother, then for father
  str1 <- 'Inheritance each parent, paralogue %d/%d: %s'
  for (i in 1:nPara){    #One paralogue at the time
    if (i %in% seq(1,nPara,10))
      message(sprintf(str1,i,nPara,rownames(paraCalls$m)[i]))
    for (j in 1:nMothers){    #Inheritance from mothers
      jm <- pData(BSSnp)$PedigreeID %in% paste(indMothers[j],'00000',sep='')#; paraCalls$m[i,jm] %in% .5
      if (paraCalls$m[i,jm] %in% .5){    #Mother j heterozygous in paralogue i
        iHetSnp <- assayData(BSSnp)$call[,jm] %in% .5   #Polymorphic markers
        j.fam <- substr(pData(BSSnp)$PedigreeID,1,4) %in% indMothers[j] &! jm
        iSamp <- !is.na(paraCalls$m[i,j.fam])
        iName <- strsplit(rownames(paraCalls$m)[i],'_Para')[[1]]
        nHits <- sum(!is.na(altC[iName[1],]))
        #sum(iSamp);nHits;filledIn[iName[1],jm]
        if (sum(iSamp)>=mO$offspringLim & nHits>0 &! filledIn[iName[1],jm]){
          snps <- as.matrix(inheritP$m[iHetSnp,j.fam][,iSamp],nrow=sum(iHetSnp))
          nDen <- apply(snps,1,function(x) sum(!is.na(x)))
          nNum <- apply(snps,1,function(x,p) max(sum(x==p,na.rm=TRUE),sum(x==1-p,na.rm=TRUE)),paraCalls$m[i,j.fam][iSamp])
          #rr <- cbind(lMap[iHetSnp,],nNum,nDen,Ratio=nNum/nDen)
          rr <- cbind(fData(BSSnp)[iHetSnp,c('Chromosome','Female')],nNum,nDen,Ratio=nNum/nDen)
          if (nHits==2){
            bestNum <- totalDen <- medianRF <- linkRatio <- rep(0,nHits)
            for (k in 1:nHits){   #Test both alternatives
              #indC <- rr$Ssa %in% paste('ssa',sprintf('%.02d',altC[iName[1],k]),sep='')#; rr[indC,]
              indC <- rr$Chromosome %in% altC[iName[1],k]#; rr[indC,]
              indGood <- rr[indC,'Ratio']>=mO$ratioLim & rr[indC,'nNum']>=mO$offspringLim#; rr[indC,][indGood,]
              while (any(indGood)){
                medGPos <- quantile(rr[indC,'Female'][indGood],na.rm=TRUE,probs=.5,type=1)#;medGPos #No averaging!
                indCand <- rr[indC,'Female']>medGPos-mO$rngLD & rr[indC,'Female']<medGPos+mO$rngLD#; rr[indC,][indCand,]
                lRat <- sum(rr[indC,'nNum'][indCand])/sum(rr[indC,'nDen'][indCand])#;lRat
                bNum <- max(rr[indC,'nNum'][indCand][rr[indC,'Ratio'][indCand]>=mO$ratioLim],na.rm=TRUE)#;bNum
                #snps[indC,][indCand,]; assayData(BSSnp)$call[iHetSnp,j.fam][indC,iSamp][indCand,]
                tDen <- sum(apply(!is.na(matrix(matrix(snps[indC,],nrow=sum(indC))[indCand,],nrow=sum(indCand))),2,any))#;tDen
                if (lRat>linkRatio[k] | (lRat==linkRatio[k] & tDen>totalDen[k])){
                  linkRatio[k] <- lRat
                  medianRF[k] <- medGPos
                  totalDen[k] <- tDen
                  bestNum[k] <- bNum
                } #if
                indGood[indGood]=!indCand[indGood]#; rr[indC,][indGood,]
              } #while
            } #for k
            #totalDen;bestNum;medianRF;linkRatio
            if (any(linkRatio>=mO$ratioLim & totalDen>=mO$offspringLim) &! (diff(linkRatio)==0 & diff(totalDen)==0)){
              rankFit <- order(linkRatio,totalDen,decreasing=TRUE)
              iMerge <- paste(iName[1],'_Chrom',altC[iName[1],rankFit],sep='')#;iMerge
              iSingle1 <- which(names(iRunning) %in% iName[1])#;iSingle1
              iOrder1 <- c(iSingle1[as.numeric(iName[2])],iSingle1[-as.numeric(iName[2])])#;iOrder1
              indSCalls$m[iMerge,jm|j.fam] <- indSCalls$running[iOrder1,jm|j.fam]
              iSingle2 <- which(sub('_Para[12]','',rownames(paraCalls$m)) %in% iName[1])
              iOrder2 <- c(iSingle2[as.numeric(iName[2])],iSingle2[-as.numeric(iName[2])])
              mergedCalls[iMerge,jm] <- paraCalls$m[iOrder2,jm]
              if (any(!is.na(positionFemale[iMerge,j])))
                stop('Attempt to overwrite positionFemale[iMerge,j]...')
              positionFemale[iMerge[1],j] <- medianRF[rankFit[1]]
             } #if any
          }else if(nHits==1){
            sumGoodC <- tapply(rr$Ratio>=mO$ratioLim & rr$nNum>=mO$offspringLim,rr$Chromosome,sum)#; sumGoodC
            sumGoodC <- sumGoodC[order(as.numeric(names(sumGoodC)))]   #order chromosomes increasingly
            if (any(sumGoodC)){
              iCand <- which(sumGoodC>0)
              totalDen <- bestNum <- medianRF <- linkRatio <- rep(0,length(iCand))
              for (k in 1:length(iCand)){   #Test all possible alternatives
                #indC <- rr$Ssa %in% paste('ssa',sprintf('%.02d',iCand[k]),sep='')#; rr[indC,]
                indC <- rr$Chromosome %in% iCand[k]#; rr[indC,]
                indGood <- rr[indC,'Ratio']>=mO$ratioLim & rr[indC,'nNum']>=mO$offspringLim#; rr[indC,][indGood,]
                while (any(indGood)){
                  medGPos <- quantile(rr[indC,'Female'][indGood],na.rm=TRUE,probs=.5,type=1)#;medGPos #No averaging!
                  indCand <- rr[indC,'Female']>medGPos-mO$rngLD & rr[indC,'Female']<medGPos+mO$rngLD#; rr[indC,][indCand,]
                  lRat <- sum(rr[indC,'nNum'][indCand])/sum(rr[indC,'nDen'][indCand])#; lRat
                  bNum <- max(rr[indC,'nNum'][indCand][rr[indC,'Ratio'][indCand]>=mO$ratioLim],na.rm=TRUE)#; bNum
                  #snps[indC,][indCand,]
                  tDen <- sum(apply(!is.na(matrix(matrix(snps[indC,],nrow=sum(indC))[indCand,],nrow=sum(indCand))),2,any))#;tDen
                  if (lRat>linkRatio[k] | (lRat==linkRatio[k] & tDen>totalDen[k])){
                    linkRatio[k] <- lRat
                    medianRF[k] <- medGPos
                    totalDen[k] <- tDen
                    bestNum[k] <- bNum
                  } #if
                  indGood[indGood]=!indCand[indGood]#; rr[indC,][indGood,]
                } #while
              } #for k
              #totalDen;bestNum;medianRF;linkRatio
              compFrame <- data.frame(cbind(linkRatio,totalDen,bestNum,medianRF),row.names=names(iCand))
              if (any(linkRatio>=mO$ratioLim & totalDen>=mO$offspringLim)){
                expChrom <- as.numeric(iCand%in%altC[iName[1],1])#;expChrom
                rankFit <- order(linkRatio,totalDen,expChrom,decreasing=TRUE)#;rankFit  #NB! Favours exp. chrom.
                #partnerCandidates[[iName[1]]][[indMothers[j]]] <- compFrame[rankFit,]
                if (iCand[rankFit[1]]%in%altC[iName[1],1]){
                  iMerge <- paste(iName[1],'_Chrom',altC[iName[1],nHits],sep='')
                  indSCalls$m[iMerge,jm|j.fam] <- indSCalls$running[iMerge,jm|j.fam]+as.numeric(iName[2])-1
                  mergedCalls[iMerge,jm] <- paraCalls$m[i,jm]
                  if (any(!is.na(positionFemale[iMerge,j])))
                    stop('Attempt to overwrite positionFemale[iMerge,j]...')
                  positionFemale[iMerge,j] <- medianRF[rankFit[1]]
                } #if iCand
              } #if any linkRatio
            } #if any sumGoodC
          } #if nHits
        } #if sum
      } #if polymorphic
    } #for j
    
    for (j in 1:nFathers){    #Inheritance from fathers
      jf <- pData(BSSnp)$PedigreeID %in% paste(substr(indFathers[j],1,1),'000',substr(indFathers[j],2,4),'00',sep='')#; paraCalls$f[i,jf] %in% .5
      if (paraCalls$f[i,jf] %in% .5){    #Father j heterozygous in paralogue i
        iHetSnp <- assayData(BSSnp)$call[,jf] %in% .5   #Polymorphic markers
        j.fam <- pp %in% indFathers[j] &! jf
        iSamp <- !is.na(paraCalls$f[i,j.fam])
        iName <- strsplit(rownames(paraCalls$f)[i],'_Para')[[1]]
        nHits <- sum(!is.na(altC[iName[1],]))
        #sum(iSamp);nHits;filledIn[iName[1],jf]
        if (sum(iSamp)>=mO$offspringLim & nHits>0 &! filledIn[iName[1],jf]){
          snps <- as.matrix(inheritP$f[iHetSnp,j.fam][,iSamp],nrow=sum(iHetSnp))
          nDen <- apply(snps,1,function(x) sum(!is.na(x)))
          nNum <- apply(snps,1,function(x,p) max(sum(x==p,na.rm=TRUE),sum(x==1-p,na.rm=TRUE)),paraCalls$f[i,j.fam][iSamp])
          #rr <- cbind(lMap[iHetSnp,],nNum,nDen,Ratio=nNum/nDen)
          rr <- cbind(fData(BSSnp)[iHetSnp,c('Chromosome','Male')],nNum,nDen,Ratio=nNum/nDen)
          if (nHits==2){
            bestNum <- totalDen <- medianRF <- linkRatio <- rep(0,nHits)
            for (k in 1:nHits){   #Test both alternatives
              #indC <- rr$Ssa %in% paste('ssa',sprintf('%.02d',altC[iName[1],k]),sep='')#; rr[indC,]
              indC <- rr$Chromosome %in% altC[iName[1],k]#; rr[indC,]
              indGood <- rr[indC,'Ratio']>=mO$ratioLim & rr[indC,'nNum']>=mO$offspringLim#; rr[indC,][indGood,]
              while (any(indGood)){
                medGPos <- quantile(rr[indC,'Male'][indGood],na.rm=TRUE,probs=.5,type=1)#;medGPos #No averaging!
                indCand <- rr[indC,'Male']>medGPos-mO$rngLD & rr[indC,'Male']<medGPos+mO$rngLD#; rr[indC,][indCand,]
                lRat <- sum(rr[indC,'nNum'][indCand])/sum(rr[indC,'nDen'][indCand])#;lRat
                bNum <- max(rr[indC,'nNum'][indCand][rr[indC,'Ratio'][indCand]>=mO$ratioLim],na.rm=TRUE)#;bNum
                #snps[indC,][indCand,]
                tDen <- sum(apply(!is.na(matrix(matrix(snps[indC,],nrow=sum(indC))[indCand,],nrow=sum(indCand))),2,any))#;tDen
                if (lRat>linkRatio[k] | (lRat==linkRatio[k] & tDen>totalDen[k]) | (lRat==linkRatio[k] & tDen==totalDen[k] & bNum>bestNum[k])){
                  linkRatio[k] <- lRat
                  medianRF[k] <- medGPos
                  totalDen[k] <- tDen
                  bestNum[k] <- bNum
                } #if
                indGood[indGood]=!indCand[indGood]; rr[indC,][indGood,]
              } #while
            } #for k
            #totalDen;bestNum;medianRF;linkRatio
            if (any(linkRatio>=mO$ratioLim & totalDen>=mO$offspringLim) &! (diff(linkRatio)==0 & diff(totalDen)==0)){
              rankFit <- order(linkRatio,totalDen,decreasing=TRUE)#;rankFit
              iMerge <- paste(iName[1],'_Chrom',altC[iName[1],rankFit],sep='')#;iMerge
              iSingle1 <- which(names(iRunning) %in% iName[1])#;iSingle1
              iOrder1 <- c(iSingle1[as.numeric(iName[2])],iSingle1[-as.numeric(iName[2])])#;iOrder1
              indSCalls$f[iMerge,jf|j.fam] <- indSCalls$running[iOrder1,jf|j.fam]
              iSingle2 <- which(sub('_Para[12]','',rownames(paraCalls$f)) %in% iName[1])#;iSingle2
              iOrder2 <- c(iSingle2[as.numeric(iName[2])],iSingle2[-as.numeric(iName[2])])#;iOrder2
              mergedCalls[iMerge,jf] <- paraCalls$f[iOrder2,jf]
              if (any(!is.na(positionMale[iMerge,j])))
                stop('Attempt to overwrite positionMale[iMerge,j]...')
              positionMale[iMerge[1],j] <- medianRF[rankFit[1]]
          } #if any
          }else if(nHits==1){
            sumGoodC <- tapply(rr$Ratio>=mO$ratioLim & rr$nNum>=mO$offspringLim,rr$Chromosome,sum)#; sumGoodC
            sumGoodC <- sumGoodC[order(as.numeric(names(sumGoodC)))]   #order chromosomes increasingly
            if (any(sumGoodC)){
              iCand <- which(sumGoodC>0)
              totalDen <- bestNum <- medianRF <- linkRatio <- rep(0,length(iCand))
              for (k in 1:length(iCand)){   #Test all possible alternatives
                #indC <- rr$Ssa %in% paste('ssa',sprintf('%.02d',iCand[k]),sep='')#; rr[indC,]
                indC <- rr$Chromosome %in% iCand[k]
                indGood <- rr[indC,'Ratio']>=mO$ratioLim & rr[indC,'nNum']>=mO$offspringLim#; rr[indC,][indGood,]
                while (any(indGood)){
                  medGPos <- quantile(rr[indC,'Male'][indGood],na.rm=TRUE,probs=.5,type=1)#;medGPos #No averaging!
                  indCand <- rr[indC,'Male']>medGPos-mO$rngLD & rr[indC,'Male']<medGPos+mO$rngLD#; rr[indC,][indCand,]
                  lRat <- sum(rr[indC,'nNum'][indCand])/sum(rr[indC,'nDen'][indCand])#; lRat
                  bNum <- max(rr[indC,'nNum'][indCand][rr[indC,'Ratio'][indCand]>=mO$ratioLim],na.rm=TRUE)#; bNum
                  #snps[indC,][indCand,]
                  tDen <- sum(apply(!is.na(matrix(matrix(snps[indC,],nrow=sum(indC))[indCand,],nrow=sum(indCand))),2,any))#;tDen
                  if (lRat>linkRatio[k] | (lRat==linkRatio[k] & tDen>totalDen[k])){
                    linkRatio[k] <- lRat
                    medianRF[k] <- medGPos
                    totalDen[k] <- tDen
                    bestNum[k] <- bNum
                  } #if
                  indGood[indGood]=!indCand[indGood]; rr[indC,][indGood,]
                } #while
              } #for k
              #totalDen;bestNum;medianRF;linkRatio
              compFrame <- data.frame(cbind(linkRatio,totalDen,bestNum,medianRF),row.names=names(iCand))
              if (any(linkRatio>=mO$ratioLim & totalDen>=mO$offspringLim)){
                expChrom <- as.numeric(iCand%in%altC[iName[1],1])#;expChrom
                rankFit <- order(linkRatio,totalDen,expChrom,decreasing=TRUE)#;rankFit  #NB! Favours exp. chrom.
                #partnerCandidates[[iName[1]]][[indFathers[j]]] <- compFrame[rankFit,]
                if (iCand[rankFit[1]]%in%altC[iName[1],1]){
                  iMerge <- paste(iName[1],'_Chrom',altC[iName[1],nHits],sep='')
                  indSCalls$f[iMerge,jf|j.fam] <- indSCalls$running[iMerge,jf|j.fam]+as.numeric(iName[2])-1
                  mergedCalls[iMerge,jf] <- paraCalls$f[i,jf]
                  if (any(!is.na(positionMale[iMerge,j])))
                    stop('Attempt to overwrite positionMale[iMerge,j]...')
                  positionMale[iMerge,j] <- medianRF[rankFit[1]]
                } #if iCand
              } #if any linkRatio
            } #if any sumGoodC
          } #if nHits
        } #if sum
      } #if polymorphic
    } #for j
  } #for i

  #Fill in for offspring within each full-sib family
  indOffspring <- !substr(pData(BSSnp)$PedigreeID,8,9)%in%'00'
  indFullsib <- unique(substr(pData(BSSnp)$PedigreeID[indOffspring],1,7))
  nFullsib <- length(indFullsib)
  str1 <- 'Merging calls for full-sib family %d/%d'
  for (j in 1:nFullsib){
    if (j %in% seq(1,nFullsib,10))
      message(sprintf(str1,j,nFullsib))
    j.fam <- substr(pData(BSSnp)$PedigreeID,1,7)%in%indFullsib[j]
    indCurr <- substr(pData(BSSnp)$PedigreeID[j.fam][1],1,7)
    jm <- pData(BSSnp)$PedigreeID %in% paste(substr(indCurr,1,4),'00000',sep='')
    jf <- pData(BSSnp)$PedigreeID %in% paste(substr(indCurr,1,1),'000',substr(indCurr,5,7),'00',sep='')
    if (any(jm) & any(jf)){
      indM <- as.vector(indSCalls$m[,j.fam])
      indF <- as.vector(indSCalls$f[,j.fam])
      mergedCalls[,j.fam] <- (paraCalls$m[indM]+paraCalls$f[indF])/2
      pCalls <- mergedCalls[,c(which(jm),which(jf))]#;pCalls[1:10,]
      indAA <- apply(pCalls==0,1,sum)%in%2  #NB! The following may mask NA's which are due to ped-error
      mergedCalls[indAA,j.fam][is.na(mergedCalls[indAA,j.fam])] <- 0  #all AA
      indBB <- apply(pCalls==1,1,sum)%in%2
      mergedCalls[indBB,j.fam][is.na(mergedCalls[indBB,j.fam])] <- 1  #all BB
      indAB <- (pCalls[,1]%in%0 & pCalls[,2]%in%1) | (pCalls[,1]%in%1 & pCalls[,2]%in%0)
      mergedCalls[indAB,j.fam][is.na(mergedCalls[indAB,j.fam])] <- .5   #all AB
      for (k in 1:nNew){   #Improve speed by looping through MSV-5's instead
        pairInd <- rNames[markerNames%in%markerNames[k]]#;pairInd
        step <- which(pairInd%in%rNames[k])#;step
        #mergedCalls[pairInd,c(which(jm),which(jf),which(j.fam))]
        #assayData(BSRed)$call[markerNames[k],c(which(jm),which(jf),which(j.fam))]
        if (length(pairInd)==2){
          indNA <- is.na(mergedCalls[pairInd[step],j.fam]) &! is.na(mergedCalls[pairInd[-step],j.fam])#;indNA
          mergedCalls[k,j.fam][indNA] <- 2*assayData(BSRed)$call[markerNames[k],j.fam][indNA]-mergedCalls[pairInd[-step],j.fam][indNA]  #NB! Possibility of calls >1 or <0 ???
        } #if length
        indNA1 <- is.na(mergedCalls[k,j.fam])#;indNA1
        if (any(indNA1)){
          ii <- is.na(pCalls[k,])#;ii
          if (any(ii)){ # &step!=2
            sInd <- c(sub('Chrom[0-9XY]+','Para1',pairInd[1]),sub('Chrom[0-9XY]+','Para2',pairInd[1]))
            #paraCalls$m[sInd,c(which(jm),which(jf),which(j.fam))]
            #paraCalls$m[indSCalls$m[pairInd,jm]]
            rankFitM <- order(indSCalls$m[pairInd,jm])#;rankFitM
            if (ii[1] | length(pairInd)==1){
              iim <- apply(matrix(matrix(paraCalls$m[sInd,j.fam],nrow=length(pairInd))[,indNA1],nrow=length(pairInd)),2,function(x) diff(x)%in%0)
            }else{
              iim <- !is.na(paraCalls$m[sInd[1],j.fam][indNA1])
            } #if ii[1]
            #paraCalls$f[sInd,c(which(jm),which(jf),which(j.fam))]
            #paraCalls$f[indSCalls$f[pairInd,jf]]
            rankFitF <- order(indSCalls$f[pairInd,jf])#;rankFitF
            if (ii[2] | length(pairInd)==1){
              iif <- apply(matrix(matrix(paraCalls$f[sInd,j.fam],nrow=length(pairInd))[,indNA1],nrow=length(pairInd)),2,function(x) diff(x)%in%0)
            }else{
              iif <- !is.na(paraCalls$f[sInd[1],j.fam][indNA1])
            } #if ii[2]
            if (any(iif&iim) & step==2) stop('Cannot include step==1 -criterion...')
            mergedCalls[k,j.fam][indNA1][iim&iif] <- (paraCalls$m[sInd[rankFitM[step]],j.fam]+paraCalls$f[sInd[rankFitF[step]],j.fam])[iim&iif]/2
          }else{
            pp <- assayData(BSRed)$call[markerNames[k],c(which(jm),which(jf))]#;pp
            rngParents <- c(max(0,max(pp)-.5),min(1,min(pp)+.5))
            if (all(pp<=.25)){
              rngParents[2] <- sum(pp)
            }else if (all(pp>=.75)){
              rngParents[1] <- sum(pp)-1
            }
            #rngParents
            indLow <- assayData(BSRed)$call[markerNames[k],j.fam][indNA1] %in% rngParents[1]#;indLow
            indHigh <- assayData(BSRed)$call[markerNames[k],j.fam][indNA1] %in% rngParents[2]#;indHigh
            if (any(indLow|indHigh) & step==2) stop('Cannot include step==1 -criterion...')
            mergedCalls[k,j.fam][indNA1][indLow] <- max(pCalls[k,])-.5  #NB! Possibility of calls >1 or <0 ???
            mergedCalls[k,j.fam][indNA1][indHigh] <- min(pCalls[k,])+.5  #NB! Possibility of calls >1 or <0 ???
          } #if any ii
        } #if any indNA1
      } #for k
    } #if any
  } #for j
  x <- NULL
  x$x <- mergedCalls
  x$chromPairs <- altC
  x$positionFemale <- positionFemale
  x$positionMale <- positionMale
  x
}



preprocessBeadSet <- function(BSData,normInd,normOpts=setNormOptions(shearInf1=!is.null(normInd))){
  if ('noiseIntensity'%in%colnames(pData(BSData)))
    stop("Column 'noiseIntensity' found in phenoData. Will not pre-process.")
  if (normOpts$shearInf1){
    indTrain <- list()
    indRed <- grepl('10[0-9]',colnames(normInd))
    indTrain$'10r' <- apply(as.matrix(normInd[,indRed]),1,any)
    indGreen <- grepl('20[0-9]',colnames(normInd))
    indTrain$'20g' <- apply(as.matrix(normInd[,indGreen]),1,any)
    #indTrain$'10r' <- normInd[,'102']|normInd[,'101']  #Inf-I using red channel
    #indTrain$'20g' <- normInd[,'202']|normInd[,'201']  #Inf-I using green channel
  }else{
    indTrain <- rep(TRUE,nrow(BSData))
  } #end if
  BSData <- shearRawSignal(BSData,normOpts=normOpts,indTrain=indTrain)
  message('Estimating origin and noise-levels...')
  noiseDist <- getNoiseDistributions(BSData,normInd=normInd,normOpts=normOpts)
  message(sprintf('Transformation of signal: %s...',normOpts$trans))
  trChannel <- transformChannels(assayData(BSData)$R,assayData(BSData)$G,normOpts=normOpts)
  message(sprintf('Transformation of red SE: %s...',normOpts$trans))
  trSE.R <- transformSEs(assayData(BSData)$R,assayData(BSData)$se.R,normOpts=normOpts)
  message(sprintf('Transformation of green SE: %s...',normOpts$trans))
  trSE.G <- transformSEs(assayData(BSData)$G,assayData(BSData)$se.G,normOpts=normOpts)
  message(sprintf('Normalize channels: %s...',normOpts$method))
  trChannel <- normalizeShearedChannels(trChannel,noiseDist,normOpts=normOpts)
  assayData(BSData)$R <- trChannel$X
  assayData(BSData)$se.R <- trSE.R$SE
  assayData(BSData)$G <- trChannel$Y
  assayData(BSData)$se.G <- trSE.G$SE
  noiseData <- data.frame(noiseIntensity=normOpts$nSD*trChannel$X.SE,row.names=sampleNames(BSData))
  noiseMetadata <- data.frame(labelDescription=paste('Estimated red channel SD*',normOpts$nSD,sep=''),row.names=colnames(noiseData))
  noiseInt <- new('AnnotatedDataFrame',data=noiseData,varMetadata=noiseMetadata,dimLabels=dimLabels(phenoData(BSData)))
  phenoData(BSData) <- combine(phenoData(BSData),noiseInt)
  validObject(BSData)
  BSData
}



plotPreprocessing <- function(BSData,normInd,normOpts=setNormOptions(shearInf1=!is.null(normInd)),plotArray=1,...){
  pp <- par()
  par(mfrow=c(2,2),mai=c(.8,.8,.8,.4))
  if (normOpts$shearInf1){
    indTrain <- list()
    indRed <- grepl('10[0-9]',colnames(normInd))
    indTrain$'10r' <- apply(as.matrix(normInd[,indRed]),1,any)
    indGreen <- grepl('20[0-9]',colnames(normInd))
    indTrain$'20g' <- apply(as.matrix(normInd[,indGreen]),1,any)
    #indTrain$'10r' <- normInd[,'102']|normInd[,'101']  #Inf-I using red channel
    #indTrain$'20g' <- normInd[,'202']|normInd[,'201']  #Inf-I using green channel
  }else{
    indTrain <- rep(TRUE,nrow(BSData))
  } #end if
  scatterArrays(BSData,arrays=plotArray,smooth=TRUE,newFigure=FALSE,probeInd=1:nrow(BSData))
  BSsheared = shearRawSignal(BSData[,plotArray],plot=TRUE,newFigure=FALSE,normOpts=normOpts,indTrain=indTrain,verbose=TRUE)
  noiseDist = getNoiseDistributions(BSsheared,normInd=normInd,normOpts=normOpts,plot=TRUE,newFigure=FALSE)
  plotEstimatedNoise(BSsheared,noiseDist,normOpts=normOpts,newFigure=FALSE,xlab='R',ylab='G',...)
  par(pp[c('mfcol','mai')])
}



makeFilenames <- function(tag,normOpts,path='.'){
  if (!substr(path,nchar(path),nchar(path)) %in% '/'){
    path <- paste(path,'/',sep='')
  }
  intFile <- paste(path,'Int_',normOpts$transf,'_',normOpts$method,'_',tag,'.txt',sep='')
  thFile <- paste(path,'Theta_',normOpts$transf,'_',normOpts$method,'_',tag,'.txt',sep='')
  seFile <- paste(path,'SE_',normOpts$transf,'_',normOpts$method,'_',tag,'.txt',sep='')
  phFile <- paste(path,'Pheno_',normOpts$transf,'_',normOpts$method,'_',tag,'.txt',sep='')
  callFile <- paste(path,'Call_',normOpts$transf,'_',normOpts$method,'_',tag,'.txt',sep='')
  resFile <- paste(path,'Res_',normOpts$transf,'_',normOpts$method,'_',tag,'.txt',sep='')
  genoFile <- paste(path,'Alleles_',normOpts$transf,'_',normOpts$method,'_',tag,'.txt',sep='')
  dataFiles <- c(intFile=intFile,thFile=thFile,seFile=seFile,phFile=phFile,callFile=callFile,resFile=resFile,genoFile=genoFile)
  ind <- substr(dataFiles,nchar(path)+1,nchar(dataFiles)) %in% dir(path)
  if (any(ind)){
    OK = FALSE
    while (!OK){
      ovrwrt <- readline(prompt='Warning! Existing files may be overwritten. Continue? (yes/no) ')
      OK <- ovrwrt %in% c('yes','no')
    }
    if (ovrwrt %in% 'yes'){
      message('Saving to \'dataFiles\' may overwrite existing files!')
    }else{
      message("Please provide a unique tag for the filenames")
      dataFiles <- NULL
    }
  }
  dataFiles
}


modifyExistingFilenames <- function(dataFiles,oldtag,newtag){
  dataFiles['callFile'] <- sub(oldtag,newtag,dataFiles['callFile'])
  dataFiles['resFile'] <- sub(oldtag,newtag,dataFiles['resFile'])
  dataFiles['genoFile'] <- sub(oldtag,newtag,dataFiles['genoFile'])
  print(dataFiles)
  dataFiles
}



#Saves data with valid row/column names
writeAlleleSet <- function(dataFiles,BSRed,append=TRUE,sep='\t',quote=FALSE){
  sampleNames(BSRed) <- make.names(sampleNames(BSRed))
  featureNames(BSRed) <- make.names(featureNames(BSRed))
  re <- gregexpr('/',dataFiles)
  nDir <- length(re[[1]])
  path <- substr(dataFiles[[1]],1,re[[1]][nDir]-1)
  if (any(c('callFile','resFile','genoFile') %in% names(dataFiles)) & any(c('call','genotype') %in% names(assayData(BSRed)))){
    if ('callFile' %in% names(dataFiles)){
      fName <- substr(dataFiles[['callFile']],re[[1]][nDir]+1,nchar(dataFiles[['callFile']]))
      if (!fName %in% dir(path,fName)){
        message(sprintf("Writing calls to '%s'...",dataFiles['callFile']))
      }else if (append){
        message(sprintf("Appending calls to '%s'...",dataFiles['callFile']))
      }else{
        message(sprintf("Overwriting '%s'...",dataFiles['callFile']))
      }
      write.table(formatC(assayData(BSRed)$call,digits=2,format='f'),file=dataFiles['callFile'],append=append,sep=sep,quote=quote,col.names=!append)
    }
    if ('resFile' %in% names(dataFiles)){
      fName <- substr(dataFiles[['resFile']],re[[1]][nDir]+1,nchar(dataFiles[['resFile']]))
      if (!fName %in% dir(path,fName)){
        message(sprintf("Writing featureData to '%s'...",dataFiles['resFile']))
      }else if (append){
        message(sprintf("Appending featureData to '%s'...",dataFiles['resFile']))
      }else{
        message(sprintf("Overwriting '%s'...",dataFiles['resFile']))
      }
      write.table(format(fData(BSRed)),file=dataFiles['resFile'],append=append,sep=sep,quote=quote,col.names=!append)
    }
    if ('genoFile' %in% names(dataFiles)){
      fName <- substr(dataFiles[['genoFile']],re[[1]][nDir]+1,nchar(dataFiles[['genoFile']]))
      if (!fName %in% dir(path,fName)){
        message(sprintf("Writing genotypes to '%s'...",dataFiles['genoFile']))
      }else if (append){
        message(sprintf("Appending genotypes to '%s'...",dataFiles['genoFile']))
      }else{
        message(sprintf("Overwriting '%s'...",dataFiles['genoFile']))
      }
       write.table(formatC(assayData(BSRed)$genotype,digits=2,format='f'),file=dataFiles['genoFile'],append=append,sep=sep,quote=quote,col.names=!append)
    }
    iExtra = !names(dataFiles) %in% c('callFile','resFile','genoFile')
    if (any(iExtra)){
      message('Will NOT write the following data to file(s) (to force - save separately):')
      message(sprintf("'%s'\t",names(dataFiles)[iExtra]))
    }
  }else if (any(c('callFile','resFile','genoFile') %in% names(dataFiles))){
    message('No calls or genotypes available. No files written')
  }else if (any(c('intFile','thFile','seFile','phFile') %in% names(dataFiles))){
    if ('intFile' %in% names(dataFiles)){
      fName <- substr(dataFiles[['intFile']],re[[1]][nDir]+1,nchar(dataFiles[['intFile']]))
      if (!fName %in% dir(path,fName)){
        message(sprintf("Writing intensities to '%s'...",dataFiles['intFile']))
      }else if (append){
        message(sprintf("Appending intensities to '%s'...",dataFiles['intFile']))
      }else{
        message(sprintf("Overwriting '%s'...",dataFiles['intFile']))
      }
      write.table(formatC(t(assayData(BSRed)$intensity),format='f'),file=dataFiles['intFile'],append=append,sep=sep,quote=quote,col.names=!append)
    }
    if ('thFile' %in% names(dataFiles)){
      fName <- substr(dataFiles[['thFile']],re[[1]][nDir]+1,nchar(dataFiles[['thFile']]))
      if (!fName %in% dir(path,fName)){
        message(sprintf("Writing thetas to '%s'...",dataFiles['thFile']))
      }else if (append){
        message(sprintf("Appending thetas to '%s'...",dataFiles['thFile']))
      }else{
        message(sprintf("Overwriting '%s'...",dataFiles['thFile']))
      }
      write.table(formatC(t(assayData(BSRed)$theta),format='f'),file=dataFiles['thFile'],append=append,sep=sep,quote=quote,col.names=!append)
    }
    if ('seFile' %in% names(dataFiles)){
      fName <- substr(dataFiles[['seFile']],re[[1]][nDir]+1,nchar(dataFiles[['seFile']]))
      if (!fName %in% dir(path,fName)){
        message(sprintf("Writing standard errors to '%s'...",dataFiles['seFile']))
      }else if (append){
        message(sprintf("Appending standard errors to '%s'...",dataFiles['seFile']))
      }else{
        message(sprintf("Overwriting '%s'...",dataFiles['seFile']))
      }
      write.table(formatC(t(assayData(BSRed)$SE),format='f'),file=dataFiles['seFile'],append=append,sep=sep,quote=quote,col.names=!append)
    }
    if ('phFile' %in% names(dataFiles)){
      fName <- substr(dataFiles[['phFile']],re[[1]][nDir]+1,nchar(dataFiles[['phFile']]))
      if (!fName %in% dir(path,fName)){
        message(sprintf("Writing phenoData to '%s'...",dataFiles['phFile']))
      }else if (append){
        message(sprintf("Appending phenoData to '%s'...",dataFiles['phFile']))
      }else{
        message(sprintf("Overwriting '%s'...",dataFiles['phFile']))
      }
      write.table(format(pData(BSRed)),file=dataFiles['phFile'],append=append,sep=sep,quote=quote,col.names=!append)
    }
  }else{
    message('No valid file specified')
  }
}

