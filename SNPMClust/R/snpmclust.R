snpmclust <- function(indata,p=1,priorfrac=.2,uncertcutoff=.01,qcutoff=0,showplots=FALSE,
  xm1=NA,xm2=NA,xm3=NA,ym1=NA,ym2=NA,ym3=NA,ranseed=1969,R.lowcutoff=.05) {

  ## compute number of prior points
  R = indata$R.trans[p,]
  theta = indata$logratio[p,]
  R[indata$R[p,]<R.lowcutoff] = NA
  R[R=="NaN"] <- NA
  R[is.element(theta,c(-Inf,Inf,"NaN"))] <- NA
  R[is.na(theta)] <- NA
  theta[is.na(R)] <- NA
  priorpoints <- round(length(theta)*priorfrac)
  ## if either data vector is empty, return null results
  if (sum(!is.na(R))==0 | sum(!is.na(theta))==0) {
    outset <- data.frame(indata$SNP[p],indata$SampleID,as.character(rep("NC",length(R))),
      as.numeric(rep(NA,length(R))),stringsAsFactors=FALSE)
    names(outset) <- c("SNP","SampleID","MClustCalls","Uncertainty")
    return(list(calls=outset,snp=indata$SNP[p],callrate=0,priorfrac=priorfrac,
      uncertcutoff=uncertcutoff,qcutoff=qcutoff))
  }

  ## compute mean and variance of homozygous clusters
  xtmp = theta
  Rtmp = R	
  keepit <- !is.na(theta) & !is.na(R)
  xtmp <- theta[keepit]
  Rtmp <- R[keepit]
  gscalltmp <- indata$GType[p,][keepit]
  lowesstmp <- lowess(xtmp,Rtmp)
  leftmean <- lowesstmp$y[order(lowesstmp$x)][1]
  rightmean <- lowesstmp$y[order(lowesstmp$x)][length(lowesstmp$y)]
  tmpresid <- Rtmp[order(xtmp)]-lowesstmp$y[order(lowesstmp$x)]
  noisesd <- sd(tmpresid)
  pseudodata <- generatepriors(xtmp,Rtmp,gscalltmp,priorpoints,xm1=xm1,xm2=xm2,
    xm3=xm3,ym1=ym1,ym2=ym2,ym3=ym3,ranseed=ranseed)
  ## add pseudodata
  xtrans <- theta
  Rtrans <- R
  if (is.na(pseudodata)[1]) {
    snpdat <- data.frame(xtrans[!is.na(xtrans) & !is.na(Rtrans)],  
      Rtrans[!is.na(xtrans) & !is.na(Rtrans)])
  } else {
    snpdat <- data.frame(c(xtrans[!is.na(xtrans) & !is.na(Rtrans)],pseudodata[,1]),  
      c(Rtrans[!is.na(xtrans) & !is.na(Rtrans)],pseudodata[,2]))
  }  
  if (showplots) {
    plot(xtrans,Rtrans)
    if (!is.na(pseudodata)[1]) {
      points(pseudodata[,1],pseudodata[,2],col=2,pch=3)
    }
  }

  ## determine number of clusters
  tmptab <- table(indata$GType[p,][is.element(indata$GType[p,],c("AA","AB","BB"))])
  tmptab <- tmptab[tmptab>0]
  numclust <- length(tmptab)
  if (numclust==0) {
    numclust <- 3
    tmptab <- table(c("AA","AB","BB"))
  }

  ## run main routine and output results
  callrate <- NA
  if (dim(unique(snpdat))[1]>5) {
    snpmclust <- Mclust(snpdat,G=numclust:numclust,modelNames=c("EEE","EEV"))
    classif <- rep(NA,length(xtrans))
    uncert <- rep(NA,length(xtrans))
    classif[!is.na(xtrans) & !is.na(Rtrans)] <- 
      snpmclust$classification[1:length(xtrans[!is.na(xtrans) & !is.na(Rtrans)])]
    uncert[!is.na(xtrans) & !is.na(Rtrans)] <- 
      snpmclust$uncertainty[1:length(xtrans[!is.na(xtrans) & !is.na(Rtrans)])]
    classif2 <- classif
    classif2[uncert>uncertcutoff] <- 4
    reorder <- rank(snpmclust$parameters$mean[1,])
    classif2 <- names(tmptab)[reorder[classif2]]
    classif2[is.na(classif2)] <- 4
    ## apply quantile cutoff
    if (qcutoff>0 && qcutoff<1) {  
      qtmp = quantile(uncert,qcutoff,na.rm=TRUE)
      uncert = sapply(uncert,function(x) max(x,qtmp))
      classif2[uncert>uncertcutoff] <- 4
    }
    callrate <- sum(table(classif2)["AA"],table(classif2)["AB"],table(classif2)["BB"],na.rm=TRUE)/
      length(classif2)
    names(callrate) = NULL
    if (showplots) {
      plot(snpmclust)
      plot(theta,R,type="n",main=paste(indata$SNP[p],
        "\n uncertcutoff=",uncertcutoff,"; qcutoff=",qcutoff,
        "; call rate=",sprintf("%.3f",callrate),sep=""),xlab="log-ratio",ylab="R.trans")
      points(theta[classif2=="AA"],R[classif2=="AA"],col=6,pch=5)
      points(theta[classif2=="AB"],R[classif2=="AB"],col=4,pch=2)
      points(theta[classif2=="BB"],R[classif2=="BB"],col=2,pch=6)
      points(theta[classif2=="4"],R[classif2=="4"],col=1,pch=4)
    }
    calltxt <- rep("NC",length(classif2))
    calltxt[classif2=="AA"] <- "AA"
    calltxt[classif2=="AB"] <- "AB"
    calltxt[classif2=="BB"] <- "BB"
  } else {
    calltxt <- rep("NC",length(theta))
    uncert <- rep(1,length(theta))
  }
  outset <- data.frame(indata$SNP[p],indata$SampleID,calltxt,uncert,stringsAsFactors=FALSE)
  names(outset) <- c("SNP","SampleID","MClustCalls","Uncertainty")
  return(list(calls=outset,snp=indata$SNP[p],callrate=callrate,priorfrac=priorfrac,
    uncertcutoff=uncertcutoff,qcutoff=qcutoff))
}