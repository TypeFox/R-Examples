
.add2matrix <-
function(values,xmat=NULL){
  
  #xmat   - n X ? matrix with one row, columns are integer values
  #values - length-n vector be added/slotted in to xvec
  
  if(is.null(xmat)){
    n    <- length(values)
    cc   <- sort(unique(values))
    xmat <- matrix(0,n,length(cc),dimnames = list(1:n,cc))
    xmat[ cbind( c(1:n),match(values,cc)) ] <- 1
    return(xmat)
  }
  
  n <- nrow(xmat)
  if(length(values) != n)stop('vector length must equal rows in xmat')
  
  all <- sort( unique( c(values,as.numeric(colnames(xmat))) ))
  nc       <- length(all)
  
  xnew <- matrix(0,n,nc,dimnames = list(1:n,all))
  xnew[,colnames(xmat)] <- xmat
  
  xnew[ cbind(c(1:n),match(values,all)) ] <- xnew[ cbind(c(1:n),match(values,all)) ] + 1
  xnew
}
.appendMatrix <-
function(m1,m2,fill=NA,SORT=F,asNumbers=F){  
  
  # matches matrices by column names
  # asNumbers: if column heads are numbers and SORT, then sort numerically

   if(length(m1) == 0){
     if(is.matrix(m2)){
       m3 <- m2
     } else {
       m3 <- matrix(m2,nrow=1)
     }
     if( !is.null(names(m2)) )colnames(m3) <- names(m2)
     return(m3)
   }
   if(length(m2) == 0){
     if(!is.matrix(m1))m1 <- matrix(m1,nrow=1)
     return(m1)
   }
   if( is.vector(m1) | (length(m1) > 0 & !is.matrix(m1)) ){
     nn <- names(m1)
     if(is.null(nn))warning('cannot append matrix without names')
     m1 <- matrix(m1,1)
     colnames(m1) <- nn
   }  
   if( is.vector(m2) | (length(m2) > 0 & !is.matrix(m2)) ){
     nn <- names(m2)
     if(is.null(nn))warning('cannot append matrix without names')
     m2 <- matrix(m2,1)
     colnames(m2) <- nn
   }

   c1 <- colnames(m1)
   c2 <- colnames(m2)
   r1 <- rownames(m1)
   r2 <- rownames(m2)
   n1 <- nrow(m1)
   n2 <- nrow(m2)

   allc <-  unique( c(c1,c2) ) 
   if(SORT & !asNumbers)allc <- sort(allc)
   if(SORT & asNumbers){
     ac <- as.numeric(allc)
     allc <- as.character( sort(ac) )
   }

   nr <- n1 + n2
   nc <- length(allc)

   if(is.null(r1))r1 <- paste('r',c(1:n1),sep='-')
   if(is.null(r2))r2 <- paste('r',c((n1+1):nr),sep='-')
   new <- c(r1,r2)

   mat1 <- match(c1,allc)
   mat2 <- match(c2,allc)

   out <- matrix(fill,nr,nc)
   colnames(out) <- allc
   rownames(out) <- new

   out[1:n1,mat1] <- m1
   out[(n1+1):nr,mat2] <- m2
   out
}
.byIndex <-
function(xx,INDICES,FUN,coerce=F,...){  
  
#INDICES is list, each same length as  x
  
#  fun <- match.fun(FUN)
  
  nl <- length(INDICES)
  
  tmp  <-  unlist(by( as.vector(xx),INDICES,FUN,...) ) 
  nd   <- dim(tmp)
  tmp  <- array(tmp,dim=nd, dimnames=dimnames(tmp))
  
  tmp[is.na(tmp)] <- 0
  
  if(!coerce)return(tmp)
  
  dname <- dimnames(tmp)
  mk    <- rep(0,length(nd))
  
  for(k in 1:length(nd))mk[k] <- max(as.numeric(dimnames(tmp)[[k]]))
  
  wk <- which(mk > nd)
  if(length(wk) > 0){
    tnew  <- array(0,dim=mk)
    if(length(dim(tnew)) == 1)tnew <- matrix(tnew,dim(tnew),1)
    for(k in wk){
      newk <- c(1:mk[k])
      mat  <- match(dimnames(tmp)[[k]],newk)
      if(k == 1){
        tnew[mat,] <- tmp
        rownames(tnew) <- 1:nrow(tnew)
      }
      if(k == 2){
        tnew[,mat] <- tmp
        colnames(tnew) <- c(1:ncol(tnew))
      }
      tmp <- tnew
    }
  }
  tmp
}
.chains2density <-
function(chainMat,labs=NULL,reverseM=F,varName=NULL){
  
  #assumes column names are varName or 'something_varname'
  
  #chainMat - MCMC output [samples,chains]
  
  chNames <- colnames(chainMat)
  
  if(!is.null(varName)){
    wc <- grep(varName,colnames(chainMat))
    if(length(wc) == 0)stop('varName not found in colnames(chainMat)')
    
    ww <- grep('_',colnames(chainMat))
    if(length(ww) > 0){
      tmp <- matrix( unlist(strsplit(colnames(chainMat),'_')),ncol=2,byrow=T) 
      wc <- which(tmp[,2] == varName)
      if(length(wc) == 0)wc <- which(tmp[,1] == varName)
    }
    chainMat <- chainMat[,wc]
    if(!is.matrix(chainMat))chainMat <- matrix(chainMat,ncol=1)
    colnames(chainMat) <- chNames[wc]
  }
  
  nj <- ncol(chainMat)
  nd <- 512
  
  clab <- colnames(chainMat)
  if(is.null(labs) & !is.null(clab))labs <- clab
  if(is.null(labs) & is.null(clab)) labs <- paste('v',c(1:nj),sep='-')
  
  xt <- yt <- matrix(NA,nj,nd)
  rownames(xt) <- rownames(yt) <- labs
  
  xrange <- signif(range(chainMat),2)
  
  for(j in 1:nj){
    
 #   lj  <- labs[j]
    xj  <- chainMat[,j]
    tmp <- density(xj,n = nd, cut=0, na.rm=T)
    xt[j,]  <- tmp$x
    yt[j,]  <- tmp$y
    
  }
  yymax <- max(yt,na.rm=T)
  
  if(reverseM){
    xt <- -t( apply(xt,1,rev) )
    yt <- t( apply(yt,1,rev) )
  }
  
  list(x = xt, y = yt, xrange = xrange, ymax = yymax, chainMat = chainMat)
}
.checkDesign <-
function(x,intName='intercept',xflag=':'){  # name of intercept column

  # xflag - indicates that variable is an interaction
  
  p <- ncol(x)
  
  if(ncol(x) < 3){
    message('no design check, x has 2 columns')
    return( list(VIF = 0, correlation = 1, rank = 2, p = 2, isFactor=character(0)) )
  }
    
  if(is.null(colnames(x))){
    colnames(x) <- paste('x',c(1:p),sep='_')
  }
  xrange      <- apply(x,2,range)
  wi          <- which(xrange[1,] == 1 & xrange[2,] == 1)
  if(length(wi) > 0)colnames(x)[wi] <- 'intercept'
  
  wx <- grep(xflag,colnames(x))
  wi <- which(colnames(x) == 'intercept')
  wi <- unique(c(wi,wx))

  xname <- colnames(x)
  
  wmiss <- which(is.na(x),arr.ind=T)
  
  if(length(wmiss) > 0){
    rowTab <- table( table(wmiss[,1]) )
    colTab <- table(wmiss[,2])
  }
    
  VIF <- rep(NA,p)
  names(VIF) <- xname
  
  isFactor <- character(0)
  
  for(k in 1:p){

    if(xname[k] %in% wi)next
    
    notk <- xname[xname != xname[k] & !xname %in% xname[wi]]
    ykk  <- x[,xname[k]]
    xkk  <- x[,notk]

    tkk    <- summary(lm(ykk ~ xkk))$adj.r.squared
    VIF[k] <- 1/(1 - tkk)
    
    xu <- sort( unique(x[,k]) )
    tmp <- identical(c(0,1),xu)
    if(tmp)isFactor <- c(isFactor,xname[k])
  }

  VIF <- VIF[-wi] 

  corx <- cor(x[,-wi])
  rankx <- qr(x)$rank
  corx[upper.tri(corx,diag=T)] <- NA
  
  findex <- rep(0,p)
  
  findex[xname %in% isFactor] <- 1
  
  designTable <- list('table' = rbind( round(VIF,1),findex[-wi],round(corx,2)) )
  rownames(designTable$table) <- c('VIF','factor',xname[-wi])
  
  designTable$table <- designTable$table[-3,]
  
  if(p == rankx)designTable$rank <- paste('full rank:',rankx,'= ncol(x)')
  if(p < rankx) designTable$rank <- paste('not full rank:',rankx,'< ncol(x)')

  list(VIF = round(VIF,1), correlation = round(corx,2), rank = rankx, p = p,
       isFactor = isFactor, designTable = designTable)
}
.clusterPlot <-
function(dcor=NULL,dist=NULL,main=' ',xlab='Species',method='complete',
                        cex=1,ncluster=2, add=F,
                        xlim=NULL, colCode = NULL,horiz=T,textSize=1){  
  
  #dcor is a correlation matrix
  #dist is a distance matrix

  
  if(!is.null(dist)){
    nn <- nrow(dist)
    diss <- as.dist( dist )
  }
  if(is.null(dist)){
    nn <- nrow(dcor)
    diss <- as.dist(.cov2Dist(dcor))
  }

  tmp    <- hclust(diss,method)
  corder <- tmp$order
  ctmp   <- cutree(tmp,k=1:ncluster)
  
  wclus <- ctmp[,ncluster]
  clusterCol <- NULL
  
  clusterIndex <- ctmp[,ncluster]
  
  clusterList <- character(0)
  
#  mycols <- mapColors(ncluster)
  
  colF   <- colorRampPalette(c('black','blue','orange','brown','red'))
  mycols   <- colF(ncluster)
  
  if(is.null(colCode)){
    colCode <- mycols[ctmp[,ncluster]]
    names(colCode) <- rownames(ctmp)
  }
  
  colLab <- function(n) {
 
    if(is.leaf(n)) {
      a <- attributes(n)
      
      attr(n, "nodePar") <-
        c(a$nodePar, list(col = colCode[n[1]],lab.col = colCode[n[1]]))
    }
    n
  }
  
  tdendro <- as.dendrogram(tmp)
  dL      <- dendrapply(tdendro,colLab)
  
 # textSize <- exp(-.01*nn)
  
  new <- F
  if(add)new <- T
  par(new = new,cex=textSize)
  tmp <- plot(dL,nodePar=list(cex=.1,lab.cex=textSize),horiz=horiz,xlim=xlim)
  title(main)
  
  invisible(list( clusterList = clusterList, colCode = colCode, clusterIndex = clusterIndex,
                  corder = corder) )

}
.colorLegend <-
function(xx,yy,ytick=NULL,scale=seq(yy[1],yy[2],length=length(cols)),
                        cols,labside='right',
                        text.col=NULL,
                        bg=NULL,endLabels=NULL){  
  # xx = (x1,x2), y = (y1,y2)
  # bg is color of border

  nn <- length(scale)
  ys <- seq(yy[1],yy[2],length=nn)

  for(j in 1:(length(scale)-1)){
    rect(xx[1],ys[j],xx[2],ys[j+1],col=cols[j],border=NA)
  }
  if(!is.null(bg))rect(xx[1],yy[1],xx[2],yy[2],border=bg,lwd=3)
  if(!is.null(ytick)){
    
    ys <- diff(yy)/diff(range(ytick))*ytick
    yt <- ys - min(ys) + yy[1]
    
    for(j in 1:length(yt)){
      lines(xx,yt[c(j,j)])
    }
  }
  if(!is.null(endLabels)){ 
    cx <- cols[c(1,nn)]
    if(!is.null(text.col))cx <- text.col
    if(!is.null(text.col))cx <- text.col
    if(labside == 'right')text(diff(xx)+c(xx[2],xx[2]),yy,endLabels,col=cx)
    if(labside == 'left')text(c(xx[1],xx[1]),yy,endLabels,pos=2,col=cx)
  }
}
.corMatCI <- function(rmat, nmat, alpha = .05){
  
  # rmat   - m by m correlation matrix
  # nmat   - m by m matrix of counts
  # lo, hi - lower, upper CI
  # sig    - rmat outside CI
  
  m <- nrow(rmat)
  if(!is.matrix(nmat))nmat <- matrix(nmat,m,m)
  
  lo <- hi <- sr <- nmat*0
  ii <- lower.tri(nmat)
  
  fz <- .5*log( (1 + rmat[ii])/(1 - rmat[ii]) )
  se <- 1/sqrt(nmat[ii] - 3)
  ds <- se*qnorm(1 - alpha/2)
  lo[ii] <- .fishz2r(fz - ds)
  hi[ii] <- .fishz2r(fz + ds)
  sr[ (rmat < lo | rmat > hi) & ii] <- 1
  
  list(lo = lo, hi = hi, sig = sr)
}

.corPlot <-
function(cmat,slim=NULL,PDIAG=F,plotScale=1,
                    makeColor=NULL,textSize=NULL,
                    corLines=T,tri='lower',colorGrad = NULL,
                    cex=1,specLabs=T,squarePlot=T,LEGEND=T,
                    widex=10.5,widey=6.5,add=F,new=T){  #correlation or covariance matrix

  # makeColor - list of matrices of indices for boxes
  #   names of matrices are colors
  # if(PDIAG)diag(cmat) <- 0
  # tri - 'lower','upper', or 'both'
  # colorGrad - mapColors, heatColors, darkColors
  # squarePlot makes symbols square
  # new means NOT NEW 

  if(tri == 'upper')cmat[lower.tri(cmat)] <- 0
  if(tri == 'lower')cmat[upper.tri(cmat)] <- 0

  dy  <- nrow(cmat)
  dx  <- ncol(cmat)
  d <- dx
  xtext <- rep(c(1,100),dx/2)
  if(length(xtext) < d)xtext <- c(xtext,1)

  if(d < 20)xtext <- xtext*0 + 1

  xtext <- xtext*0 + 1

  ncol <- 200
  colF   <- colorRampPalette(c('darkblue','blue','green','white',
                               'yellow','red','brown'))
  colseq <- colF(ncol)

  if(is.null(slim))slim = range(cmat)
  slim  <- signif(slim,1)
  scale <- seq(slim[1],slim[2],length.out=ncol)
  
  if(slim[1] < 0 & slim[2] > 0){
    dp <- slim[2] - 0
    dm <- 0 - slim[1]
    ncol <- 200
    
    colseq <- colF(ncol)
    if(dp < dm)colseq <- colseq[101 + c(-100:round(dp/dm*100))]
    if(dp > dm)colseq <- colseq[ round((1 - dm/dp)*100):200 ]
    ncol  <- length(colseq)
  }

  ww   <- as.matrix(expand.grid(c(1:dy),c(1:dx)))  # note reverse order

  if(tri == 'upper'){
    ww  <- ww[ww[,1] <= ww[,2],]
    ww  <- ww[order(ww[,1]),]
  }
  if(tri == 'lower'){
    ww  <- ww[ww[,1] >= ww[,2],]
    ww  <- ww[order(ww[,1]),]
  }

  icol <- findInterval(cmat[ww],scale,all.inside=T)
  coli <- colseq[icol]

  if(PDIAG)coli[ww[,1] == ww[,2]] <- 'white'
  
  ss <- max(c(dx,dy))/5/plotScale
  
  if(squarePlot).mapSetup(c(0,dx),c(0,dy),ss,widex=widex,widey=widey)
  
  if(new)par(new = new)
  
  symbols(ww[,2],dy - ww[,1] + 1,squares=rep(1,nrow(ww)),
          xlim=c(0,dx+4),ylim=c(0,dy+4),
          fg=coli,bg=coli,inches=F,xlab=' ',ylab=' ',xaxt='n',yaxt='n',add=add)
  
  if(!is.null(makeColor)){
    
    for(k in 1:length(makeColor)){
      mm <- makeColor[[k]]
      if(length(mm) == 0)next
      if(tri == 'upper')mm <- mm[mm[,1] <= mm[,2],]
      if(tri == 'lower')mm <- mm[mm[,1] >= mm[,2],]
      ss <- matrix(0,dy,dx)
      ss[mm] <- 1
      wk <- which(ss[ww] == 1)
      ccc <- names(makeColor)[[k]]
      symbols(ww[wk,2],dy - ww[wk,1]+1,squares=rep(1,length(wk)),
              fg=ccc,bg=ccc,inches=F,xaxt='n',yaxt='n',add=T)
    }
  }

  if(corLines & tri == 'lower'){
    for(kk in 1:d){
      kb <- kk - .5
      ke <- d - kk + .5
      
      if(kk <= d)lines(c(kb,kb),c(0,ke),col='grey',lwd=1.5)          #verticle
      if(kk > 1) lines( c(.5,kb),c(ke,ke),col='grey',lwd=1.5)        #horizontal
      if(kk > 1) lines(c(kb,kb+.5),c(ke,ke+.5),col='grey',lwd=1.5)    #diagonal
    }
  }
  rect(0,-1,d+1,.5,col='white',border=NA)
  
  if(is.null(textSize))textSize <- exp(-.02*ncol(cmat))
  labels   <- rev(rownames(cmat))
  if(!specLabs)labels <- F
  
  if(tri == 'lower' & specLabs)text( c(d:1)+.1*xtext, c(1:d)+.5, rev(colnames(cmat)),pos=4,srt=45,cex=textSize)
  if(tri == 'both'){
    if(specLabs)text( c(dx:1)-.1*xtext, xtext*0+dy+.8, rev(colnames(cmat)),
                      pos=4,srt=55,cex=textSize)
    par(las = 1)
    axis(side=2,at=c(1:dy),labels=labels,tick=F,lwd=0,pos=.5,cex.axis=textSize)
    par(las = 0)
  }
  
  labside <- 'right'
  
  if(LEGEND).colorLegend(c(dx+1,dx+2),c(.5*dy,dy),ytick=c(slim[1],0,slim[2]),
                        scale,cols=colseq,labside=labside,
                        endLabels=range(slim),text.col='black')
}
.cov2Cor <-
function(covmat){  #covariance matrix to correlation matrix

  d    <- nrow(covmat)
  di   <- diag(covmat)
  s    <- matrix(di,d,d)
  covmat/sqrt(s*t(s))
}
.cov2Dist <-
function(sigma){ #distance induced by covariance
	
	n <- nrow(sigma)
	matrix(diag(sigma),n,n) + matrix(diag(sigma),n,n,byrow=T) - 2*sigma
}
.dMVN <-
function(xx,mu,smat=NULL,sinv=NULL,log=F){          #MVN density for mean 0
  
  if(!is.matrix(xx))xx <- matrix(xx,1)
  if(!is.matrix(mu))mu <- matrix(mu,1)
  
  xx <- xx - mu
  
  if(!is.null(sinv)){
    distval <- diag( xx%*%sinv%*%t(xx) )
    ev      <- eigen(sinv, only.values = TRUE)$values
    logd    <- -sum(log(ev))
  }
  
  if(is.null(sinv)){
    testv <- try(chol(smat),T)
    if(inherits(testv,'try-error')){
       tiny  <- min(abs(xx))/100 + 1e-5
       smat  <- smat + diag(diag(smat + tiny))
       testv <- try(chol(smat),T)
    }
    covm    <- chol2inv(testv)
    distval <- rowSums((xx %*% covm) * xx)
    ev      <- eigen(smat, only.values = TRUE)$values 
    logd    <- sum(log( ev ))
  }

    z <- -(ncol(xx) * log(2 * pi) + logd + distval)/2
    if(!log)z <- exp(z)
    z
}
.fishz2r <- function(z){ (exp(2*z) - 1)/(exp(2*z) + 1) }
.getScoreNorm <-
function(x,mu,xvar){  #Gneiting and Raftery's proper scoring rule

  #outcome x, prediction mean variance (mu, xvar)

  - ( (x - mu)^2)/xvar - log(xvar)

}
.gjamA2B <-
function(alpha,sigma){
  S <- nrow(sigma)
  q <- nrow(alpha)
  alpha/matrix( sqrt(diag(sigma)),q,S,byrow=T)
}
.gjamBaselineHist <-
function(y1,ylo = 0,nclass=20){
  
  # add histogram to base of current plot
  
  hh <- hist(y1,nclass=nclass,plot=F)
  dx <- diff(hh$breaks)[1]
  yo <- par('usr')[3]
  dy <- diff(par()$usr[3:4])
  hx <- hh$mids
  hy <- hh$density
  hy <- .3*hy*dy/max(hy) + yo
  nh <- length(hy)
  px <- hx[nh] + diff(hy)[nh-1]
  hx <- c(par()$xaxp[1],hx,rep(px,2))
  hy <- c(hy[1],hy,hy[nh],0)
  
  lines(hx,hy+yo,type='s',lwd=4,col='black')
  lines(hx,hy+yo,type='s',lwd=2,col='white')
}
.gjamCensorSetup <-
function(y,w,z,plo,phi,wm,censorMat){
  
  nc <- ncol(censorMat)
  br <- numeric(0)
  nk <- length(wm)
  n  <- nrow(y)
  
  zk <- y[,wm]*0
  blast <- -Inf
  
  for(j in 1:nc){
    
    valuej <- censorMat[1,j]
    bj     <- censorMat[2:3,j]
    names(bj) <- paste('c-',names(bj),j,sep='')
    
    if(j > 1){
      if(censorMat[2,j] < censorMat[3,j-1] )stop('censor intervals must be unique')
      if(bj[1] == br[length(br)])bj <- bj[2]
    }
    br <- c(br,bj)
    nb <- length(br)
    
    zk[ y[,wm] > blast & y[,wm] < bj[1] ] <- nb - 2
    zk[ y[,wm] == valuej ] <- nb - 1
    blast <- br[length(br)]
  }
  
  if(nc == 1){
    zk[zk == 0] <- 2
    br <- c(br,Inf)
  }
  
  zk[zk == 0] <- 1
  br <- matrix(br,nk,length(br),byrow=T)
  
  censk    <- which(y[,wm] %in% censorMat[1,])
  z[,wm]   <- zk
  
  tmp   <- .gjamGetCuts(z,wm)
  cutLo <- tmp$cutLo
  cutHi <- tmp$cutHi
  
  plo[,wm] <- br[cutLo]
  phi[,wm] <- br[cutHi]
  
  tmp <-  .tnorm(nk*n,plo[,wm],phi[,wm],w[,wm],1)
  
  w[,wm][censk] <- tmp[censk]
  
  imat <- w*0                    #location in full matrix
  imat[,wm][censk] <- 1
  censValue <- which(imat == 1)
  
  list(w = w, z = z, cutLo = cutLo, cutHi = cutHi, plo = plo, phi = phi, censValue = censValue,
       breakMat = br)
}
.gjamCuts2theta <-
function(tg,ss){   # variance to correlation scale
  nc   <- ncol(tg)
  sr   <- nrow(ss)
  tg/matrix( sqrt(diag(ss)),sr,nc)
}
.gjamGetCuts <-
function(zz,wk){
  
  nk <- length(wk)
  n  <- nrow(zz)

  cutLo <- cbind( rep(1:nk,each=n), as.vector(zz[,wk]) )
  cutHi <- cbind( rep(1:nk,each=n), as.vector(zz[,wk]) + 1 )
  
  list(cutLo = cutLo, cutHi = cutHi)
}
.gjamGetTypes <-
function(typeNames=NULL){
  
  TYPES <- c('PA','CA','DA','FC','CC','OC')
  FULL  <- c('presenceAbsence','continuous','discrete','fracComp','countComp','ordinal')
  
  if(is.null(typeNames)){
    names(FULL) <- TYPES
    return( list(typeCols = NULL, TYPES = TYPES, typeFull = FULL ) )
  }
  
  typeCols <- match(typeNames,TYPES)
  
  ww <- which(is.na(typeCols))
  if(length(ww) > 0)stop( paste('type code error',typeNames[ww],sep=' ') )
  
  list(typeCols = typeCols, TYPES = TYPES, typeFull = FULL[typeCols] )
}
.gjamHoldoutSetup <-
function(holdoutIndex,holdoutN,n){
  
  #holdout samples
  if(length(holdoutIndex) > 0)holdoutN <- length(holdoutIndex)
  if(holdoutN > (n/5))stop('too many holdouts')
  
  inSamples <- c(1:n)
  if(holdoutN > 0){
    if(length(holdoutIndex) == 0)holdoutIndex <- sort( sample(n,holdoutN) )
    inSamples <- inSamples[-holdoutIndex]
  }
  nIn <- length(inSamples)
  
  list(holdoutIndex = holdoutIndex, holdoutN = holdoutN, inSamples = inSamples, nIn = nIn)
}
.gjamMissingValues <-
function(x,y){
  
  # missing values in x
  xmiss  <- which(!is.finite(x),arr.ind=T)
  nmiss  <- length(xmiss)
  missX  <- missX2 <- xprior <- yprior <- numeric(0)
  
  xbound <- apply(x,2,range,na.rm=T)
  
  if(nmiss > 0){         #initialize missing values with means
    xmean    <- colMeans(x,na.rm=T)
    x[xmiss] <- xmean[xmiss[,2]]
    xprior   <- x[xmiss]
    warning( paste(nmiss,' missing values in x will be imputed',sep='') )
    nmiss <- nrow(xmiss)
    missX <- missX2 <- rep(0,nmiss)
  }
  
  # missing values in y
  ymiss <- which(!is.finite(y),arr.ind=T)
  mmiss <- length(ymiss)
  missY <- missY2 <- numeric(0)
  
  if(mmiss > 0){         #initialize missing values with means by TYPEs
    ymean    <- colMeans(y,na.rm=T)
    y[ymiss] <- ymean[ymiss[,2]]
    yprior   <- y[ymiss]
    warning( paste(mmiss,' missing values in y will be imputed',sep='') )
    mmiss <- nrow(ymiss)
    missY <- missY2 <- rep(0,mmiss)
  }
  
  list(xmiss = xmiss, xbound = xbound, missX = missX, missX2 = missX2,
       ymiss = ymiss, missY = missY, xprior = xprior, yprior = yprior)
}
.gjamPlotPars <-
function(type='CA',y1,yp,censm=NULL){
  
  if(!is.matrix(y1))y1 <- matrix(y1)
  if(!is.matrix(yp))yp <- matrix(yp)
  
  n       <- nrow(y1)
  nk      <- ncol(y1)
  nbin    <- NULL
  nPerBin <- n*nk/15
  breaks  <- NULL
  xlimit  <- range(y1,na.rm=T)
  ylimit  <- quantile(yp,c(.01,.99),na.rm=T)
  vlines  <- NULL
  wide    <- NULL
  MEDIAN  <- T
  LOG     <- F
  yss     <- sd(as.vector(y1))
  
  if(type == 'PA'){
    breaks  <- c(-.5,.5,1.5)
    wide    <- c(.04,.04)
    nPerBin <- NULL
    ylimit  <- xlimit <- c(0,1)
  } 
  if(type == 'OC'){
    breaks  <- seq(-.5,max(y1) + .5,by=1)
    wide    <- 1/max(y1)
    nPerBin <- NULL
    ylimit  <- range(yp,na.rm=T)
    xlimit  <- c( min(floor(y1)), max(ceiling(y1)) )
  } 
  if(type == 'DA')MEDIAN <- F
  if(type %in% c('DA','CA')){
    if(yss > 5){
      xlimit <- range(y1)
      LOG <- T
    }
  }
  if(type %in% c('FC','CC')){
    xlimit <- range(y1)
    MEDIAN <- F
    nPerBin <- round( n*nk/15,0 )
    if(type  == 'CC')LOG <- T
  } 
  if( !is.null(censm) ){
    
    cc  <- censm$partition
    vlines  <- numeric(0)
    breaks  <- NULL
    nPerBin <- n*nk/30
    xlimit  <- range(y1)
    ylimit  <- quantile(yp,c(.01,.99),na.rm=T)
    
    if(ncol(cc) > 1){
      cm     <- unique( as.vector(cc[-1,]) )
      vlines <- cm[is.finite(cm)]
      breaks <- vlines
      nbin   <- nPerBin <- NULL
      uncens <- cbind(cc[3,-ncol(cc)],cc[2,-1])
      wu     <- which( uncens[,1] != uncens[,2] )
      for(m in wu){
        sm <- seq(uncens[m,1],uncens[m,2],length=round(10/length(wu),0))
        if(type == 'DA') sm <- c(uncens[m,1]:uncens[m,2])
        breaks <- c(breaks,sm)
      }
      breaks <- c(breaks,max(y1) + 1)
      breaks <- sort( unique(breaks) )
    }
  }
  
  if(LOG){
    xlimit[1] <- ylimit[1] <- 1
    w0     <- which(y1 == 0)
    y1[w0] <- ylimit[1]
    w0     <- which(yp == 0)
    yp[w0] <- ylimit[1]
    nPerBin <- nPerBin/4
    ylimit[1] <- xlimit[1] <- 1
  }
  
  list( y1 = y1, yp = yp, nbin=nbin, nPerBin=nPerBin, vlines=vlines,
        xlimit=xlimit,ylimit=ylimit,breaks=breaks,wide=wide,LOG=LOG,
        POINTS=F,MEDIAN=MEDIAN )
}
.gjamPredictTraits <-
function(w,traitMat,traitTypes){
  
  M  <- nrow(traitMat)
  tn <- rownames(traitMat)
  
  ww <- w
  ww[ww < 0] <- 0
  
  tt <- ww%*%t(traitMat)
 # wf <- grep('FC',traitTypes)
 # if(length(wf) > 0){
 #   w0 <- which(tt[,wf] < 0)
 #   tt[tt[,wf] < 0,wf] <- 0
 #   tsum <- colSums(tt)
 #   tt   <- sweep(tt,1,tsum,'/')
 # }
  tt
}
.gjamSetup <-
function(typeCols,x,y,breakList=NULL,holdoutN, holdoutIndex,
                       censor=NULL,effort=NULL,
                       maxBreaks=100){
  
  Q <- ncol(x)
  n <- nrow(y)
  S <- ncol(y)
  allTypes <- unique(typeCols)
  typeNames <- names(typeCols)
  
  tmp <- .gjamGetTypes(typeNames)
  typeFull <- tmp$typeFull
  
  cuts <- cutLo <- cutHi <- numeric(0)
  minOrd <- maxOrd <- breakMat <- numeric(0)
  
  ordCols  <- which(names(typeCols) == 'OC')
  disCols  <- which(names(typeCols) == 'DA')
  compCols <- which(names(typeCols) == 'CC')
  
  wo <- grep('others',colnames(y))
  if(length(wo) > 0)colnames(y)[wo] <- 'other'
  
  other <- which(colnames(y) == 'other')
  
  colnames(y) <- .replaceString(colnames(y),now=' ',new='')
  colnames(x) <- .replaceString(colnames(x),now=' ',new='')
  
  w  <- y 
  z  <- w*0
  z[y == 0] <- 1
  z[y > 0]  <- 2
  plo <- phi <- y*0
  plo[z == 1] <- -Inf
  phi[z == 2] <- Inf
  
  censorCA <- censorDA <- numeric(0)      # indicates CA and DA values that will be sampled
  
  for(k in allTypes){
    
    wk <- which(typeCols == k)
    nk <- length(wk)
    
    if( typeFull[wk[1]] == 'presenceAbsence' ){       
      
      w[,wk] <- .tnorm(nk*n,plo[,wk],phi[,wk],0,1)
      br <- c(-Inf,0,Inf)
      br <- matrix(br,nk,length(br),byrow=T)
      colnames(br) <- as.character(c(1:ncol(br)))
      rownames(br) <- paste('PA',wk,sep='-')
      breakMat     <- .appendMatrix(breakMat,br,SORT=T,asNumbers=T)
    }
    
    if( typeFull[wk[1]] == 'continuous' ){       
      
      bb <- solve(crossprod(x))%*%crossprod(x,y[,wk])
      mu <- x%*%bb
      w0 <- which(y[,wk] == 0)
      w[,wk][w0] <- .tnorm(length(w0),plo[,wk][w0],phi[,wk][w0],mu[w0],1)
      
      br <- c(-Inf,0,Inf)
      br  <- matrix(br,nk,length(br),byrow=T)
      colnames(br) <- as.character(c(1:ncol(br)))
      rownames(br) <- paste('CA',wk,sep='-')
      
      if( !is.null(censor) & 'CA' %in% names(censor) ){
        
        wc     <- which(names(censor) == 'CA')
        bc     <- censorCA <- numeric(0)
        
        for(m in wc){
          
          wm     <- censor[[m]]$columns
          tmp    <- .gjamCensorSetup(y,w,z,plo,phi,wm,censorMat=censor[[m]]$partition)
          w[,wm] <- tmp$w[,wm]
          z[,wm] <- tmp$z[,wm]
          plo[,wm] <- tmp$plo[,wm]
          phi[,wm] <- tmp$phi[,wm]
          censorCA <- c(censorCA,tmp$censValue)
          bt       <- tmp$breakMat
          colnames(bt) <- as.character(c(1:ncol(bt)))
          rownames(bt) <- paste('CA',wm,sep='-')
          
          bc <- .appendMatrix(bc,bt,SORT=T,asNumbers=T)
        }
        mm <- match(rownames(bc),rownames(br))
        
        bb <- br[-mm,]
        tmp <- .appendMatrix(bc,bb,SORT=T,asNumbers=T)
        o   <- as.numeric( matrix( unlist(strsplit(rownames(tmp),'-')),ncol=2,byrow=T)[,2] ) 
        br <- tmp[order(o),]
      }
      breakMat <- .appendMatrix(breakMat,br,SORT=T,asNumbers=T)
    }
    
    if( typeFull[wk[1]] == 'discrete' ){
      
      if(is.null(effort)){      # if effort, w is on effort scale
        wide <- .5
      } else {
        we       <- wk[wk%in% effort$columns]
        w[,we]   <- y[,we]/effort$values
        wide     <- .5/effort$values
      }
      
      disCols <- wk
      
      lo <- w[,wk] - wide
      hi <- w[,wk] + wide
      lo[lo < 0] <- -Inf
      
      plo[,wk] <- lo
      phi[,wk] <- hi
      
      z[,wk] <- y[,wk] + 1
      w[,wk] <- .tnorm(nk*n,lo,hi,w[,wk],1)
      
      n <- nrow(y)
      S <- ncol(y)
      
      br <- c(-Inf,seq(0,(max(y[,wk])-1)),Inf)
      if(length(br) > maxBreaks){
        warning('breaks created')
        br <- c(br[1:maxBreaks],Inf)
      }
      br <- matrix(br,nk,length(br),byrow=T)
      colnames(br) <- as.character(c(1:ncol(br)))
      rownames(br) <- paste('DA',wk,sep='-')
      
      if( !is.null(censor) & 'DA' %in% names(censor) ){
        
        wc     <- which(names(censor) == 'DA')
        bc     <- censorDA <- numeric(0)
        
        for(m in wc){
          wm     <- censor[[m]]$columns
          tmp    <- .gjamCensorSetup(y,w,z,plo,phi,wm,censorMat=censor[[m]]$partition)
          w[,wm] <- tmp$w[,wm]
          z[,wm] <- tmp$z[,wm]
          plo[,wm] <- tmp$plo[,wm]
          phi[,wm] <- tmp$phi[,wm]
          censorDA <- c(censorDA,tmp$censValue)
          bt       <- tmp$breakMat
          colnames(bt) <- as.character(c(1:ncol(bt)))
          rownames(bt) <- paste('DA',wm,sep='-')
          
          bc <- .appendMatrix(bc,bt,SORT=T,asNumbers=T)
        }
        mm <- match(rownames(bc),rownames(br))
        
        bb <- br[-mm,]
        tmp <- .appendMatrix(bc,bb,SORT=T,asNumbers=T)
        o   <- as.numeric( matrix( unlist(strsplit(rownames(tmp),'-')),ncol=2,byrow=T)[,2] ) 
        br <- tmp[order(o),]
      }

      breakMat <- .appendMatrix(breakMat,br,SORT=T,asNumbers=T)
    }
      
    if(typeFull[wk[1]] == 'fracComp'){
      
      lo <- plo[,wk]
      hi <- phi[,wk]
      lo[lo < -1] <- -1
      hi[hi > 2]  <- 2
      plo[,wk] <- lo
      phi[,wk] <- hi
      
      w0 <- which(y[,wk] == 0)
      w[,wk][w0] <- .tnorm(length(w0),lo[w0],hi[w0],0,1)
      
      br <- c(-1,0,1)
      br <- matrix(br,nk,length(br),byrow=T)
      colnames(br) <- as.character(c(1:ncol(br)))
      rownames(br) <- paste('FC',wk,sep='-')
      
      breakMat <- .appendMatrix(breakMat,br,SORT=T,asNumbers=T)
    }
    
    if(typeFull[wk[1]] == 'countComp'){
      
      z[,wk] <- y[,wk] + 1
      ee     <- rowSums(y[,wk])  + 1
      
      lo <- (y[,wk] - .5)/ee          #zeros might not be zeros
      hi <- (y[,wk] + .5)/ee
      
      lo[lo < 0]  <- -2
      
      plo[,wk] <- lo
      phi[,wk] <- hi
      
      tmp <- matrix( .tnorm(nk*n,plo[,wk],phi[,wk],0,1),n,nk )
      
      tt <- tmp
      tt[tt < 0] <- 0
      tsum <- rowSums(tt)
      tt   <- sweep(tt,1,tsum,'/')
      tt[tmp < 0] <- tmp[tmp < 0]
      
      w[,wk] <- tt
      
      br <- c(-1,0,1)
      br <- matrix(br,nk,length(br),byrow=T)
      colnames(br) <- as.character(c(1:ncol(br)))
      rownames(br) <- paste('CC',wk,sep='-')
      
      breakMat <- .appendMatrix(breakMat,br,SORT=T,asNumbers=T)
    }
    
    if(typeFull[wk[1]] == 'ordinal'){
      
      ordCols <- wk
      nc <- apply(y[,wk],2,max)
      
      # more than one obs needed in last cell to estimate partition
      ii  <- list(spec = as.vector(matrix(c(1:nk),n,nk,byrow=T)), ss = as.vector(y[,wk]))
      ctmp <- .byIndex(as.vector(y[,wk])*0+1,ii,sum)
      
      ncc <- nc + 1
      if(max(ncc) > ncol(ctmp))ncc <- nc
      
      maxOne <- which(ctmp[ cbind(1:nk,ncc) ] == 1)
      
      if(length(maxOne) > 0){

        for(m in 1:length(maxOne)){
          mc <- wk[maxOne[m]]
          y[y[,mc] == nc[maxOne[m]],mc] <- nc[maxOne[m]] - 1
        }
        nc <- apply(y[,wk],2,max)
        message('note: single values in last ordinal class moved down one class')
      }
      
      ncut <- max(y[,wk])
      crow <- c(0:ncut)
      cuts <- t( matrix(crow,(ncut+1),nk) )
      cuts[ cbind((1:nk),nc+1) ] <- Inf
      
      call <- t( apply(cuts,1,cumsum) )
      cuts[call == Inf] <- Inf
      cuts <- cbind(-Inf,cuts)
      if(!is.matrix(cuts))cuts <- matrix(cuts,1)
    
      tmp   <- .gjamGetCuts(y + 1,wk)
      cutLo <- tmp$cutLo
      cutHi <- tmp$cutHi
      
      ss   <- seq(0,(nk-1)*n,by=n)
      wh <- as.vector( outer(holdoutIndex,ss,'+') )
      c1 <- cutLo
      if(length(wh) > 0)c1 <- cutLo[-wh,]

      otab   <- table(c1[,1],c1[,2])
      
      oo <- cbind(0,t( apply(otab,1,cumsum) ))
      wo <- which(oo == 0,arr.ind=T)
      minOrd <- .byIndex(wo[,2],wo[,1],max)
      
      oo <- cbind(0,t( apply( t(apply(otab,1,rev)),1,cumsum) ))
      wo <- which(oo == 0,arr.ind=T)
      maxOrd <- ncut - .byIndex(wo[,2],wo[,1],max) + 2
      
      
      plo[,wk] <- cuts[cutLo]
      phi[,wk] <- cuts[cutHi]
      
      z[,wk] <- y[,wk] + 1
      w[,wk] <- matrix( .tnorm(nk*n,plo[,wk],phi[,wk],y[,wk],1),n,nk )
      
      colnames(cuts) <- c(1:ncol(cuts))
      rownames(cuts) <- paste('OC',wk,sep='-')
      breakMat <- .appendMatrix(breakMat,cuts,SORT=T,asNumbers=T)
    }
  }
  
  ii <- list(spec = as.vector(matrix(c(1:S),n,S,byrow=T)), discrete_class = as.vector(z))
  classBySpec <- .byIndex(as.vector(z)*0+1,ii,sum)
  rownames(classBySpec) <- colnames(y)
  
  ncc <- min(20,ncol(classBySpec))
  nrr <- min(20,nrow(classBySpec))
 # print( classBySpec[1:nrr,1:ncc] )
  
    list(w = w, z = z, y = y, other = other, cuts = cuts, cutLo = cutLo, cutHi = cutHi, 
         plo = plo, phi = phi, ordCols=ordCols, disCols = disCols, compCols = compCols,
         classBySpec = classBySpec, breakMat = breakMat, minOrd = minOrd, maxOrd = maxOrd,
         censorCA = censorCA, censorDA = censorDA )
 }
.gjamTheta2cuts <-
function(tg,ss){
  nc   <- ncol(tg)
  sr    <- nrow(ss)
  tg/matrix( sqrt(diag(ss)),sr,nc)
}
.gjamTrueVest <-
function(chains,true,typeCode,allTypes,xlim=NULL,ylim=NULL,
                          label=NULL,colors=NULL,add=F,legend=T){
  
  true   <- as.vector(true)
  ntypes <- length(allTypes)
  
  if(is.null(ylim))ylim <- range(chains,na.rm=T)
  if(is.null(xlim))xlim <- range(true,na.rm=T)
  
  if(!is.matrix(chains)){
    chains <- matrix(chains,ncol=1)
    bCoeffTable <- c(mean(chains),sd(chains),quantile(chains,c(.025,.975)),true)
    bCoeffTable <- matrix(bCoeffTable,1)
  } else {
    bCoeffTable <- .processPars(chains,xtrue=true )
  }
  
  if(is.null(colors)){
    colors <- 1
    if(ntypes > 1)colors <- typeCode
  }
  if(length(colors) == 1) colors <- rep(colors,ntypes)
  
  .predVsObs(true,p=chains,xlab='true',xlim=xlim,ylim=ylim,ylab='estimated',
            colors=colors,add=add)
  
  if(ntypes > 1 & legend)legend('topleft',allTypes,text.col=colors,bty='n')
  if(!is.null(label)).plotLabel(label,above=T)
  
  bCoeffTable
}
.gjamUpdateBetaNoPrior <-
function(WIX,IXX,sg,...){
  
  bg <- matrix( .rMVN(1,as.vector(WIX),kronecker(sg,IXX)),nrow(IXX),ncol(WIX) )
  
  list(bg = bg, bzero = bg)
}

.conditionalMVNRcpp <- function(xx, mu, sigma, cdex, p=ncol(mu)){  
  # xx, mu are matrices
  
  # if(!length(xx) > nrow(sigma))return( conditionalMVN(xx,mu,sigma,cdex) )
  
  gdex <- (1:p)[-cdex] - 1
  cdex <- cdex - 1
  .conditionalMVNRcppCpp(cdex, gdex, xx, mu, sigma) 
}

.byRcpp <- function(x, i, j, summat=matrix(0,max(i),max(j)), 
                   totmat=summat, fun='mean'){  #
  
  nn <- length(x)
  if( nn != length(i) | nn != length(j) )stop('vectors unequal in byFunctionRcpp')
  if( nrow(summat) < max(i) | ncol(summat) < max(j) )stop('matrix too small')
  
  ww <- which(is.na(x))
  if(length(ww) > 0){
    x <- x[-ww]
    i <- i[-ww]
    j <- j[-ww]
  }
  
  frommat <- cbind(i,j,x)
  
  nr  <- nrow(frommat)
  maxmat <- summat*0 - Inf
  minmat <- summat*0 + Inf
  
  tmp <- .byRccpCpp(nr, frommat, totmat, summat, minmat, maxmat)
  
  if(fun == 'sum')return(tmp$sum)
  if(fun == 'mean'){
    mu <- tmp$sum/tmp$total
    mu[is.na(mu)] <- 0
    return(mu)
  }
  if(fun == 'min'){
    return( tmp$min )
  }
  tmp$max
}
.tnormMVNmatrixRcpp <- function(avec, muvec, smat, 
                                lo=matrix(-1000,nrow(muvec),ncol(muvec)), 
                                hi=matrix(1000,nrow(muvec),ncol(muvec)),
                                whichSample = c(1:nrow(smat))){
  
  #lo, hi must be same dimensions as muvec,avec
  
  lo[lo < -1000] <- -1000
  hi[hi > 1000]  <- 1000
  
  if(max(whichSample) > length(muvec))stop('whichSample outside length(muvec)')
  
  r <- avec
  a <- .tnormMVNmatrixRcppCpp(avec, muvec, smat, lo, hi, whichSample, 
                               0:(nrow(smat)-1))  
  r[,whichSample] <- a[,whichSample]
  r
}

.gjamUpdateBetaPrior <-
function(WIX,IXX,sg,alpha,loBeta,hiBeta,...){
  
  S <- ncol(WIX)
  Q <- nrow(WIX)
  tmp <- .tnormMVNmatrixRcpp(avec=alpha,muvec=WIX,smat=sg,
                      lo=matrix(loBeta,Q,S),hi=matrix(hiBeta,Q,S))
  
  tmp[!is.finite(tmp)] <- alpha[!is.finite(tmp)]
  list(bg = tmp, alpha = tmp, tau = NULL)
}
.gjamUpdateTheta <-
function(w,cutg,cutLo,cutHi,ordCols,
                                holdoutN,holdoutIndex,minOrd,maxOrd){
  
  word <- w[,ordCols]
  ncut <- ncol(cutg)
  nc   <- ncut - 1
  n    <- nrow(w)
  nk   <- length(ordCols)
  
  c1 <- cutLo[,1]
  c2 <- cutLo[,2]
  c3 <- cutHi[,1]
  c4 <- cutHi[,2]
  
  if(holdoutN > 0){
    word <- word[-holdoutIndex,]
    ss   <- seq(0,(nk-1)*n,by=n)
    wh <- as.vector( outer(holdoutIndex,ss,'+') )
    c1 <- c1[-wh]
    c2 <- c2[-wh]
    c3 <- c3[-wh]
    c4 <- c4[-wh]
  }
  
  cmin <- .byRcpp(as.vector(word),c1,c2,fun='min')
  cmax <- .byRcpp(as.vector(word),c1,c2,fun='max')
  
  cmin[!is.finite(cmin[,1]),1] <- -10
  cmin[,2] <- 0
  cmax[,1] <- 0
  cmax[cmax == -Inf] <- Inf
  
  tmp <- .interpRows(cmax,startIndex=minOrd+1,endIndex=maxOrd-1,
             INCREASING=T,minVal=0,maxVal=Inf,
             defaultValue=NULL,tinySlope=.001)
  
  cmax[!is.finite(cmax)] <- tmp[!is.finite(cmax)]
  
  ww <- which(!is.finite(cmin) & is.finite(cmax),arr.ind=T)
  if(length(ww) > 0){
    w0 <- ww
    w0[,2] <- w0[,2] - 1
    cmin[ww] <- runif(nrow(ww),cmax[w0],cmax[ww])
  }
  
  ww <- which(is.finite(cmax))
  
  clo <- cmax[,-nc]
  chi <- cmin[,-1]
  clo[,1] <- -1
  
  ww <- which(is.finite(clo))
  
 # wna <- which(chi < clo)
 # if(length(wna) > 0){
 #   chi[wna] <- clo[wna]
 #   print(chi - clo)
 # }
  
  chi[ww] <- .tnorm(length(ww),clo[ww],chi[ww],clo[ww],10)
  chi[,1] <- 0
  cmax <- cbind(-Inf,chi,Inf)
  
  cmax[,ncut] <- Inf
  cmax[ cbind(1:nk,maxOrd+1) ] <- Inf
  
  wmin <- which(minOrd > 1)
  if(length(wmin) > 0){
    for(j in wmin)cmax[j,2:c(minOrd[j]+1)] <- 0:(minOrd[j] - 1)
  }
  
  cmax
}
.gjamUpdateW <-
function(w,muw,x,y,sg,alpha,cutg,plo,phi,effort = effort,
         typeNames,holdoutN, holdoutIndex, censorCA, censorDA,
         notOther,pg=NULL,pgPrior=c(1000,50),breakMat = breakMat){
  
  Q <- ncol(x)
  n <- nrow(y)
  S    <- nrow(sg)
  
  tmp <- .gjamGetTypes(typeNames)
  typeFull <- tmp$typeFull
  typeCols <- tmp$typeCols
  allTypes <- unique(typeCols)
  
  miny <- breakMat[,2]              
  yPredict  <-  w*0
  
  if( 'PA' %in% typeNames | 'OC' %in% typeNames ){    #expanded w must be drawn on this scale
    rss  <- sqrt(diag(sg))
    wss  <- w/matrix(rss,n,S,byrow=T)
    css  <- .cov2Cor(sg)
    muss <- x%*%alpha
  } 
  
  for(k in allTypes){
    
    wk <- which(typeCols == k)
    nk <- length(wk)
    wo <- which(wk %in% notOther)
    wu <- which(typeCols[notOther] == k)
    yp <- wp <- matrix(w[,wk]*0,ncol=length(wk))
    
    if( typeFull[wk[1]] %in% c('presenceAbsence','ordinal') ) {
      wp[,wo] <- .gjamWLoop(wss[,notOther],muss[,notOther],css[notOther,notOther],
                            wu,plo[,notOther],phi[,notOther])[,wu]
      yp[,wo] <- .rMVN(n,muss[,notOther],css[notOther,notOther])[,wu]
    } else {
      wp[,wo] <- .gjamWLoop(w[,notOther],muw[,notOther],sg[notOther,notOther],
                            wu,plo[,notOther],phi[,notOther])[,wu]
      yp[,wo] <- .rMVN(n,muw[,notOther],sg[notOther,notOther])[,wu]
    }

    if( typeFull[wk[1]] == 'ordinal' ){
      for(s in 1:nk)yp[,s] <- findInterval(yp[,s],cutg[s,]) - 1
    }
    if( typeFull[wk[1]] == 'presenceAbsence' ){
      yp[yp > 0]  <- 1 
      yp[yp <= 0] <- 0
    }
    
    yp[yp < matrix(miny[wk],n,nk,byrow=T)] <- 0
    
    if( typeFull[wk[1]] == 'continuous' ){
      if(length(censorCA) > 0){
        wp[-censorCA] <- y[,wk][-censorCA]           # if censoring, then defined for all
      } else {
        wy <- which(y[,wk] > 0)
        wp[wy] <- y[,wk][wy]
      }
    }
        
    if( typeFull[wk[1]] == 'discrete' ){
      if(length(censorDA) > 0) wp[-censorDA] <- y[,wk][-censorDA]
      if(!is.null(effort))yp <- yp*effort$values
    }
    
    if( typeFull[wk[1]] %in% c('countComp', 'fracComp') ){
      
      if( typeFull[wk[1]] == 'fracComp'){  #fractional data do not impute positive values
        wy <- which(y[,wk[wo]] > 0)
        wp[wy] <- y[,wk[wo]][wy]
      }
      
      ww <- wp
      ww[ww < 0] <- 0              
      
      tmp <- .gjamCompW2Y(ww,pg,pgPrior,notOther=wo)
      wc  <- tmp$ww
      pg  <- tmp$pg
      
      wc[wp < 0] <- wp[wp < 0]
      wp <- wc
      
      yp[yp < 0] <- 0
      
      yp <- .gjamCompW2Y(yp,pg,notOther=wo)$ww
      
      if( typeFull[wk[1]] == 'fracComp')wp[ww > 0] <- wc[,wo][ww > 0]     #normalized pos with negative values
      
      if( typeFull[wk[1]] == 'countComp' ){
        ysum <- rowSums(y[,wk])
        yp   <- round( sweep(yp,1,ysum,'*'), 0)  # prediction conditioned on actual sample size
      }
    }  
    w[,wk] <- wp
    yPredict[,wk] <- yp
  }
  
  list(w = w, yp = yPredict,pg = pg )
}
.gjamWLoop <-
function(ws,mus,sgs,wkk,lo,hi){
  
  n <- nrow(lo)
  tiny <- .00001

  for(s in wkk){
    
    ls <- lo[,s]
    hs <- hi[,s]
    
    tmp <- .conditionalMVNRcpp(ws,mus,sgs,s)
    mu  <- tmp$mu
    vr  <- max(tmp$vr,tiny)
    
    tmp <- .tnorm(n,ls,hs,mu,sqrt(vr))
    
    wl  <- which(tmp == ls)
    if(length(wl) > 0) tmp[wl] <- ls[wl] + tiny*(ls[wl])
    
    wl  <- which(tmp == hs)
    if(length(wl) > 0) tmp[wl] <- hs[wl] - (1 - tiny)*hs[wl]
    
    ws[,s] <- tmp
  }
  ws
}
.gjamXY <-
  function(x,y,typeCode,standardX=NULL,facNames){
    
    Q <- ncol(x)
    S <- ncol(y)
    n <- nrow(x)
    
    colnames(y) <- .replaceString(colnames(y),':','x')  # ':' reserved for interactions
    colnames(y) <- .replaceString(colnames(y),' ','')  
    
    #check design
    checkInt <- range(x[,1])
    if(checkInt[1] != 1 | checkInt[2] != 1)stop( paste('x[,1] must be intercept (ones)') )
    colnames(x)[1] <- 'intercept'
    
    xnames <- colnames(x)
    snames <- colnames(y)
    predXcols <- 2:Q
    VIF <- isFactor <- factorList <-isNonLinX <- isInt <- intMat <- isSquare <- NULL
    
    if(Q > 2){
      tmp <- .checkDesign(x)
      if(tmp$rank < tmp$p)stop( 'x not full rank' )
      isFactor    <- tmp$isFactor
      VIF         <- tmp$VIF
      designTable <- tmp$designTable
      wi <- grep(':',isFactor)
      if(length(wi) > 0)isFactor <- isFactor[-wi]
      
      if(length(isFactor) > 0){
        factorList        <- vector('list',length(facNames))
        names(factorList) <- facNames
        for(g in 1:length(factorList)){
          factorList[[g]] <- isFactor[grep(facNames[g],isFactor)]
        }
      }  
      
      wx <- grep('2)',colnames(x))
      if(length(wx) > 0){
        mm <- unique(unlist(strsplit(colnames(x)[wx],'^2)',fixed=T)))
        mm <- unlist( strsplit(mm,'I(',fixed=T) )
        mm <- mm[nchar(mm) > 0]
        wx <- c( which(colnames(x)%in%mm),wx )
        isSquare <- wx
      }
      
      wx <- grep(':',colnames(x))
      if(length(wx) > 0){
        mm <- matrix(unlist(strsplit(colnames(x)[wx],':')),ncol=2,byrow=T)
        mat <- matrix( match(mm,colnames(x)), ncol=2)
        mat <- cbind(wx,mat)
        colnames(mat) <- c('int','main1','main2')
        wx <- c( which(colnames(x)%in%mm),wx )
        isInt <- wx
        intMat <- mat
      }
      
      if(!is.null(isSquare))isNonLinX <- isSquare
      if(!is.null(isInt))isNonLinX <- sort(unique( c(isNonLinX,isInt)))
    }
    
    sc  <- which( xnames %in% standardX )
    xmean <- colMeans(x[,sc])
    xsd   <- apply(x[,sc],2,sd)
    xscale <- rbind(xmean,xsd)
    if(!is.null(standardX)){
      x[,sc] <- (x[,sc] - matrix(xmean,n,length(sc),byrow=T))/matrix(xsd,n,length(sc),byrow=T)
    }
    
    if(is.null(snames))snames <- paste('S',1:S,sep='-')
    if(is.null(xnames))xnames <- paste('x',1:Q,sep='-')
    
    snames <- sub('_','-',snames)
    xnames <- sub('_','-',xnames)
    
    colnames(y) <- snames
    colnames(x) <- xnames
    
    list(x = x, y = y, snames = snames, xnames = xnames, predXcols = predXcols,
         isInt = isInt, intMat = intMat, isSquare  = isSquare, 
         factorList = factorList, isFactor = isFactor, isNonLinX = isNonLinX, 
         designTable = designTable, xscale = xscale)
  }
.gjamCompW2Y <-
function(ww,pg,pgPrior=c(1000,50),notOther=c(1:(ncol(ww)-1))){
  
  n  <- nrow(ww)
  W  <- rowSums(ww[,notOther])
  wh <- which(W > pg)
  
  if(length(wh) > 0){
    contract <- (1 - (1 - pg)^(W[wh]/pg))/W[wh]
    ww[wh,]  <- ww[wh,]*contract        
  }
  pg <- rbeta(1,pgPrior[1] + n - length(wh),pgPrior[2] + length(wh) )
  
  ww[,-notOther] <- 1 - rowSums(ww[,notOther])
  
  list(pg = pg, ww = ww )
}
.imputX_MVN <-
function(xx,yy,beta,wmiss,sinv,xprior=0,xbound=NULL,priorWT=1){
  
  # priorWT is inverse of variance
  
  wcol <- unique(wmiss[,2])
  S    <- nrow(sinv)
  Q    <- nrow(beta)
  
  if(is.null(xbound))xbound <- apply(xx,2,range,na.rm=T)
  
  for(j in wcol){
    
    wj <- wmiss[wmiss[,2] == j,]
    if(!is.matrix(wj))wj <- matrix(wj,1,2)
    wr <- wj[,1]
    xp <- xprior[wmiss[,2] == j]
    
    bj <- matrix(beta[j,],1)
    bn <- matrix(beta[-j,],Q - 1)
    
    xn <- matrix(xx[wr,-j],length(wr))
    z <- yy[wr,] - xn%*%bn
    datwt <- bj%*%sinv%*%t(bj)
    V     <- 1/( datwt + priorWT*datwt )
    v     <- z %*%sinv%*%t(bj) + xp*priorWT
    
    xx[wj] <- .tnorm(length(wr),xbound[1,j],xbound[2,j],v%*%V,sqrt(V))
  }
  xx
}
.interp <-
function(y,INCREASING=F,minVal=-Inf,maxVal=Inf,defaultValue=NULL,
                   tinySlope=NULL){  #interpolate vector x

  if(is.null(defaultValue))defaultValue <- NA

  tiny <- .00001
  if(!is.null(tinySlope))tiny <- tinySlope

  y[y < minVal] <- minVal
  y[y > maxVal] <- maxVal

  n  <- length(y)
  wi <- which(is.finite(y))

  if(length(wi) == 0)return(rep(defaultValue,n))
  if(length(wi) == 1)ss <- tiny

  xx  <- c(1:n)
  z  <- y

  if(wi[1] != 1) wi <- c(1,wi)
  if(max(wi) < n)wi <- c(wi,n)

  ss <- diff(z[wi])/diff(xx[wi])

  ss[is.na(ss)] <- 0

  if(length(ss) > 1){
    if(length(ss) > 2)ss[1] <- ss[2]
    ss[length(ss)] <- ss[length(ss)-1]
  }
  if(INCREASING)ss[ss < tiny] <- tiny

  if(is.na(y[1]))  z[1] <- z[wi[2]] - xx[wi[2]]*ss[1]
  if(z[1] < minVal)z[1] <- minVal
  if(z[1] > maxVal)z[1] <- maxVal

  for(k in 2:length(wi)){

     ki <- c(wi[k-1]:wi[k])
     yk <- z[wi[k-1]] + (xx[ki] - xx[wi[k-1]])*ss[k-1]
     yk[yk < minVal] <- minVal
     yk[yk > maxVal] <- maxVal
     z[ki] <- yk
  }
  z
}
.interpRows <-
function(x,startIndex=rep(1,nrow(x)),endIndex=rep(ncol(x),nrow(x)),
                       INCREASING=F,minVal=-Inf,maxVal=Inf,
                       defaultValue=NULL,tinySlope=.001){  
  #interpolate rows of x subject to increasing

  nn  <- nrow(x)
  p  <- ncol(x)
  xx <- c(1:p)

  if(length(minVal) == 1)minVal <- rep(minVal,nn)
  if(length(maxVal) == 1)maxVal <- rep(maxVal,nn)

  ni   <- rep(NA,nn)
  flag <- numeric(0)

  z <- x

  for(i in 1:nn){
    if(startIndex[i] == endIndex[i]){
      z[i,-startIndex[i]] <- NA
      next
    }
    z[i,startIndex[i]:endIndex[i]] <- .interp(x[i,startIndex[i]:endIndex[i]],
                                             INCREASING,minVal[i],maxVal[i],
                                             defaultValue,tinySlope)
  }
  
  z
}
.invMatZero <-
function(sgibbs,nsim=2000,knames,index=NULL,SS = length(knames)){   # return conditional independence
  
  if(is.null(index))index <- c(1:nrow(sgibbs))
  
  S1 <- sqrt(ncol(sgibbs))
  
  if(is.null(colnames(sgibbs))){
    cnames <- paste('S',c(1:S1),sep='')
    cnames <- outer( cnames,cnames,paste,sep='_')
    colnames(sgibbs) <- cnames
    knames <- cnames
  }
  
#  dnames <- as.vector( outer( knames,knames,paste,sep='_') )
 # gnames <- matrix( unlist(strsplit(colnames(sgibbs),'_')),ncol=2,byrow=T)
 # keepS  <- which(dnames %in% colnames(sgibbs))
  
 # mnames <- gnames[1:S1,1]
 # keep   <- which(mnames %in% knames)
  
 # SS <- length(keep)
  
  jj   <- sample(index,nsim,replace=T)
  sj   <- rj <- matrix(NA,nsim,S1*S1)
  
  for(j in 1:nsim){
    ss     <- matrix(sgibbs[jj[j],],S1,S1) 
    rr     <- .cov2Cor(ss)
    rj[j,] <- as.vector( rr )
    sj[j,] <- as.vector( solve( ss ) )
  }
  
  tmp    <- apply(rj,2,quantile,c(.5,.025,.975))
  loMar  <- which(tmp[2,] < 0 & tmp[3,] < 0)
  hiMar  <- which(tmp[2,] > 0 & tmp[3,] > 0)
  inMar  <- which(tmp[2,] < 0 & tmp[3,] > 0)
  
  tmp    <- apply(sj,2,quantile,c(.5,.025,.975))
  loCon  <- which(tmp[2,] < 0 & tmp[3,] < 0)
  hiCon  <- which(tmp[2,] > 0 & tmp[3,] > 0)
  inCon  <- which(tmp[2,] < 0 & tmp[3,] > 0)
  
  ss <- matrix(0,S1,S1)
  ss[inCon] <- 1
#  ss[upper.tri(ss)] <- 0
  inConMat <- which(ss == 1,arr.ind=T)
  ss <- matrix(0,S1,S1)
  ss[inMar] <- 1
 # ss[upper.tri(ss)] <- 0
  inMarMat <- which(ss == 1,arr.ind=T)
  
  list( inMar = inMar, inCon = inCon, inMarMat = inMarMat, inConMat = inConMat )
}
.joinCharVec <-
function(charVec,sep=''){
  
  if(length(charVec) == 1) return(charVec)
  xx <- charVec[1]
  for(j in 2:length(charVec))xx <- paste(xx,charVec[j],sep=sep)
  xx
}
.mapSetup <-
function(xlim,ylim,scale=NULL,widex=10.5,widey=6.5){  
  
  #scale is x per inch
  
  if(is.null(scale))scale <- 1

  px   <- diff(xlim)/scale
  py   <- diff(ylim)/scale
  
  if(px > widex){
    dx <- widex/px
    px <- widex
    py <- py*dx
  }
  if(py > widey){
    dx <- widey/py
    py <- widey
    px <- px*dx
  }
    
  par(pin=c(px,py))
  
  invisible( c(px,py) )

}
.multivarChainNames <-
function(rowNames,colNames){
  as.vector( t(outer(colNames,rowNames,paste,sep='_')) )
}
.multivarEmat <-
function(bchains,covx,snames,orderB){
  
  #SB - columns in bchains, excludes other
  
  Q <- ncol(covx)
  SM <- length(snames)
  SB <- length(orderB)
  
  onames <- snames[orderB]
  o1     <- paste(onames,'_',sep='')
  o2     <- paste('_',onames,sep='')
  
  bcols <- ccols <- numeric(0)
  ccols <- character(0)
  for(j in 1:SB){
    wj    <- grep(o1[j],colnames(bchains))
    if(length(wj) == 0)wj    <- grep(o2[2],colnames(bchains))
    if(length(wj) == 0)next
    bcols <- c(bcols,min( wj ) )
    ccols <- c(ccols,colnames(bchains)[wj])
  }
  
  cnames <- colnames(bchains[,bcols])
  cnames <- .replaceString(cnames,now='intercept_',new='')
  cnames <- .replaceString(cnames,now='_intercept',new='')
  
  nsim  <- 2000
  ng    <- nrow(bchains)  
  index <- c(round(ng/10,0):ng)
  
  jj <- sample(index,nsim,replace=T)
  rj <- matrix(NA,nsim,SB*SB)
  
  for(j in 1:nsim){
    
    bb     <- matrix( bchains[jj[j],ccols],Q,SB)  #includes omitSpec, but not other
    ss     <- t(bb)%*%covx%*%bb
    rr     <- .cov2Cor(ss)
    rj[j,] <- as.vector( rr )
  }
  
  tmp    <- apply(rj,2,quantile,c(.5,.025,.975))
  bm     <- matrix(tmp[1,],SB,SB)
  rownames(bm) <- colnames(bm) <- cnames
  loMar  <- which(tmp[2,] < 0 & tmp[3,] < 0)
  hiMar  <- which(tmp[2,] > 0 & tmp[3,] > 0)
  inMar  <- which(tmp[2,] < 0 & tmp[3,] > 0)
  
  ss <- matrix(0,SB,SB)
  ss[inMar] <- 1
  whichZero <- which(ss == 1,arr.ind=T)
  
  list(bm = bm, whichZero = whichZero)
}
.rMVN <-
function (nn, mu, sigma){
  
  # nn - no. samples from one mu vector or nrow(mu) for matrix
  
  if(!is.matrix(mu)) mu <- matrix(mu,1)
  if(length(mu) == 1)mu <- matrix(mu,1,nrow(sigma))
  if(ncol(mu) == 1)  mu <- t(mu)
  
  m <- ncol(sigma)
  
  if(ncol(mu) != m)stop('dimension mismatch mu, sigma')
  
  if(nn > 1 & nrow(mu) == 1)mu <- matrix(mu,nn,m,byrow=T)
  
  if(nn != nrow(mu))stop('sample size does not match mu')
  
  testv <- try(svd(sigma),T)
  
  if( inherits(testv,'try-error') ){
    ev <- eigen(sigma, symmetric = TRUE)
    testv <- t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
  } else {
    testv <- t(testv$v %*% (t(testv$u) * sqrt(testv$d)))
  }
  
  p <- matrix(rnorm(nn * m), nn) %*% testv
  p + mu
}
.omitChainCol <-
function(cmat,omitCols){
  
  #omitCols - characterVector
  
  keep <- c(1:ncol(cmat))
  ocol <- numeric(0)
  for(j in 1:length(omitCols)){
    ocol <- c(ocol,grep(omitCols[j],colnames(cmat)))
  }
  if(length(ocol) > 0)keep <- keep[-ocol]
  list(keep = keep, omit = ocol)
}
.plotChainDensity <-
function(cMat,vnames=NULL,ncolPlot=2,xlab=' ',ylab=' ',
                             LOG=F,inColor=F){
  
  if(is.null(vnames)){
    mp <- 1
  } else {
    mp  <- length(vnames)/ncolPlot
  }
  nrr <- floor(mp) + 1
  
  graphics.off()
  par(mfrow=c(nrr,ncolPlot),bty='n', 
      oma=c(1,1,0,0), mar=c(1,1,1,0), mgp=c(0,0,0))
  
  knames <- vnames
  if(is.null(vnames))knames <- 'all'
  
  for(k in 1:length(knames)){
    
    kk <- knames[k]
    if(is.null(vnames))kk <- NULL
    
    if(!is.null(kk)){
      gg <- grep(kk,colnames(cMat))
      if(length(gg) == 0)next
    }
    
    tmp <- .chains2density(cMat,varName=kk)
    xt  <- tmp$x
    yt  <- tmp$y
    chainMat <- tmp$chainMat
    
    cols <- rep('black',nrow(xt))
    if(inColor){
      colF   <- colorRampPalette(c('black','brown','orange'))
      cols <- colF(nrow(xt))
    }
    
    nn <- nrow(chainMat)
    
    xrange <- quantile(xt,c(.05,.95))
    xlim <- range(xt)
    if(!LOG) plot(10,10,xlim=xlim,ylim=c(0,1.8*max(yt)),xlab=xlab,ylab=ylab,cex=.01)
    if(LOG)  plot(10,10,xlim=xlim,ylim=c(0,1.8*max(yt)),xlab=xlab,ylab=ylab,cex=.01,log='x')
    
    .plotLabel(kk,above=T)
    
    j1 <- 1
    if(knames[1] == 'intercept')j1 <- 2
    
    for(j in j1:nrow(xt)){
      
      cj <- cumsum(yt[j,])
      cj <- cj/max(cj)
      ci <- xt[j, findInterval(c(.02,.98),cj) ]
      
      label <- rownames(xt)[j]
      
      wm <- which.max(yt[j,])
      lines(xt[j,],yt[j,],lwd=2)
      lines(range(xt[j,]),c(0,0),col=cols[j],lwd=2)
      if(ncol(chainMat) < 25)text(xt[j,wm],1.1*yt[j,wm],label,srt=55,pos=4,cex=1)
    }
  }
}
.plotLabel <-
function(label,location='topleft',cex=1.3,font=1,above=F,below=F,bg=NULL){
  
  if(above){
    adj <- 0
    if(location == 'topright')adj=1
    title(label,adj=adj, font.main = font, font.lab =font)
    return()
  }
  if(below){
    adj <- 0
    if(location == 'bottomright')adj=1
    mtext(label,side=1,adj=adj, outer=F,font.main = font, font.lab =font,cex=cex)
    return()
  }
    
  if(is.null(bg)){
    tmp <- legend(location,legend=' ',bty='n')
  } else {
    tmp <- legend(location,legend=label,bg=bg,border=bg,text.col=bg,bty='o')
  }
  
  xt <- tmp$rect$left # + tmp$rect$w
  yt <- tmp$text$y
  
  pos <- 4
  tmp <- grep('right',location)
  if(length(tmp) > 0)pos <- 2
  
  XX <- par()$xlog
  YY <- par()$ylog
  
  if(XX)xt <- 10^xt
  if(YY)yt <- 10^yt
  
  text(xt,yt,label,cex=cex,font,pos=pos)
  
}
.plotObsPred <-
function(obs,yMean,ySE=NULL,nbin=NULL,nPerBin=NULL,breaks=NULL,
                        LOG=F,xlimit=NULL,ylimit=NULL,xlabel='Observed',ylabel='Predicted',
                        ptcol=NULL,boxPerc = .6826895, whiskerPerc = .95,
                        fill=NULL,add=F,box.col='black',wide=NULL,POINTS=T,
                        MEDIAN=T){
  
  aa <- (1 - boxPerc)/2
  boxQuant <- c(aa, 1 - aa )
  aa <- (1 - whiskerPerc)/2
  whiskerQuant <- c(aa, 1 - aa )
  
  if(is.null(ptcol)){
    ptcol <- 'black'
    ptcol <- 'grey'
    if(!is.null(nbin))ptcol <- 'grey'
  }
  if(length(ptcol) == 1)ptcol <- rep(ptcol,length(obs))
  
  if(is.null(xlimit))xlimit <- range(obs[is.finite(obs)],na.rm=T)
  if(is.null(ylimit))ylimit <- range(yMean[is.finite(yMean)],na.rm=T)
  
  xxx <- obs
  yyy <- yMean
  
  if(LOG){
    if(is.null(xlimit))xlimit <- range( obs[obs > 0],na.rm=T )
    if(is.null(ylimit))ylimit <- range( yMean[yMean > 0],na.rm=T )
    if(xlimit[1] <= 0)xlimit[1] <- .001
  }
  
  if(!POINTS){
    xxx <- xlimit[1]
    yyy <- ylimit[1]
  }
  
  if(!add){
    if(is.null(ylimit)){
      if(!LOG)plot(xxx,yyy,col=ptcol,cex=.3,xlab=xlabel,ylab=ylabel)
      if(LOG) plot(xxx,yyy,col=ptcol,cex=.3,xlab=xlabel,ylab=ylabel,log='xy')
    }
    if(!is.null(ylimit)){
      if(!LOG)plot(xxx,yyy,col=ptcol,cex=.3,xlab=xlabel,ylab=ylabel,xlim=xlimit,ylim=ylimit)
      if(LOG) plot(xxx,yyy,col=ptcol,cex=.3,xlab=xlabel,ylab=ylabel,xlim=xlimit,log='xy',ylim=ylimit)
    }
  }
  if(!is.null(ySE)){
    ylo <- yMean - 1.96*ySE
    yhi <- yMean + 1.96*ySE
    for(i in 1:length(obs))lines(c(obs[i],obs[i]),c(ylo[i],yhi[i]),col='grey',lwd=2)
  }
    
  if(is.null(breaks)){
    
    if(is.null(nbin))nbin <- 20
    
    br   <- range(obs[is.finite(obs)],na.rm=T)
    bins <- seq(br[1],br[2],length=nbin)
    if(LOG){
      oo <- min( obs[obs > 0],na.rm=T )
      br[1] <- .5*oo
      bins <- 10^seq(log10(br[1]),log10(br[2]),length=nbin)
    }
    if(!is.null(nPerBin)){
      nbb <- nPerBin/length(obs)
      nbb <- seq(0,1,by=nbb)
      if(max(nbb) < 1)nbb <- c(nbb,1)
      bins <- quantile(obs,nbb,na.rm=T)
      bins <- sort(unique(bins))
      nbin <- length(bins)
    }
  } else {
    bins <- breaks
    nbin <- length(bins)
  }
    
    if(is.null(wide))wide <- diff(bins)/2
    if(length(wide) == 1)wide <- rep(wide,nbin)
  
    for(k in 1:(nbin-1)){
      
      ok <- which(obs >= bins[k] & obs <= bins[k+1])
      qk <- which(is.finite(yMean) & obs >= bins[k] & obs <= bins[k+1])
      q  <- quantile(yMean[qk],c(.5,whiskerQuant[1],boxQuant[1],
                                 boxQuant[2],whiskerQuant[2]),na.rm=T)
      if(LOG){
        q[q <= 0] <- .001
      }
      ym <- mean(yMean[qk])
      xx <- mean(bins[k:(k+1)])
      if(!is.null(nPerBin) & MEDIAN)xx <- median(obs[ok],na.rm=T)
      points(xx,q[1],pch=3,col=box.col)
      yy <- q[c(2,5)]
      yy[1] <- max(yy[1],ylimit[1],na.rm=T) + .0000001
      yy[2] <- max(yy)
      lines(c(xx,xx),yy,lwd=2,col=box.col)
   
      yy1 <- q[3]
      yy1 <- max(yy1,ylimit[1],na.rm=T) + .00000001
      yy2 <- max(yy1,q[4])
      if(is.null(nPerBin)){
        
        minmax <- par('usr')[1:2]
        minx   <- max(c(minmax[1],xx-.5*wide[k])) + .00001
        maxx   <- min(c(minmax[2],xx+.5*wide[k])) - .00001
        rect(minx,yy1,maxx,yy2,col=fill,border=box.col)
        lines(c(minx,maxx),c(ym,ym),lwd=2,col=box.col)
      }
      if(!is.null(nPerBin)){
        qo <- quantile(obs[ok],c(.3,.7,.25,.75),na.rm=T)
        if(!MEDIAN)qo <- c(xx-.2*wide[k],xx+.2*wide[k],xx-.3*wide[k],xx+.3*wide[k])
        rect(qo[1],yy1,qo[2],yy2,col=fill,border=box.col)
        lines(c(qo[3],qo[4]),c(ym,ym),lwd=2,col=box.col)
      }
    }

  invisible( bins )
}
.predictY2X_linear <-
function(xx,yy,bb,ss,priorIV = diag(1e-10,ncol(xx)), 
                              priorX=matrix(0,ncol(xx)),predCols=c(2:ncol(xx))){
  
  #inverse prediction for multivariate linear in x
  
  priorX <- priorX[predCols]
  if(!is.matrix(priorX))priorX <- matrix(priorX)
  
  nn <- nrow(yy)
  notPred <- c(1:ncol(xx))[-predCols]
  
  bn <- matrix(bb[notPred,],length(notPred))
  bp <- matrix(bb[predCols,],length(predCols))
  
  yi <- yy - xx[,notPred]%*%bn
  pp <- length(predCols)
  
  bs <- bp%*%solve(ss)
  
  V <- solve( bs%*%t(bp) + priorIV[predCols,predCols] )
  v <- yi%*%t(bs) + matrix( priorIV[predCols,predCols] %*% priorX,nn,pp,byrow=T)
  mu <- v%*%V
  if(ncol(mu) > 1) xn <- .rMVN(nn,mu,V)
  if(ncol(mu) == 1)xn <- rnorm(nn,mu,sqrt(V))
  xx[,predCols] <- xn
  xx
}
.predictY2X_nonLinear <-
function(xx,yy,bb,ss,priorIV = diag(1e-10,ncol(xx)), 
                                 priorX=matrix(0,ncol(xx)),
                                 predCols=c(2:ncol(xx)),isInt=NULL,intMat=NULL,
                                 isFactor=NULL,factorList=NULL){
  
  #inverse prediction for multivariate nonlinear in x and factors, metropolis
  
  iFcol  <- NULL
  priorX <- priorX[predCols]
  if(!is.matrix(priorX))priorX <- matrix(priorX)
  
  nn <- nrow(yy)
  intercept <- xx[,1]
  
  xnew <- xx
  xnew[,predCols] <- .rMVN(nn,xx[,predCols],diag(.01,length(predCols)))
  
  if(!is.null(isFactor)){          # all factors, main effects
    xnew[,isFactor] <- 0
    for(k in 1:length(factorList)){
      nf <- length(factorList[[k]]) + 1
      wk <- sample(1:nf,nn,replace=T)
      tm <- cbind( 1:nn, wk )
      xnew[,c('intercept',factorList[[k]])][tm] <- 1
    }
    iFcol <- match(isFactor,colnames(xx))
  }
  
  if(length(intMat) > 0){                     #are some of the nlin terms interactions?
    pindex <- unique(as.vector(intMat[,-1]))
    if(length(iFcol) > 0)pindex <- pindex[!pindex %in% iFcol]  # those that are not factors
    if(length(pindex) > 0) xnew[,pindex] <- .rMVN(nn,xx[,pindex],diag(.01,length(pindex)))
    xnew[,intMat[,1]] <- xnew[,intMat[,2]]*xnew[,intMat[,3]]
  }
  
  pnow <- .dMVN(yy,xx%*%bb,smat=ss)
  pnew <- .dMVN(yy,xnew%*%bb,smat=ss)
  
  a  <- exp(pnew - pnow)
  z  <- runif(nn,0,1)
  wa <- which(z < a)
  xx[wa,] <- xnew[wa,]
  
  list(x = xx, accept = length(wa))
}
.predVsObs <-
function(true,p,xlim=range(true),ylim=range(p,na.rm=T),xlab=' ',ylab=' ',
                      colors=rep(1,length(true)),lwd=2,add=F){ 
	
  #true  - length n vector of obs or true values
  #p - ng by n matrix of estimates
  
  if(!is.matrix(p))p <- matrix(p,ncol=1)
  
  nn <- length(true)
  y  <- apply(p,2,quantile,c(.5,.025,.975))
  
  if(!add)plot(true,y[1,],xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,col=colors,pch=3,lwd=lwd)
  points(true,y[1,],col=colors,pch=3,lwd=lwd)

  for(j in 1:nn)lines(c(true[j],true[j]),y[2:3,j],col=colors[j],lwd=lwd)
  abline(0,1,lty=2)
  
  invisible(y)
}
.processPars <-
function(xgb,xtrue=numeric(0),CPLOT=F,DPLOT=F,
                        sigOnly = F,burnin=1,xlimits = NULL){  

  #xg      - matrix of gibbs chains
  #xtrue   - true values (simulated data)
  #CPLOT   - if T, plot chains
  #DPLOT   - if T, plot density
  #burnin  - analyze chains > burnin
  #xlimits - xlimits for plot
  #sigOnly - plot only parameters that 95% CI does not include 0
  
  if(!is.matrix(xgb))xgb <- matrix(xgb,ncol=1)
  if(is.null(colnames(xgb)))colnames(xgb) <- paste('V',c(1:ncol(xgb)),sep='-')
  
  NOPARS <- F
  
  if(sigOnly){
    wi   <- grep('intercept',colnames(xgb))      #extract covariates for plotting
    btmp <- xgb
    if(length(wi) > 0){
    	btmp <- xgb[,-wi]
      if(length(xtrue) > 0)xtrue <- xtrue[-wi]
    }

    wq   <- apply(btmp,2,quantile,c(.025,.975),na.rm=T)  #extract parameters != 0
    wq   <- which(wq[1,] < 0 & wq[2,] > 0)
    
    if(length(wq) == ncol(btmp))NOPARS <- T
    if(NOPARS) warning('no significant pars to plot')
    if(length(wq) > 0 & !NOPARS){
      xgb  <- btmp[,-wq]
      if(length(xtrue) > 0)xtrue <- xtrue[-wq]
    }
   }

  if(!is.matrix(xgb))xgb <- as.matrix(xgb)
  if(burnin > 1){
  	     if(burnin > (nrow(xgb) + 100))stop("burnin too large")
  	     xgb <- xgb[-c(1:burnin),]
  }
  if(!is.matrix(xgb))xgb <- as.matrix(xgb)
  nc <- ncol(xgb)
  nf <- round(sqrt(nc),0)

  out <- t(rbind(apply(xgb,2,mean,na.rm=T),apply(xgb,2,sd,na.rm=T),
                 apply(xgb,2,quantile,c(.025,.975),na.rm=T)))
  if(!is.null(colnames(xgb)))rownames(out) <- colnames(xgb)
  colnames(out) <- c('estimate','se','0.025','0.975')
  if(length(xtrue) > 0){
    out <- cbind(out,xtrue)
    colnames(out) <- c('estimate','se','0.025','0.975','true value')
  }

  if(CPLOT | DPLOT)par(mfrow=c((nf+1),nf),mar=c(4,2,2,2))
  if(CPLOT & DPLOT)par(mfrow=c((nf+1),nc),mar=c(4,2,2,2))

  if(CPLOT & !NOPARS){
      for(j in 1:nc){
       plot(xgb[,j],type='l')
       abline(h=out[j,],lty=2)
       if(length(xtrue) > 0)abline(h=xtrue[j],col='red')
       abline(h = 0, col='grey',lwd=2)
       title(colnames(xgb)[j])
     }
  }
  xlims <- xlimits
  if(DPLOT & !NOPARS){
      for(j in 1:nc){
        xj <- density(xgb[,j])
        if(is.null(xlimits))xlims <- range(xj$x)
        plot(xj$x,xj$y,type='l',xlim=xlims)
        abline(v=out[j,],lty=2)
        if(length(xtrue) > 0)abline(v=xtrue[j],col='red')
        title(colnames(xgb)[j])
     }
  }
  list(summary = signif(out,4)
)

}
.replaceString <-
function(xx,now='_',new=' '){  #replace now string in vector with new
  
  ww <- grep(now,xx)
  if(length(ww) == 0)return(xx)
  
  for(k in ww){
    s  <- unlist( strsplit(xx[k],now) )
    ss <- s[1]
    if(length(s) == 1)ss <- paste( ss,new,sep='')
    if(length(s) > 1)for(kk in 2:length(s)) ss <- paste( ss,s[kk],sep=new)
    xx[k] <- ss
  }
  xx
}
.tnorm <-
function(n,lo,hi,mu,sig){   

  #normal truncated lo and hi

  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))

  q1 <- pnorm(lo,mu,sig)
  q2 <- pnorm(hi,mu,sig) 

  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == Inf]  <- lo[z == Inf]
  z[z == -Inf] <- hi[z == -Inf]
  z
}
.traitLabel <-
function(tname){
  
  tname <- .replaceString(tname,now='soilFactor',new='')
  tname[tname == 'gmPerSeed'] <- 'Seed mass'
  tname[tname == 'gmPerCm']   <- 'Wood dens'
  tname[tname == 'woodSG']    <- 'Wood dens (green)'
  tname[tname == 'maxHt']     <- 'Max ht'
  tname[tname == 'leafN']     <- 'leaf [N]'
  tname[tname == 'leafP']     <- 'leaf [P]'
  tname[tname == "other"]  <- 'Deciduous'
  tname[tname == "broaddeciduous"]  <- 'Deciduous'
  tname[tname == "broadevergreen"]  <- 'BL evergrn'
  tname[tname == "needleevergreen"] <- 'NL evergrn'
  tname[tname == "dioecious"] <- 'Dioecious'
  tname[tname == "u1"] <- 'Slope'
  tname[tname == "u2"] <- 'Aspect 1'
  tname[tname == "u3"] <- 'Aspect 2'
  tname[tname == "ringPorous"] <- 'RP xylem'
  tname[tname == "temp"] <- 'Winter temperature'
  tname[tname == "stdage"] <- 'Stand age'
  for(j in length(tname)){
    tname[j] <- paste(toupper(substring(tname[j], 1, 1)), substring(tname[j], 2),sep = "", collapse = " ")
  }
  tname
}
.updateWishartNoPrior <-
function(xx,yy,df,beta=NULL,IXX=NULL,WX=NULL,WIX=NULL,TRYPRIOR=F){
  
  # more stable without prior
  # TRYPRIOR includes non-informative prior if cholesky fails 
  
  index <- 0
  
  if(is.null(IXX)) IXX <- solve( crossprod(xx) )
  if(is.null(WX))  WX  <- crossprod(xx,yy)
  if(is.null(WIX)) WIX <- IXX%*%WX
  
  SS   <- crossprod(yy) - t(WX)%*%WIX
  testv <- try(chol(SS),T)
  
  if( inherits(testv,'try-error') ){

    message('warning: updateWishartNoPrior')
 #   df <- nrow(SS) + nrow(xx)
    SS <- crossprod(yy - xx%*%beta) +  diag(diag(SS)*.001)*nrow(SS)
    testv <- try(chol(SS),T)
  }
  
  SI <- chol2inv(testv)
  
  testChol <- try(chol(SI),T)
    
  if( inherits(testChol,'try-error') ){
    message('warning: prior used in updateWishartNoPrior')
    if(TRYPRIOR){
      index  <- 1
      SI     <- SI + diag(diag(SI)*.01)
      df     <- nrow(SI) + nrow(xx)
      testChol <- try(chol(SI),T)
    }
  }

  z     <- matrix(rnorm(df*nrow(SS)),df,nrow(SS))%*%testChol
  sinv  <- crossprod(z)
  
  testSolve <- try( solve(sinv),T )
  if( !inherits(testSolve,'try-error') )sigma <- testSolve
  
  if( inherits(testSolve,'try-error') ){
    message('warning: prior used in updateWishartNoPrior')
    if(TRYPRIOR){
      sinv   <- sinv + diag(diag(sinv)*.001)
      df     <- nrow(sinv) + nrow(xx)
      testSolve <- try(chol(sinv),T)
      sigma   <- chol2inv(testSolve)
    }
  }
  
  list( sigma = sigma, sinv = sinv, indicator = index )
}
.w2z <-
function(w,ss){
  
  rr <- sqrt(diag(ss))
  w/matrix(rr,nrow(w),ncol(w),byrow=T)
}

.yaxisHorizLabs <-
function( labels, at=c(1:length(labels)), xshift=.05 ){
  
  #add horizontal y axis labels to existing plot
  
  text(par('usr')[3] - xshift*par('usr')[4] - par('usr')[3],y=at,labels,xpd=T)
}
