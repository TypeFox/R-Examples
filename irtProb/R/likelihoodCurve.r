`likelihoodCurve` <-
function(x, s, b, c, d,
         limitT=c(min=-4,max=4), limitS=c(min=0,max=4),
         limitC=c(min=0,max=1), limitD=c(min=0,max=1),
         grain=150, annotate=TRUE, logLikelihood=FALSE, color=TRUE,
         main="Likelihood Curve", xlab=expression(theta), ylab=NULL, zlab="P(X)",
         type="levelplot", m=0) {
 # Verification of the coherence of the item parameters
 if (length(b) <= 1) stop("The length of the item parameters vector must be greater than 1")
 if (length(s) == 1 && length(s) != length(b)) stop("Item parameters not well defined")
 if (length(c) == 1 && length(c) != length(b)) stop("Item parameters not well defined")
 if (length(d) == 1 && length(d) != length(b)) stop("Item parameters not well defined")
 nItems   <- max(length(s),length(b),length(c),length(d))

 if (color == FALSE) {col1 <- "black"; col2 ="black"; colPalette=gray(0:8/8); drape=FALSE}
 if (color == TRUE)  {col1 <- "black"; col2 ="red";   drape=TRUE}
 
 # Person parameters estimation
 parameters <- function(x, s, b, c, d, m, model="T") {
  uniform <- m4plEstimate(x , s = s, b = b, c = c, d = d, m = m,
                          model = model, prior = "uniform")
  normal  <- m4plEstimate(x , s = s, b = b, c = c, d = d, m = m,
                          model = model, prior = "normal")
  none    <- m4plEstimate(x , s = s, b = b, c = c, d = d, m = m,
                          model = model, prior = "none")
  res     <- data.frame(uniform=uniform, normal=normal, none=none)
  return(res)
  }
 parT  <- round(parameters(x, s, b, c, d, m, model="T"),2)
 parS  <- round(parameters(x, s, b, c, d, m, model="S"),2)
 parC  <- round(parameters(x, s, b, c, d, m, model="C"),2)
 parD  <- round(parameters(x, s, b, c, d, m, model="D"),2)
 param <- list(T=parT,S=parS,C=parC,D=parD)
 
 # Log Likelihood of each model
 llT     <- pggrm4pl(x=x,rep=1,n=dim(x)[2],N=dim(x)[1],theta=parT$none,
                     S=0,            C=0,            D=0,
                     s=s, b=b, c=c, d=d, log.p=TRUE)$prob
 llS     <- pggrm4pl(x=x,rep=1,n=dim(x)[2],N=dim(x)[1],theta=parS$none[1],
                     S=parS$none[2], C=0,            D=0,
                     s=s, b=b, c=c, d=d, log.p=TRUE)$prob
 llC     <- pggrm4pl(x=x,rep=1,n=dim(x)[2],N=dim(x)[1],theta=parC$none[1],
                     S=0,            C=parC$none[2], D=0,
                     s=s, b=b, c=c, d=d, log.p=TRUE)$prob
 llD     <- pggrm4pl(x=x,rep=1,n=dim(x)[2],N=dim(x)[1],theta=parD$none[1],
                     S=0,            C=0,            D=parD$none[2],
                     s=s, b=b, c=c, d=d, log.p=TRUE)$prob
 llParam <- data.frame(T=llT, S=llS, C=llC, D=llD)
 
 # Duplication of the response vector and initialisation of theta parameter
 nValues <- grain
 X       <- data.frame(matrix(as.numeric(rep(x,nValues)),ncol=nItems, byrow=TRUE))
 pTHETA  <- seq(limitT[1],limitT[2],length=nValues)
 
 # Likelihood function for the theta only model
 pS      <- seq(0,0,length=nValues)
 pC      <- seq(0,0,length=nValues)
 pD      <- seq(0,0,length=nValues)
 prob <- pggrm4pl(x=X,rep=1,n=dim(x)[2],N=dim(x)[1],theta=pTHETA,
                  S=pS, C=pC, D=pD, s=s, b=b, c=c, d=d, log.p=logLikelihood)
 plotT <-
  xyplot((prob$prob) ~ prob$theta, type="l", xlab=expression(theta), main=main, ylab="P(X)", col=col1,
   panel = function(x,y,  ...) {
    panel.xyplot(x,y, ...)
    if (annotate == TRUE) {
     vU <- param$T$uniform; vN <- param$T$normal; v <- param$T$none
     panel.abline(v=v,  lty=1, col=col1)
     posV <- ifelse (v > 0,4,2) # Right position if v < 0, else Left
     panel.text(x=v,y=max(prob$prob),paste("(",v,")", sep=""),
                pos=posV, cex=0.7)
     panel.abline(v=vN, lty=2, col=col2)
     posVN <- ifelse (v > 0,2,4) # Left position if v < 0, else Right
     panel.text(x=vN,y=max(prob$prob),paste("(",vN,")", sep=""),
                pos=posVN, cex=0.7, col=col2)
     panel.abline(v=vU, lty=3, col=col1)
     }
    }
  )
  
 # 3D likelihood function
 plot3D <-
  function(x, model="S", nValues, type="levelplot", main=NULL, color=TRUE,
           pS=seq( 0,0, length=nValues),
           pC=seq( 0,0, length=nValues),
           pD=seq( 0,0, length=nValues),
           limitT, limitS, limitC, limitD,
           s, b, c, d, logLikelihood=FALSE ) {
 pTHETA <- seq( limitT[1],limitT[2], length=nValues)
 X <- data.frame(matrix(as.numeric(rep(x,nValues^2)),ncol=length(b), byrow=TRUE))
 ## S  MODEL ...
 if (model == "S") {
  pS      <- pS + seq( limitS[1],limitS[2], length=nValues)
  grid    <- expand.grid(pTHETA=pTHETA, pS=pS)
  grid$pC <- grid$pD <- 0
  prob    <- pggrm4pl(x=X,rep=1,n=dim(x)[2],N=dim(x)[1],theta=grid$pTHETA,
                      S=grid$pS, C=0, D=0, s=s, b=b, c=c, d=d, log.p=logLikelihood)
  res <- switch(type,
  levelplot = levelplot(prob ~ theta * S, data=prob, drape=color, labels=list(cex=0.7),
                        xlab=expression(theta), ylab="Fluctuation", zlab="P(X)", main=main,
                        panel = function(x,y, subscripts, ...) {
                         panel.levelplot(x,y,subscripts, ...)
                         if (annotate == TRUE) {
                          v <- param$S$none[1]; h <- param$S$none[2]
                          panel.abline(v=v, h=h, lty=1, col=col1)
                          vN <- param$S$normal[1]; hN <- param$S$normal[2]
                          panel.abline(v=vN, h=hN, lty=3, col=col1)
                          panel.text(x=v,y=h,labels=paste("(",v,", ",h,")",sep=""), pos=3, cex=0.7)}
                          }
                         ),
  contourplot = contourplot(prob ~ theta * S, data=prob, labels=list(cex=0.7), main=main,
                            xlab=expression(theta), ylab="Fluctuation", zlab="P(X)",
                            region = FALSE,
                            panel = function(x,y, subscripts, ...) {
                             panel.levelplot(x,y,subscripts, ...)
                             if (annotate == TRUE) {
                              v <- param$S$none[1]; h <- param$S$none[2]
                              panel.abline(v=v, h=h, lty=1, col=col1)
                              vN <- param$S$normal[1]; hN <- param$S$normal[2]
                              panel.abline(v=vN, h=hN, lty=3, col=col2)
                              panel.text(x=v,y=0,labels=paste(v),                 pos=3, col=col1, cex=0.7)
                              panel.text(x=(limitT[1]+0.1),y=h,labels=paste(h),   pos=3, col=col1, cex=0.7)
                              panel.text(x=vN,y=(limitS[1]+.3),labels=paste(vN),  pos=2, col=col2, cex=0.7)
                              panel.text(x=(limitT[1]+0.3),y=hN,labels=paste(hN), pos=1, col=col2, cex=0.7)
                              }
                              }
                            ),
  wireframe  = wireframe(prob ~ theta * S, data=prob, main=main,
                         xlab=expression(theta), ylab="S", zlab="P(X)",
                         region = color, drape=color, colorkey=FALSE,
                         scales=list(arrows=FALSE, cex=0.7))
                         )
  }
 ## C MODEL ...
 if (model == "C") {
  pC      <- pC + seq( limitC[1],limitC[2], length=nValues)
  grid    <- expand.grid(pTHETA=pTHETA, pC=pC)
  grid$pS <- grid$pD <- 0
  prob    <- pggrm4pl(x=X,rep=1,n=dim(x)[2],N=dim(x)[1],theta=grid$pTHETA,
                      S=0, C=grid$pC, D=0, s=s, b=b, c=c, d=d, log.p=logLikelihood)
  res <- switch(type,
  levelplot = levelplot(prob ~ theta * C, data=prob, drape=color, labels=list(cex=0.7),
                        xlab=expression(theta), ylab="Pseudo-Guessing", zlab="P(X)", main=main,
                        panel = function(x,y, subscripts, ...) {
                         panel.levelplot(x,y,subscripts, ...)
                         if (annotate == TRUE) {
                          v <- param$C$none[1]; h <- param$C$none[2]
                          panel.abline(v=v, h=h, lty=1, col=col1)
                          vN <- param$C$normal[1]; hN <- param$C$normal[2]
                          panel.abline(v=vN, h=hN, lty=3, col=col1)
                          panel.text(x=v,y=h,labels=paste("(",v,", ",h,")",sep=""), pos=3, cex=0.7)}
                          }
                         ),
  contourplot = contourplot(prob ~ theta * C, data=prob, labels=list(cex=0.7), main=main,
                            xlab=expression(theta), ylab="Pseudo-Guessing", zlab="P(X)",
                            region = FALSE,
                            panel = function(x,y, subscripts, ...) {
                             panel.levelplot(x,y,subscripts, ...)
                             if (annotate == TRUE) {
                              v <- param$C$none[1]; h <- param$C$none[2]
                              panel.abline(v=v, h=h, lty=1, col=col1)
                              vN <- param$C$normal[1]; hN <- param$C$normal[2]
                              panel.abline(v=vN, h=hN, lty=3, col=col2)
                              panel.text(x=v,y=0,labels=paste(v),                 pos=3, col=col1, cex=0.7)
                              panel.text(x=(limitT[1]+0.1),y=h,labels=paste(h),   pos=3, col=col1, cex=0.7)
                              panel.text(x=vN,y=(limitS[1]+.01),labels=paste(vN), pos=2, col=col2, cex=0.7)
                              panel.text(x=(limitT[1]+0.3),y=hN,labels=paste(hN), pos=1, col=col2, cex=0.7)
                              }
                              }
                            ),
  wireframe  = wireframe(prob ~ theta * C, data=prob, main=main,
                         xlab=expression(theta), ylab="C", zlab="P(X)",
                         region = color, drape=color, colorkey=FALSE,
                         scales=list(arrows=FALSE, cex=0.7))
                         )
  }
 ## D MODEL ...
 if (model == "D") {
  pD      <- pD + seq( limitD[1],limitD[2], length=nValues)
  grid    <- expand.grid(pTHETA=pTHETA, pD=pD)
  grid$pC <- grid$pS <- 0
  prob    <- pggrm4pl(x=X,rep=1,n=dim(x)[2],N=dim(x)[1],theta=grid$pTHETA,
                      S=0, C=0, D=grid$pD, s=s, b=b, c=c, d=d, log.p=logLikelihood)
  res <- switch(type,
  levelplot = levelplot(prob ~ theta * D, data=prob, drape=color, labels=list(cex=0.7),
                        xlab=expression(theta), ylab="Inattention", zlab="P(X)",  main=main,
                        panel = function(x,y, subscripts, ...) {
                         panel.levelplot(x,y,subscripts, ...)
                         if (annotate == TRUE) {
                          v <- param$D$none[1]; h <- param$D$none[2]
                          panel.abline(v=v, h=h, lty=1, col=col1)
                          vN <- param$D$normal[1]; hN <- param$D$normal[2]
                          panel.abline(v=vN, h=hN, lty=3, col=col1)
                          panel.text(x=v,y=h,labels=paste("(",v,", ",h,")",sep=""), pos=3, cex=0.7)}
                          }
                         ),
  contourplot = contourplot(prob ~ theta * D, data=prob, labels=list(cex=0.7), main=main,
                            xlab=expression(theta), ylab="Inattention", zlab="P(X)",
                            region = FALSE,
                            panel = function(x,y, subscripts, ...) {
                             panel.levelplot(x,y,subscripts, ...)
                             if (annotate == TRUE) {
                              v <- param$D$none[1]; h <- param$D$none[2]
                              panel.abline(v=v, h=h, lty=1, col=col1)
                              vN <- param$D$normal[1]; hN <- param$D$normal[2]
                              panel.abline(v=vN, h=hN, lty=3, col=col2)
                              panel.text(x=v,y=0,labels=paste(v),                 pos=3, col=col1, cex=0.7)
                              panel.text(x=(limitT[1]+0.1),y=h,labels=paste(h),   pos=3, col=col1, cex=0.7)
                              panel.text(x=vN,y=(limitS[1]+.01),labels=paste(vN), pos=2, col=col2, cex=0.7)
                              panel.text(x=(limitT[1]+0.3),y=hN,labels=paste(hN), pos=1, col=col2, cex=0.7)
                              }
                              }
                            ),
  wireframe  = wireframe(prob ~ theta * D, data=prob, main=main,
                         xlab=expression(theta), ylab="D", zlab="P(X)",
                         region = color, drape=color, colorkey=FALSE,
                         scales=list(arrows=FALSE, cex=0.7))
                         )
  }

 return(res)
 }
 ## APPLYING MODELS
 plotC  <- plot3D(x=x, model="C", nValues=nValues, type=type, main=main, color=color,
                  pS=seq( 0,0, length=nValues),
                  pC=seq( 0,0, length=nValues),
                  pD=seq( 0,0, length=nValues),
                  limitT=limitT, limitS=limitS, limitC=limitC, limitD=limitD,
                  s=s, b=b, c=c, d=d, logLikelihood=logLikelihood ) #; plotC
                  
 plotS  <- plot3D(x=x, model="S", nValues=nValues, type=type, main=main, color=color,
                  pS=seq( 0,0, length=nValues),
                  pC=seq( 0,0, length=nValues),
                  pD=seq( 0,0, length=nValues),
                  limitT=limitT, limitS=limitS, limitC=limitC, limitD=limitD,
                  s=s, b=b, c=c, d=d, logLikelihood=logLikelihood  ) #; plotS
                  
 plotD  <- plot3D(x=x, model="D", nValues=nValues, type=type, main=main, color=color,
                  pS=seq( 0,0, length=nValues),
                  pC=seq( 0,0, length=nValues),
                  pD=seq( 0,0, length=nValues),
                  limitT=limitT, limitS=limitS, limitC=limitC, limitD=limitD,
                  s=s, b=b, c=c, d=d, logLikelihood=logLikelihood  ) #; plotD
 ## ...............
 return(list(parameters=param, logLikelihood=llParam, plotT=plotT, plotS=plotS, plotC=plotC, plotD=plotD))
 }







