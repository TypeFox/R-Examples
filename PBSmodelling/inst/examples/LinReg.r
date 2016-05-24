# R Code to run Linear Regression demo
# ***********************************************
# Function for viewing linear regression examples
# ***********************************************

local(envir=.PBSmodEnv,expr={
locale = sys.frame(sys.nframe() - 1) # local environment

pairsLR <- function() {
   # ----------------------------
   # Pairs plot of chosen dataset
   # ----------------------------
   getWinVal(scope="L")
   if (dset=="sim") {
      x   <- runif(nsim,xmin,xmax)
      eps <- rnorm(nsim,0,1)
      y   <- asim + bsim * x + (eps * ssim)
      sim <- data.frame(X=x,Y=y); tput(sim)
   }
   file <- get(dset,pos=.PBSmodEnv); writeList(file,"LinRegDat.txt");
   flds <- names(file); nfld <- ncol(file)
   msg  <- paste("# Records =",nrow(file))
   msg  <- paste(msg,paste(1:nfld,"=",flds,collapse="\n"),sep="\n")
   setWinVal(list(allflds=msg))
   pairs(file,gap=0)

}

plotLR <- function(addBRugs=NUL, ...) {
   # ----------------------------------------------------------------
   # Compute classical fit with predictions and a confidence interval
   # ----------------------------------------------------------------
   getWinVal(scope="L")
   resetGraph(); par(ask=FALSE);
   if (is.null(addBRugs)) addBRugs <- as.logical(getWinAct()[1])

   file <- get(dset,pos=.PBSmodEnv)
   flds <- names(file); nfld <- ncol(file)
   msg  <- paste("# Records =",nrow(file))
   msg  <- paste(msg,paste(1:nfld,"=",flds,collapse="\n"),sep="\n")
   if (xfld>nfld) {
      msg <- paste(msg,paste("Choose X from 1:",nfld,sep=""),sep="\n")
      setWinVal(list(allflds=msg))
      return()
   }
   if (yfld>nfld) {
      msg <- paste(msg,paste("Choose Y from 1:",nfld,sep=""),sep="\n")
      setWinVal(list(allflds=msg))
      return()
   }
   setWinVal(list(allflds=msg))
   n    <- nrow(file)
   x    <- file[,xfld]; xlim <- range(x,na.rm=T); xdiff <- diff(xlim)
   y    <- file[,yfld]; ylim <- range(y,na.rm=T); ydiff <- diff(ylim)

   clrL     <- c("green","palevioletred","cornflowerblue","darkgray")
   clrH     <- c("forestgreen","red","blue","darkviolet")

   xlo   <- min(x,na.rm=T); xhi <- max(x,na.rm=T)
   lrfit <- lm(y ~ x); a <- lrfit$coeff[1]; b <- lrfit$coeff[2];
   up    <- ifelse(b>0,T,F)
   xp    <- data.frame( x=seq(xlo,xhi,length=100) );
   yyc   <- predict.lm(lrfit,xp,interval="confidence",level=plev);
   yyp   <- predict.lm(lrfit,xp,interval="prediction",level=plev);
   yc1   <- yyc[,2]; yc2 <- yyc[,3];  # confidence lower and upper bounds
   yp1   <- yyp[,2]; yp2 <- yyp[,3];  # prediction lower and upper bounds
   xp    <- t(xp);                    # convert to vector
   yrng  <- c(min(y,yp1),max(y,yp2)); # y range, allowing for intervals

   xlim <- c(xlo,xhi); ylim <- range(yyp)
   xlim <- xlim + (c(-1,1) * .05 * xdiff);
   #ylim <- ylim + (c(-1,1) * .05 * ydiff);

   plot(x,y,xlim=xlim,ylim=ylim,type="n",
     xlab=flds[xfld],ylab=flds[yfld],
     mgp=c(2.5,.75,0),cex.axis=1.2,cex.lab=1.5);
   if (addBRugs) modSplash()
   if (dset=="sim") abline(asim,bsim,col="blue",lwd=2)
   abline(lrfit$coefficients,lwd=2,col=clrL[3]);
   points(x,y,col=clrL[1],pch=16,cex=1.5);
   points(x,y,col=clrH[1],pch=1,cex=1.5);
   lines(xp,yc1,col=clrH[2],lwd=2,lty=1);
   lines(xp,yc2,col=clrH[2],lwd=2,lty=1);
   lines(xp,yp1,col=clrH[4],lwd=2,lty=5); 
   lines(xp,yp2,col=clrH[4],lwd=2,lty=5);
   addLabel(ifelse(up,.75,.05),.10,paste(
     paste(c("a =","b ="),signif(c(a,b),3)),collapse="\n"),
     cex=1.5,adj=0)
   box()
}

# **************
#  WinBUGS model
# **************

modCompile <- function() {
   # ----------------------------------------
   # Initialize and compile the WinBUGS model
   # ----------------------------------------
   getWinVal(scope="L");
   file <- get(dset,pos=.PBSmodEnv)
   flds <- names(file); nfld <- ncol(file)
   n    <- nrow(file)
   x    <- file[,xfld]
   y    <- file[,yfld]
   dat  <- list(n=n,x=x,y=y);
   bugsData(dat,"LRDat.txt");     # write data file (BRugs function)
   modelCheck("LinRegMod.txt");   # check model syntax
   modelData("LRDat.txt");        # load current data
   modelCompile(nc);              # compile with nc chains
   modelGenInits();               # generate randoms inits
   samplesSet(c("a","b","sig")[pset]);  # parameters to monitor
   setWinVal(list(clen=1000,cthin=1,ctot=0,s1=1,s2=1000,sthin=1,chn2=nc))
}

modUpdate <- function() {
   # -------------------------------------------------------------
   # Update the model and save complete history in global "LRhist"
   # -------------------------------------------------------------
   getWinVal(scope="L");
   modelUpdate(clen,cthin);
   LRhist <- as.data.frame( samplesHistory("*",beg=0,plot=F) );
   if (nc==1) names(LRhist) <- paste(names(LRhist),1,sep=".");
   LRhist <- LRhist; tput(LRhist)
   ctot <- dim(LRhist)[1];    # total length so far
   setWinVal(list(ctot=ctot,s1=ctot-clen+1,s2=ctot)); par(ask=F);
}

# ---------------------------
# Functions to report results
# ---------------------------

panel.hist <- function(x, ...) {
   usr <- par("usr"); on.exit(par(usr))
   h <- hist(x, breaks="FD", plot=FALSE)
   breaks <- h$breaks; nB <- length(breaks)
   y <- h$counts; y <- y/sum(y)
   par(usr = c(usr[1:2], 0, max(y)*1.5) )
   rect(breaks[-nB], 0, breaks[-1], y, col="salmon")
   box()
}

modHist <- function() {
   getWinVal(scope="L"); resetGraph(); nr <- sum(pset);
   samplesHistory("*",beg=s1-1,end=s2-1,thin=sthin,mfrow=c(nr,1),
     ask=F,firstChain=chn1,lastChain=chn2); };

modDens <- function() {
   getWinVal(scope="L"); resetGraph(); nr <- sum(pset);
   samplesDensity("*",beg=s1-1,end=s2-1,thin=sthin,mfrow=c(nr,1),
     ask=F,firstChain=chn1,lastChain=chn2); };

modACF <- function() { # default chain 1
   getWinVal(scope="L");
   resetGraph(); nr <- sum(pset);
   samplesAutoC("*",chain=chn1,beg=s1-1,end=s2-1,thin=sthin,
     mfrow=c(nr,1),ask=F); };

modPairs <- function() {
   getWinVal(scope="L"); resetGraph();
   i1   <- max(s1,1); i2 <- min(s2,ctot); # ensure valid range
   idx  <- seq(i1,i2,by=sthin);
   tget(LRhist); file <- LRhist[idx,]
   if (dset=="sim") {
      pars <- rep(c(asim,bsim,ssim),each=nc)
      file <- rbind(file,pars)
   }
   nams <- names(file); dot <- regexpr("\\.",nams);
   z <- is.element(substring(nams,dot+1),chn1:chn2)
   file <- file[,z]
   par(ask=FALSE);
   pairs(file, diag.panel=panel.hist, gap=0, cex.labels=1.5,
      panel=function(x,y,s=ifelse(dset=="sim",TRUE,FALSE)) {
         n <- length(x); nn <- n-1;
         points(x[1:ifelse(s,nn,n)],y[1:ifelse(s,nn,n)],
           pch=16,cex=0.6,col="darkgray");
         if (s) { abline(h=y[n],v=x[n],col="blue",lty=3);
                  points(x[n],y[n],col="cornflowerblue",pch=16,cex=1.5);
                  points(x[n],y[n],col="blue",pch=1,cex=1.5); }
      })
}

modSub <- function(chains=samplesGetFirstChain():samplesGetLastChain()) {
   if (!exists("LRhist",where=.PBSmodEnv)) {
      return(FALSE)
   }
   getWinVal(scope="L");
   i1  <- max(s1,1); i2 <- min(s2,ctot); # ensure valid range
   idx <- seq(i1,i2,by=sthin);

   Pfld <- c("a","b","sig")[pset]; nP <- length(Pfld)
   nch  <- length(chains)
   pfld <- paste(rep(Pfld,each=nch),chains,sep=".")
   tget(LRhist); file <- LRhist[idx,pfld]

   temp <- NULL
   for (i in chains) {
      junk <- file[,paste(Pfld,i,sep=".")]
      names(junk) <- Pfld
      temp <- rbind(temp,junk)
   }
   LRsub <- temp; tput(LRsub)
   return(TRUE)
}

modFreq <- function() {
   getWinVal(scope="L");
   chains <- chn1:chn2
   resetGraph(); par(ask=FALSE);

   Pfld <- c("a","b","sig")[pset]; nP <- length(Pfld)
   modSub(chains)
   tget(LRsub); file <- LRsub

   iclr <- c("aquamarine","gold","plum")
   par(mfrow=c(nP,1),mai=c(.25,0.6,.05,.05),omi=c(0,0,0,0),ask=F)

   for (i in Pfld) {
      ii   <- match(i,c("a","b","sig"))
      x    <- file[,i]
      q95  <- quantile(x,c(.025,.975))
      x95  <- c(rep(q95[1],2),NA,rep(q95[2],2))
      xmn  <- mean(x)

      xf   <- hist(x,nclass=25,plot=F)
      xx   <- xf$mids
      yy   <- xf$counts/sum(xf$counts)
      xoff <- diff(xf$breaks)[1]/2
      nn   <- length(xx)
      xxx  <- as.vector(rbind(xx-xoff,xx+xoff,xx+xoff,xx-xoff,rep(NA,nn)))
      yyy  <- as.vector(rbind(rep(0,nn),rep(0,nn),yy,yy,rep(NA,nn)))
      xlim <- range(xxx,na.rm=T); ylim <- c(0,max(yy)); ydiff <- diff(ylim)
      ylim <- ylim + c(0,.05*ydiff)

      plot(0,0,xlim=xlim,ylim=ylim,type="n",mgp=c(3.5,.5,0),xaxt="n",yaxt="n",
        axes=F,xlab="",ylab="Proportion",cex.lab=1.5)
      axis(1,at=xx-xoff,pos=0,mgp=c(0,.5,0),tck=-.01)
      axis(2,at=pretty(ylim),cex.axis=1.2,adj=1,mgp=c(0,.6,0),las=1)
      polygon(xxx,yyy,col=iclr[ii])

      y75 <- ylim[2]*.75; y75t <- y75*1.05
      y95<- ylim[2]*.95;  y95t <- y95*1.05

      ycl <- c(0,y75,NA,0,y75)
      lines(x95,ycl,col="red",lty=5)
      text(q95,rep(y75t,2),round(q95,2))
      lines(rep(xmn,2),c(0,y95),col="red",lwd=2)
      text(xmn,y95t,round(xmn,2))
      addLabel(.90,.90,i,cex=1.5,col="blue")
   }
}

modSplash <- function() {
   getWinVal(scope="L");
   chains <- chn1:chn2

   Pfld <- c("a","b","sig")[pset]; nP <- length(Pfld)
   isOK <- modSub(chains)
   if (isOK) {
      tget(LRsub); file <- LRsub
      iclr <- c("goldenrod","aquamarine","plum")
      addab <- function(x) {
         abline(a=x[1],b=x[2],col="gold")
      }
      apply(file,1,addab)
   }
}

# ********************************
# Load libraries and start the GUI
# ********************************

if (!require(BRugs, quietly=TRUE)) stop("The BRugs package is required for this example")
if (!require(PBSmodelling, quietly=TRUE)) stop("The PBSmodelling package is required for this example")
createWin("LinRegWin.txt"); pairsLR()

}) # end local scope
