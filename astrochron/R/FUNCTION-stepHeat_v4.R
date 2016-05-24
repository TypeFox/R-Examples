### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### stepHeat: caculate weighted mean, with MSWD and other associated statistics/plots
###           for incremental heating Ar/Ar data 
###                              (SRM: February 23-28, 2015; March 1, 2015; 
###                                    September 12-13, 2015)
###
###########################################################################


stepHeat <- function (dat,unc=1,lambda=5.463e-10,J=NULL,Jsd=NULL,CI=2,cull=-1,del=NULL,output=F,idPts=T,size=NULL,unit=1,setAr=95,color="black",genplot=T,verbose=T)
  {

   if(verbose) cat("\n----- GENERATE Ar/Ar AGE SPECTRUM AND CALCULATE WEIGHTED MEAN -----\n")

   if(is.logical(cull))
    {
      if(cull) 
       {
         cat("\n**** WARNING: option cull has been modified, and now requires a value of -1,0 or 1.\n")
         cat("              A value of -1 will be used.\n")
         cull = -1
       }
     if(!cull) 
       {
         cat("\n**** WARNING: option cull has been modified, and now requires a value of -1,0 or 1.\n")
         cat("              A value of 0 will be used.\n")
         cull = 0
       }
    }

# force cull=0 if del=!NULL
  if(!is.null(del)) cull=0
# force verbose=T if cull=-1 or 1
  if(cull != 0) verbose=T
# force genplot=T if cull=-1 or 1
  if(cull !=0 ) genplot=T

# check on number of columns
  if(dim(data.frame(dat))[2]!=7) stop("**** ERROR: input data frame must have 7 columns. Did you intend to use the function wtMean?")  

  if(!is.null(J) && is.null(Jsd)) stop("**** ERROR: Jsd must be specified.")  
  if(is.null(J) && !is.null(Jsd)) stop("**** ERROR: J must be specified.")  

  if(unit==1) unitLabel=c("Ma")
  if(unit==2) unitLabel=c("Ka")
  perAr=dat[,1]
  x=dat[,2]
  if(unc==2) dat[,3]=dat[,3]/2
  sd=dat[,3]
  KCa=dat[,4]
  Ar40=dat[,5]
  FAr=dat[,6]
  if(unc==2) dat[,7]=dat[,7]/2
  FAr_sd=dat[,7]

  if(!is.null(Jsd) && unc == 2) Jsd=Jsd/2

# Chauvenet's criterion, using raw mean. see Taylor (1982, pg. 142). 
# this should only be applied once, to the total data set.
  absX <- abs(x-mean(x))/sd(x)
  chauv <- pnorm(absX,lower.tail=TRUE) - pnorm(absX,lower.tail=FALSE)
  chauv <- (1-chauv) * length(x)

# save info for plotting  
  tag <- double(length(x))

  if (min(chauv) < 0.5)
   {
     if(verbose) cat("\n The following dates should considered for removal according to Chauvenet's criterion: \n")
     for (i in 1:length(x))
       {
        if(chauv[i] < 0.5) 
         {
           if(verbose) cat(" Datum=",i,"; Value=",x[i],"; Std. devs.=",absX[i],"; Chauvenet's value=",chauv[i],"\n")
#           tag[i] <- 1
          }
       }
   }
  
   if (min(perAr) <= 1)
    {
      if(verbose) cat("\n The following dates should be considered for removal, as they contribute <= 1% of the Ar39: \n")
      for (i in 1:length(x))
        {
          if(perAr[i] <= 1) 
           {
              if(verbose) cat(" Datum=",i,"; Value=",x[i],"; Sigma=",sd[i],"; Percent Ar39=",perAr[i],"\n")
           }
        }
    }
  
   if(!is.null(del))
    {
       if(verbose) cat(" Deleting dates:",del,"\n")
       perAr[del] <- NA
       x[del] <- NA
       sd[del] <- NA
       KCa[del] <- NA
       Ar40[del] <- NA
       FAr[del] <- NA
       FAr_sd[del] <- NA

       perAr <- subset(perAr, !(perAr == "NA"))
       x <- subset(x, !(x == "NA"))
       sd <- subset(sd, !(sd == "NA"))
       KCa <- subset(KCa, !(KCa == "NA"))
       Ar40 <- subset(Ar40, !(Ar40 == "NA"))
       FAr <- subset(FAr, !(FAr == "NA"))
       FAr_sd <- subset(FAr_sd, !(FAr_sd == "NA"))
       tag <- double(length(x))
       if(genplot) cat("\n**** WARNING: plotting disabled when del is specified!\n")
       genplot <- F
    }
   
   dat2=cbind(perAr,x,sd,KCa,Ar40,FAr,FAr_sd)

   if(genplot) 
    {
       dev.new(height=7,width=5)
       mat <- matrix(c(1,2,3), nrow=3,ncol=1)
       layout(mat, heights=c(1.5,1.5,3.5))
       if(is.null(size)) size=1.4
       sizeAxis=1.3
    }

wtMeanFun <- function(aa,bb,J95,verbose)
  {
# error checking
    if(length(aa) < 2) stop("*** ERROR: only one age available!")
    x=aa
    sd=bb
    
    npts=length(x)
    w <- sd^2
    w <- 1/w
# note that weights do not sum to one, they are 'raw' weights.
# this follows Taylor (1982, pg. 142) and McDougall and Harrison (1999, pg. 134)
    sum.w <- sum(w)
    mean.w <- sum(w*x)/sum.w
    sigma.w <-(sum.w)^(-0.5)
    sigma2.w <- 2*sigma.w
  
# calculate reduced chi-sqaure, or MSWD. This is chi-sq/(degrees of freedom)    
    x0 <- x - mean.w
# this follows McDougall and Harrison (1999, pg. 135)
    mswd3 <- (sum(x0^2/sd^2)) / (npts-1)
     
# use IsoPlot correction formula for 95% CL.
    if(CI==1) mswd95 <- qt(p=0.025,df=npts-1,lower.tail=F)*sigma.w*sqrt(mswd3)
# use ArArCALC correction formula    
    if(CI==2) mswd95 <-1.96*sigma.w*sqrt(mswd3)
    
    sigma95 <- sqrt( (1.96*sigma.w)^2 + J95^2)
    mswd95 <- sqrt( mswd95^2 + J95^2)

# calculate probabilty-of-fit, following Wendt and Carl (1991, EQ. 17. pgs. 278-279)
#   1 standard deviation of the mean MSWD is SQRT(2/f), where f is npts-1
#   the MSWD follows a chisquare distribution: it is defined as chisq/df, where df = npts-1
    pfit=pchisq(mswd3*(npts-1),df=(npts-1),lower.tail=F)

    if(verbose)
     {
       
       if(verbose) cat("\n --- WEIGHTED MEAN ESTIMATES ---")
       cat("\n Weighted mean =",mean.w,"\n")
       cat(" sigma =", sigma.w,"\n")
       cat(" 2*sigma =", sigma2.w,"\n")
       cat(" (J uncertainty NOT included) \n\n")
       
       cat(" --- 95% CONFIDENCE INTERVAL ---\n")
       if(CI==1)
        {
          if(pfit < 0.15) cat(" 95% CI, 1.96*sigma =", sigma95," \n")
          if(pfit < 0.15) cat(" 95% CI, t*sigma*sqrt(MSWD) =",mswd95,"<---- RECOMMENDED VALUE FOLLOWING ISOPLOT CONVENTION\n")   
          if(pfit >= 0.15) cat(" 95% CI, 1.96*sigma =", sigma95,"<---- RECOMMENDED VALUE FOLLOWING ISOPLOT CONVENTION \n")
          if(pfit >= 0.15) cat(" 95% CI, t*sigma*sqrt(MSWD)  =",mswd95," \n")
        }
       if(CI==2)
        { 
          if(mswd3 <= 1) cat(" 95% CI, 1.96*sigma =", sigma95,"<---- RECOMMENDED VALUE FOLLOWING ArArCALC CONVENTION \n")
          if(mswd3 <= 1) cat(" 95% CI, 1.96*sigma*sqrt(MSWD) =", mswd95," \n")
          if(mswd3 > 1) cat(" 95% CI, 1.96*sigma =", sigma95," \n")
          if(mswd3 > 1) cat(" 95% CI, 1.96*sigma*sqrt(MSWD) =", mswd95,"<---- RECOMMENDED VALUE FOLLOWING ArArCALC CONVENTION \n")
        }

       if(J95>0) cat(" (J uncertainty included) \n")
       if(J95==0) cat(" (J uncertainty NOT included) \n")

       cat("\n MSWD (of wt. mean)=",mswd3,"\n")
       cat(" Probability of fit (0-1)=",pfit,"\n")
     }

     out1=data.frame(cbind(mean.w,sigma.w,sigma2.w,sigma95,mswd95,mswd3,pfit))
     colnames(out1) <- c("wt_mean","sigma","2sigma","1.96sigma_95%","mswd_95%","mswd","prob_fit")
     return(out1)
# end wtMeanFun
  }


plotit <- function(perAr,x,sd,KCa,Ar40,mean.w,sigma.w,sigma2.w,sigma95,mswd95,mswd3,pfit,CI,cull,tag,size,sizeAxis,unitLabel,setAr,color)
  {
    npts=length(x)
# percent 39ArK released
    cumulativeAr=append(0,cumsum(dat[,1]))    
    medAr=cumulativeAr[1:npts]+dat[,1]/2

# plot Ar40
   par(mar=c(0.5, 4.5, 1, 1.5))
# set up parameters for line plots
    x0<-double(npts)
    y0<-double(npts)
    x1<-double(npts)
    y1<-double(npts)
    for (i in 1:npts)
     {
       x0[i]<-cumulativeAr[i] 
       y0[i] <- Ar40[i]
       x1[i]<-cumulativeAr[i+1]
       y1[i] <- Ar40[i]
     }
#    yrange=c(min(Ar40),max(Ar40))
    yrange=c(min(Ar40),100)
    plot(0,0,ylab="",xlab="",xlim=c(0,100),ylim=yrange,cex.axis=sizeAxis,xaxt="n",xaxs="i",bty="l")
    mtext(expression(paste("%"^40,"Ar*")),side=2,line=2.5)
    segments(x0,y0,x1,y1,col="red",lwd=2)
    abline(h = setAr, lty = 4, col = "gray", lwd = 2)
    if(idPts) 
       {
#         points(medAr,Ar40,col=color,cex=1.5*size,pch=19,xpd=TRUE)
#         text(medAr,Ar40,1:npts,col="white",cex=0.6*size,xpd=TRUE) 
         text(medAr,Ar40,1:npts,col="black",cex=0.7*size,pos=3,offset=0.2,xpd=TRUE) 
       }

# plot KCa
# set up parameters for line plots
    x0<-double(npts)
    y0<-double(npts)
    x1<-double(npts)
    y1<-double(npts)
    for (i in 1:npts)
     {
       x0[i]<-cumulativeAr[i] 
       y0[i] <- KCa[i]
       x1[i]<-cumulativeAr[i+1]
       y1[i] <- KCa[i]
     }
    yrange=c(min(KCa),max(KCa))
    plot(0,0,ylab="",xlab="",xlim=c(0,100),ylim=yrange,cex.axis=sizeAxis,xaxt="n",xaxs="i",bty="l")
    mtext("K/Ca",side=2,line=2.5)
    segments(x0,y0,x1,y1,col="red",lwd=2)
    if(idPts) 
       {
#         points(medAr,KCa,col=color,cex=1.5*size,pch=19,xpd=TRUE)
#         text(medAr,KCa,1:npts,col="white",cex=0.6*size,xpd=TRUE) 
         text(medAr,KCa,1:npts,col="black",cex=0.7*size,pos=3,offset=0.2,xpd=TRUE) 
       }


# plot ages  
   par(mar=c(4, 4.5, 3, 1.5))  
# find minimum (for plotting polygon)
    tmin=min(x-(2*sd))
    tmax=max(x+(2*sd))   
    yrange=c(tmin,tmax)
    plot(medAr,x,type="p",xlab="",ylab="",xlim=c(0,100),ylim=yrange,cex.axis=sizeAxis,col="white",xaxs="i",bty="l")
    mtext(paste("Age (",unitLabel,")",sep=""),side=2,line=2.5)
    mtext(expression(paste("%"^39,"Ar"["K"])),side=1,line=2.5)
# set up parameters for boxes
    xleft<-double(npts)
    xright<-double(npts)
    ybottom<-double(npts)
    ytop<-double(npts)
    for (i in 1:npts)
     {
       xleft[i]<-cumulativeAr[i] 
       xright[i]<-cumulativeAr[i+1]
       ybottom[i] <- x[i]+2*sd[i]
       ytop[i] <- x[i]-2*sd[i]
     }
    rect(xleft,ybottom,xright,ytop,col="gray",lwd=0.5,border="black")

# identify steps that have been selected
    if(sum(tag)>0)
        {
          ii=which(tag==1)
          xleft <- xleft[ii]
          xright <- xright[ii]
          ybottom <- ybottom[ii]
          ytop <- ytop[ii]
          rect(xleft,ybottom,xright,ytop,col="white",lwd=0.5)
        }

# plot weighted mean
    if(pfit >= 0.15 && CI == 1) rect(0,mean.w-(sigma95),100,mean.w+(sigma95), col="#FF000050",border="NA")
    if(pfit < 0.15 && CI == 1) rect(0,mean.w-mswd95,100,mean.w+mswd95, col="#FF000050",border="NA")
    if(mswd3 <= 1 && CI == 2) rect(0,mean.w-(sigma95),100,mean.w+(sigma95), col="#FF000050",border="NA")
    if(mswd3 > 1 && CI == 2) rect(0,mean.w-mswd95,100,mean.w+mswd95, col="#FF000050",border="NA")

    points(medAr,x,col=color,cex=1.5*size,pch=19,xpd=TRUE)
    if(sum(tag)>0)
        {
          points(medAr[ii],x[ii],cex=1.5*size,col="white",pch=19,xpd=TRUE) 
          points(medAr[ii],x[ii],cex=1.5*size,col="black",pch=1,xpd=TRUE) 
        }

    if(idPts) 
       {
         text(medAr,x,1:npts,col="white",cex=0.6*size,xpd=TRUE) 
         if(sum(tag)>0) text(medAr[ii],x[ii],ii,col="black",cex=0.6*size,xpd=TRUE) 
       }

    textsize=1
    if(pfit >= 0.15 && CI == 1) mtext(c(paste("Wt. Mean=",round(mean.w,digits=4),"+/-",round(sigma95,digits=4),"; MSWD=",round(mswd3,digits=2))),side=3,line=1.5,col="red",cex=textsize)
    if(pfit < 0.15 && CI == 1) mtext(c(paste("Wt. Mean=",round(mean.w,digits=4),"+/-",round(mswd95,digits=4),"; MSWD=",round(mswd3,digits=2))),side=3,line=1.5,col="red",cex=textsize)
    if(mswd3 <= 1 && CI == 2) mtext(c(paste("Wt. Mean=",round(mean.w,digits=4),"+/-",round(sigma95,digits=4),"; MSWD=",round(mswd3,digits=2))),side=3,line=1.5,col="red",cex=textsize)
    if(mswd3 > 1 && CI == 2) mtext(c(paste("Wt. Mean=",round(mean.w,digits=4),"+/-",round(mswd95,digits=4),"; MSWD=",round(mswd3,digits=2))),side=3,line=1.5,col="red",cex=textsize)

    if(pfit < 0.15) mtext(c("Warning: probability of fit is < 0.15"),col="blue",side=3,line=-1.2,cex=textsize*0.85)

    if(cull == 0) mtext(expression(paste("Boxes are 2",sigma," ; Red is wt. mean with 95% CI")),side=3,line=0,col="red")
    if(cull == -1) mtext("Click on dates for exclusion",side=3,line=0)
    if(cull == 1) mtext("Click on dates to retain",side=3,line=0)
    
# end plotit
  }    


## this script modified from '?identify' in R
identifyPch <- function(x, y=NULL, n=length(x), pch=19, cex, ...)
    {
     xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
     sel <- rep(FALSE, length(x)); res <- integer(0)
     while(sum(sel) < n) {
         ans <- identify(x[!sel], y[!sel], n=1, plot=F, ...)
         if(!length(ans)) break
         ans <- which(!sel)[ans]
         points(x[ans], y[ans], pch = pch, cex = cex, col="white",xpd=TRUE)
         points(x[ans], y[ans], pch = 1, cex = cex, col="black",xpd=TRUE)
         sel[ans] <- TRUE
         res <- c(res, ans)
       }
    res
    }

  if(!is.null(J))
   {
     wtMeanFAr <- wtMeanFun(FAr,FAr_sd,J95=0,verbose=F)  
# solve for J uncertainty
     Fval <- wtMeanFAr[,1]
     Cval <- 1-(J*Fval)
     J95 <- (Fval/(Cval*lambda))^2
     J95 <- J95 * Jsd^2
     J95 <- 1.96 * sqrt(J95)
     if(unit==1) J95 <- J95/1000000
     if(unit==2) J95 <- J95/1000
   }
  if(is.null(J)) J95=0
  
  wtMeanOut <- wtMeanFun(x,sd,J95=J95,verbose=verbose)

  if(genplot) plotit(perAr=perAr,x=x,sd=sd,KCa=KCa,Ar40=Ar40,wtMeanOut[,1],wtMeanOut[,2],wtMeanOut[,3],wtMeanOut[,4],wtMeanOut[,5],wtMeanOut[,6],wtMeanOut[,7],CI=CI,cull=cull,tag=tag,size=size,sizeAxis=sizeAxis,unitLabel=unitLabel,setAr=setAr,color=color)

  if(cull != 0) 
    {
      cat("\n * Select dates by clicking on points in the plot on the BOTTOM.\n")
      cat("   Stop by pressing ESC-key (Mac) or STOP button (Windows)\n \n")
# percent 39ArK released
      cumulativeAr=append(0,cumsum(dat[,1]))    
      medAr=cumulativeAr[1:length(dat[,1])]+dat[,1]/2
      perArCull <- perAr
      xCull <- x
      sdCull <- sd
      KCaCull <- KCa
      Ar40Cull <- Ar40
      FArCull <- FAr
      FAr_sdCull <- FAr_sd
      pts <- identifyPch(medAr,x, cex=1.5*size)
      tag[pts] <- 1

      if(cull == -1) ipts=pts   
      if(cull == 1) ipts=-pts
           
      perArCull[ipts] <- NA
      xCull[ipts] <- NA
      sdCull[ipts] <- NA
      KCaCull[ipts] <- NA
      Ar40Cull[ipts]<- NA
      FArCull[ipts] <- NA
      FAr_sdCull[ipts] <- NA
      perArCull <- subset(perArCull, !(perArCull == "NA"))
      xCull <- subset(xCull, !(xCull == "NA"))
      sdCull <- subset(sdCull, !(sdCull == "NA"))
      KCaCull <- subset(KCaCull, !(KCaCull == "NA"))
      Ar40Cull <- subset(Ar40Cull, !(Ar40Cull == "NA"))
      FArCull <- subset(FArCull, !(FArCull == "NA"))
      FAr_sdCull <- subset(FAr_sdCull, !(FAr_sdCull == "NA"))
      dat2=cbind(perArCull,xCull,sdCull,KCaCull,Ar40Cull,FArCull,FAr_sdCull)

     if(!is.null(J))
       {
         wtMeanFAr2 <- wtMeanFun(FArCull,FAr_sdCull,J95=0,verbose=F)
# solve for J uncertainty
         Fval <- wtMeanFAr2[,1]
         Cval <- 1-(J*Fval)
         J95 <- (Fval/(Cval*lambda))^2
         J95 <- J95 * Jsd^2
         J95 <- 1.96 * sqrt(J95)
         if(unit==1) J95 <- J95/1000000
         if(unit==2) J95 <- J95/1000
        }
      if(is.null(J)) J95=0   

      wtMeanOut2 <- wtMeanFun(xCull,sdCull,J95=J95,verbose=verbose)
      if(genplot) 
        {
           dev.new(height=7,width=5)
           mat <- matrix(c(1,2,3), nrow=3,ncol=1)
           layout(mat, heights=c(1.5,1.5,3.5))  
           plotit(perAr=perAr,x=x,sd=sd,KCa=KCa,Ar40=Ar40,wtMeanOut2[,1],wtMeanOut2[,2],wtMeanOut2[,3],wtMeanOut2[,4],wtMeanOut2[,5],wtMeanOut2[,6],wtMeanOut2[,7],CI=CI,cull=0,tag=tag,size=size,sizeAxis=sizeAxis,unitLabel=unitLabel,setAr=setAr,color=color)
        }
    }

   if(output && cull == 0) return(wtMeanOut)
   if(output && cull != 0) return(wtMeanOut2)
# end function stepHeat
}
