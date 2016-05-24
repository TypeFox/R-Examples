### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### wtMean: caculate weighted mean, with MSWD and other associated statistics/plots. 
###                              (SRM: March 16, 2012; Sept. 22-25, 2014; 
###                               February 6-7, 2015; February 10-15, 2015;
###                               February 24-28, 2015; March 1, 2015;
###                               September 11-13, 2015)
###
###########################################################################


wtMean <- function (dat,sd=NULL,unc=1,lambda=5.463e-10,J=NULL,Jsd=NULL,CI=2,cull=-1,del=NULL,sort=1,output=F,idPts=T,size=NULL,unit=1,setAr=95,color="black",genplot=T,verbose=T)
  {

   if(verbose) cat("\n----- CALCULATE WEIGHTED MEAN AND GENERATE SUMMARY PLOTS-----\n")

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

# check on dimensions of x
  if(dim(data.frame(dat))[2] == 1) 
    {
      pltype=1 
      if(is.null(sd)) stop("**** ERROR: sd must be specified!")  
      dat=cbind(dat,sd)
    } 

# standard deviation is in second column of input
  if(dim(data.frame(dat))[2] == 2) pltype=1

# standard deviation is in second column, K/Ar is in third, %Ar40* is in fourth, F is in fifth, F std. dev. is in sixth
  if(dim(data.frame(dat))[2] == 6) pltype=2

  if(dim(data.frame(dat))[2] > 6) stop("**** ERROR: your data frame has too many columns. Did you intend to use the function stepHeat?")  

  if(!is.null(J) && is.null(Jsd)) stop("**** ERROR: Jsd must be specified.")  
  if(is.null(J) && !is.null(Jsd)) stop("**** ERROR: J must be specified.")  

  if(unit==1) unitLabel=c("Ma")
  if(unit==2) unitLabel=c("Ka")

# sort by date
  if(sort==1) dat <- dat[order(dat[,1], na.last = NA, decreasing = F), ]
  if(sort==2) dat <- dat[order(dat[,1], na.last = NA, decreasing = T), ]

  x=dat[,1]
  if(unc==2) dat[,2]=dat[,2]/2
  sd=dat[,2]
  if(pltype==2)
   {      
     KCa=dat[,3]
     Ar40=dat[,4]
     FAr=dat[,5]
     if(unc==2) dat[,6]=dat[,6]/2
     FAr_sd=dat[,6]
   }

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
           tag[i] <- 1
          }
         }
   }
  
   if(!is.null(del))
    {
       if(verbose) cat(" Deleting dates:",del,"\n")
       x[del] <- NA
       sd[del] <- NA
       x <- subset(x, !(x == "NA"))
       sd <- subset(sd, !(sd == "NA"))
       if(pltype==2)
        {
          KCa[del] <- NA
          Ar40[del] <- NA
          FAr[del] <- NA
          FAr_sd[del] <- NA
          KCa <- subset(KCa, !(KCa == "NA"))
          Ar40 <- subset(Ar40, !(Ar40 == "NA"))
          FAr <- subset(FAr, !(FAr == "NA"))
          FAr_sd <- subset(FAr_sd, !(FAr_sd == "NA"))
        }   
       tag <- double(length(x))
    }

   if(pltype==1) 
     {
       dat2=cbind(x,sd)
       if(is.null(size)) size=1
       sizeAxis=1
     }
   if(pltype==2)
     {
       dat2=cbind(x,sd,KCa,Ar40,FAr,FAr_sd)
       if(is.null(size)) size=1.2
       sizeAxis=1.5
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


plotit <- function(dat2,wtMeanIn,CI,cull,tag,pltype,idPts,size,sizeAxis,unitLabel,setAr,color)
  {
# error checking
    if(length(dat2[,1]) < 2) stop("*** ERROR: only one age available!")  
    
    x=dat2[,1]
    sd=dat2[,2]
    if(pltype==2)
     {      
       KCa=dat2[,3]
       Ar40=dat2[,4]
       FAr=dat2[,5]
       FAr_sd=dat2[,6]
     }    
    npts=length(x)
  
    mean.w=wtMeanIn[,1]
    sigma.w=wtMeanIn[,2]
    sigma2.w=wtMeanIn[,3]
    sigma95=wtMeanIn[,4]
    mswd95=wtMeanIn[,5]
    mswd3=wtMeanIn[,6]
    pfit=wtMeanIn[,7]

    if(pltype==1)
     {
       dev.new(height=4.8,width=9)
       par(mfrow=c(1,2),mar=c(4,4,3,1))
     }
      if(pltype==2)
     {
       dev.new(height=4.5,width=11.5)
       par(mfrow=c(1,3),mar=c(4,4,4,1))
     }
    
# generate sum of gaussians
    tmin=min(x-(2*sd))
    tmax=max(x+(2*sd))
    gaus=double(1000*npts)
    dim(gaus) <- c(1000,npts)
    for (i in 1:npts)
       {
         gaus[,i]=dnorm(seq(tmin,tmax,length=1000),x[i],sd[i])
       }  
# add up gaussians
    totalGaus=rowSums(gaus)
# find minimum (for plotting polygon)
    minGaus=min(totalGaus)    
    
    resQQ=qqnorm(x, main="" ,xlab="",ylab="",cex=1.5*size,col="red",pch=19,cex.axis=sizeAxis)
    mtext("Normal Q-Q plot",side=3,line=1)
    mtext(paste("Age (",unitLabel,")",sep=""),side=2,line=2.5)
    mtext("Theoretical Quantiles",side=1,line=2.5)
    qqline(x, col = "red")
    if(idPts) text(resQQ$x,resQQ$y,1:npts,col="white",cex=0.6*size)
    xrange=c(tmin,tmax)
    yrange=c(1,npts+2)
    if(pltype==2)
        {
          plot(KCa,1:npts,type="p",ylab="",xlab="",xlim=c(min(KCa),max(KCa)),ylim=yrange,col="forestgreen",lwd=2,yaxt="n",xaxt="n",cex=1.5*size,pch=19,cex.axis=sizeAxis)
          lines(KCa,1:npts,lty=3,col="forestgreen")
          axis(1, xlim=c(min(KCa),max(KCa)),lwd=1,col="forestgreen",cex.axis=sizeAxis,col.axis="forestgreen")
          mtext("K/Ca",side=1,line=2.5,col="forestgreen")
          if(idPts) text(KCa,1:npts,1:npts,col="white",cex=0.6*size)
          par(new=T)
          plot(Ar40,1:npts,type="p",ylab="",xlab="",xlim=c(0,100),ylim=yrange,col="red",lwd=2,yaxt="n",axes=F,cex=1.5*size,pch=19)
          lines(Ar40,1:npts,lty=3,col="red")
          axis(3, xlim=c(0,100),lwd=1,col="red",cex.axis=sizeAxis,col.axis="red")
          mtext(expression(paste("%"^40,"Ar*")),side=3,line=2.2,col="red")
          if(idPts) text(Ar40,1:npts,1:npts,col="white",cex=0.6*size)
          abline(v=setAr,lty=4,col="red",lwd=2)
        }
    plot(seq(tmin,tmax,length=1000),totalGaus,type="l",xlab="",ylab="",xlim=xrange,ylim=c(0,1.1*max(totalGaus)),yaxt="n",col="black",lwd=1,cex.axis=sizeAxis)
    polyGaus=append(minGaus,totalGaus)
    polyGaus=append(polyGaus,minGaus)
    polyAge=append(tmin-(tmin*10^-10),seq(tmin,tmax,length=1000))
    polyAge=append(polyAge,tmax+(tmax*10^-10))
    polygon(polyAge,polyGaus,col="#BEBEBE5A",border=NA)
    par(new = TRUE)
#    plot(x,1:npts,xlim=xrange,ylim=yrange,xlab="",ylab="",main="",yaxt="n",xaxt="n",cex=1.5*size,col="red",pch=19)
# initialize plot
    plot(0,0,xlim=xrange,ylim=yrange,xlab="",ylab="",main="",yaxt="n",xaxt="n",cex=1.5*size,col="white",pch=19)
# set up parameters for weighted mean box (with 95% CI)
    if(pfit >= 0.15 && CI == 1) 
      {
        xleft=mean.w-sigma95
        xright=mean.w+sigma95
      } 
    if(pfit < 0.15 && CI == 1) 
      {
        xleft=mean.w-mswd95
        xright=mean.w+mswd95
      }     
    if(mswd3 <= 1 && CI == 2)  
      {
        xleft=mean.w-sigma95
        xright=mean.w+sigma95      
      }
    if(mswd3 > 1 && CI == 2) 
      {
        xleft=mean.w-mswd95
        xright=mean.w+mswd95
      }
    rect(xleft,1,xright,npts+1.5,col="#0000FF1E",lwd=0.5,border=NA)
    segments(xleft,npts+1.5,xright,npts+1.5,col="blue",lwd=2)
    
# now plot points
    points(mean.w,npts+1.5,col="blue",cex=0.5*size,pch=19)
    points(x,1:npts,cex=1.5*size,col="red",pch=19)

#       atitle=bquote(Age~(.(unit))~2*sigma)
    mtext(paste("Age (",unitLabel,")",sep=""),side=1,line=2.5)
    xa<-double(npts)
    xb<-double(npts)
    for (i in 1:npts)
        {
          xa[i]=x[i]-2*sd[i]
          xb[i]=x[i]+2*sd[i]
        }
#       points(x,1:npts,cex=1.5*size,col="red",pch=19)
    segments(xa,1:npts,xb,1:npts,col="red",lwd=2)
    if(sum(tag)>0)
         {
           ii=which(tag==1)
           points(x[ii],ii,cex=1.5*size,col="black",pch=19) 
         }
    if(idPts) text(x,1:npts,1:npts,col="white",cex=0.6*size) 

#    points(mean.w,npts+1,col="blue",cex=0.5*size,pch=19)
#    xa=mean.w-sigma2.w
#    xb=mean.w+sigma2.w
#    segments(xa,npts+1,xb,npts+1,col="blue",lwd=2)
    
#    if(pltype==1) textsize=1
#    if(pltype==2) textsize=1.5
     textsize=1

#    if(pfit >= 0.15 && CI == 1) text(xrange[1]+(xrange[2]-xrange[1])/2,npts+2,c(paste("Wt. Mean=",round(mean.w,digits=2),"+/-",round(sigma95,digits=2),"; MSWD=",round(mswd3,digits=2))),col="blue",cex=textsize)
#    if(pfit < 0.15 && CI == 1) text(xrange[1]+(xrange[2]-xrange[1])/2,npts+2,c(paste("Wt. Mean=",round(mean.w,digits=2),"+/-",round(mswd95,digits=2),"; MSWD=",round(mswd3,digits=2))),col="blue",cex=textsize)
#    if(mswd3 <= 1 && CI == 2)  text(xrange[1]+(xrange[2]-xrange[1])/2,npts+2,c(paste("Wt. Mean=",round(mean.w,digits=2),"+/-",round(sigma95,digits=2),"; MSWD=",round(mswd3,digits=2))),col="blue",cex=textsize)
#    if(mswd3 > 1 && CI == 2) text(xrange[1]+(xrange[2]-xrange[1])/2,npts+2,c(paste("Wt. Mean=",round(mean.w,digits=2),"+/-",round(mswd95,digits=2),"; MSWD=",round(mswd3,digits=2))),col="blue",cex=textsize)
    if(pfit >= 0.15 && CI == 1) mtext(c(paste("Wt. Mean=",round(mean.w,digits=4),"+/-",round(sigma95,digits=4),"; MSWD=",round(mswd3,digits=2))),col="blue",side=3,line=1.5,cex=textsize)
    if(pfit < 0.15 && CI == 1) mtext(c(paste("Wt. Mean=",round(mean.w,digits=4),"+/-",round(mswd95,digits=4),"; MSWD=",round(mswd3,digits=2))),col="blue",side=3,line=1.5,cex=textsize)
    if(mswd3 <= 1 && CI == 2)  mtext(c(paste("Wt. Mean=",round(mean.w,digits=4),"+/-",round(sigma95,digits=4),"; MSWD=",round(mswd3,digits=2))),col="blue",side=3,line=1.5,cex=textsize)
    if(mswd3 > 1 && CI == 2) mtext(c(paste("Wt. Mean=",round(mean.w,digits=4),"+/-",round(mswd95,digits=4),"; MSWD=",round(mswd3,digits=2))),col="blue",side=3,line=1.5,cex=textsize)

    if(pfit < 0.15) mtext(c("Warning: probability of fit is < 0.15"),col="blue",side=3,line=-1.2,cex=textsize*0.85)

    if(cull == 0) mtext(expression(paste("Red = 2",sigma," ; Blue = wt. mean with 95% CI")),side=3,line=0,col="black")
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
         points(x[ans], y[ans], pch = pch, cex = cex, col="white")
         points(x[ans], y[ans], pch = 1, cex = cex, col="black")
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

  if(genplot) plotit(dat2=dat2,wtMeanIn=wtMeanOut,CI=CI,cull=cull,tag=tag,pltype=pltype,idPts=idPts,size=size,sizeAxis=sizeAxis,unitLabel=unitLabel,setAr=setAr,color=color)

  if(cull != 0) 
    {
     cat("\n * Select dates by clicking on points in the plot on the BOTTOM.\n")
     cat("   Stop by pressing ESC-key (Mac) or STOP button (Windows)\n \n")
# percent 39ArK released
     xCull <- x
     sdCull <- sd
     if(pltype==2)
      {      
       KCaCull <- KCa
       Ar40Cull <- Ar40
       FArCull <- FAr
       FAr_sdCull <- FAr_sd
      }

     pts <- identifyPch(x,1:length(x), cex=1.5*size)

     tag[pts] <- 1
     if(cull == -1) ipts=pts   
     if(cull == 1) ipts=-pts
    
     xCull[ipts] <- NA
     sdCull[ipts] <- NA
     if(pltype==2)
      {      
        KCaCull[ipts] <- NA
        Ar40Cull[ipts]<- NA
        FArCull[ipts] <- NA
        FAr_sdCull[ipts] <- NA
      }      
       
     xCull <- subset(xCull, !(xCull == "NA"))
     sdCull <- subset(sdCull, !(sdCull == "NA"))
     if(pltype==2)
      {            
        KCaCull <- subset(KCaCull, !(KCaCull == "NA"))
        Ar40Cull <- subset(Ar40Cull, !(Ar40Cull == "NA"))
        FArCull <- subset(FArCull, !(FArCull == "NA"))
        FAr_sdCull <- subset(FAr_sdCull, !(FAr_sdCull == "NA"))
      }  
     if(pltype==1) dat2=cbind(xCull,sdCull)
     if(pltype==2) dat2=cbind(xCull,sdCull,KCaCull,Ar40Cull,FArCull,FAr_sdCull)

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
 
     tag <- double(length(x))    
     wtMeanOut2 <- wtMeanFun(xCull,sdCull,J95=J95,verbose=verbose)
     if(genplot) plotit(dat2=dat2,wtMeanIn=wtMeanOut2,CI=CI,cull=0,tag=tag,pltype=pltype,idPts=idPts,size=size,sizeAxis=sizeAxis,unitLabel=unitLabel,setAr=setAr,color=color)
# end cull
   }

   if(output && cull == 0) return(wtMeanOut)
   if(output && cull != 0) return(wtMeanOut2)
# end function wtMean    
}
