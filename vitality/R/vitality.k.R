## Three parameter simple model r s k, no childhood hook

#' Fitting routine for the 2-process, 3-parameter vitality model. 
#' 
#' Based on code by D.H. Salinger, J.J. Anderson and O. Hamel (2003).
#' "A parameter fitting routine for the vitality based survival model."
#' Ecological Modeling 166(3): 287--294.
#' 
#' @param time Vector. Time component of data: Defaults to \code{0:(1-length(sdata))}.
#' @param sdata Required. Survival or mortality data.  The default expects cumulative 
#'        survival fraction.  If providing incremental mortality fraction 
#'        instead, use option: datatype = "INC". 
#'        The default also expects the data to represent full mortality. 
#'        Otherwise, use option: rc.data = T to indicate right censored data.
#' @param rc.data Optional, boolean. Specifies Right Censored data.   If the data does not 
#'        represent full mortality, it is probably right censored.  The default 
#'        is rc.data = F.  A third option is rc.data = "TF".  Use this case to add
#'        a near-term zero survival point to data which displays nearly full 
#'        mortality ( <.01 survival at end).  If rc.data = F but the data does
#'        not show full mortality, rc.data = "TF" will be 
#'        invoked automatically. 
#' @param se Optional, boolean. Calculates the standard errors for the MLE parameters.
#'        Default is FALSE. Set equal to the initial study population to 
#'        compute standard errors. 
#' @param gfit Provides a Pearson C type test for goodness of fit.  
#'        Default is \code{gfit=F}. Set to initial study population as with \code{se} for 
#'	      computing goodness of fit. 
#' @param datatype Optional. Defaults to \code{"CUM"} for cumulative survival fraction data.
#'        Use \code{"INC"} - for incremental mortality fraction data. 
#' @param ttol Optional. Stopping criteria tolerance.  Default is 1e-6.
#'        Specify as ttol = .0001. If one of the liklihood plots (esp. for "k") does not look optimal, 
#'        try decreasing ttol.   If the program crashes, try increasing ttol.
#' @param init.params Optional. Please specify the initial param values.
#'        specify \code{init.params = c(r, s, k)} in that order
#'        (eg. init.params = c(.1, .02, .3)).
#' @param pplot Optional, boolean. Plots of cumulative survival for both data and fitted curves?
#'        Default \code{TRUE}. \code{FALSE} Produce no plots. Note:  the incremental 
#'        mortality plot is a continuous representation of the appropriately-
#'        binned histogram of incremental mortalities.
#' @param Iplot Optional, boolean. Incremental mortality for both data and fitted curves?
#'        Default: \code{FALSE}.
#' @param lplot Provides likelihood function plotting (default=\code{FALSE}).  
#'        Note:  these plots are not "likelihood profiles" in that while one 
#'        parameter is varied, the others are held fixed, rather than 
#'        re-optimized. (must also have \code{pplot=T}.)
#' @param cplot Provides a likelihood contour plot for a range of r and s values 
#'        (can be slow so default is \code{FALSE}).   Must also have lplot=T (and pplot=T) 
#'         to get contour plots.
#' @param tlab Optional, character. specifies units for x-axis of plots.  Default is "days".
#' @param silent Optional, boolean. Stops all print and plot options (still get most warning and all 
#'          error messages) Default is \code{FALSE}.  A third option, \code{"verbose"} also 
#'          enables the trace setting in the ms (minimum sum) S-Plus routine.
#' @export 
#' @return vector of final MLE r, s, k parameter estimates.
#'     standard errors of MLE parameter estimates (if se = <population> is specified).
vitality.k <- function(time,sdata,rc.data=F,se=F,gfit=F,datatype="CUM",ttol=.000001,init.params=F,lower = c(0, -1, 0),
                       upper = c(100,50,50), pplot=T,tlab="days",lplot=F,cplot=F,Iplot=F,silent=F)
  #
  #  
  #
  #     Vitality based survival model: parameter fitting routine:    VERSION: 11/14/2014
  #
  # REQUIRED PARAMETERS: 
  #     time - time component of data: time from experiment start.  Time should 
  #	     start after the imposition of a stressor is completed.
  #     sdata - survival or mortality data.  The default expects cumulative 
  #	     survival fraction.  If providing incremental mortality fraction 
  #          instead, use option: datatype="INC".
#          The default also expects the data to represent full mortality.   
#          Otherwise, use option: rc.data=T to indicate right censored data.
#
# OPTIONAL PARAMETERS:
#     rc.data =T  - specifies Right Censored data.   If the data does not 
#          represent full mortality, it is probably right censored.  The default 
#          is rc.data=F.  A third option is rc.data="TF".  Use this case to add
#          a near-term zero survival point to data which displays nearly full 
#          mortality ( <.01 survival at end).  If rc.data=F but the data does
#	     not show full mortality, rc.data="TF" will be 
#          invoked automatically. 
#     se =<population>  calculates the standard errors for the MLE parameters.  
#          Default is se=F.  The initial study population is necessary for 
#          computing these standard errors. 
#     gfit =<population>  provides a Pearson C type test for goodness of fit.  
#	     Default is gfit=F. The initial study population is necessary for 
#	     computing goodness of fit. 
#     datatype ="CUM" -cumulative survival fraction data- is the default.   
#          Other option: datatype="INC" - for incremental mortality fraction 
#          data.  ttol (stopping criteria tolerence.)  Default is .000001 . 
#          specify as ttol=.0001.
#          If one of the liklihood plots (esp. for "k") does not look optimal, 
#	     try decreasing ttol.   If the program crashes, try increasing ttol.
#     init.params =F  has the routine choose initial parameter estimates for 
#          r,s,k  (default: =F). If you wish to specify initial param values 
#          rather than have the routine choose them, specify 
#          init.params=c(r,s,k) in that order (eg. init.params=c(.1,.02,.003)).
#     pplot =T provides plots of cumulative survival and incremental mortality - 
#	     for both data and fitted curves (default: =T).  pplot=F provides no 
#	     plotting.  A third option:  pplot=n  (n>=1) extends the time axis of 
#	     the fitting plots 	(beyond the max time in data).  For example: 
#          pplot=1.2 extends the time axis by 20%.  (Note:  the incremental 
#          mortality plot is a continuous representation of the appropriately-
#          binned histogram of incremental mortalities.)
#     tlab ="<time units>" specifies units for x-axis of plots.  Default is 
#          tlab="days".
#     lplot =T provides likelihood function plotting (default =T).  
#          Note:  these plots are not "likelihood profiles" in that while one 
#          parameter is varied, the others are held fixed, rather than 
#          re-optimized. (must also have pplot=T.)
#     cplot =T provides a likelihood contour plot for a range of r and s values 
#          (can be slow so default is F).   Must also have lplot=T (and pplot=T) 
#          to get contour plots.
#     silent =T stops all print and plot options (still get most warning and all 
#          error messages) Default is F.  A third option, silent="verbose" also 
#          enables the trace setting in the ms (minimum sum) S-Plus routine. 
#
# RETURN:
#     vector of final MLE r,s,k parameter estimates.
#     standard errors of MLE parameter estimates (if se=<population> is
#     specified).
#
{
    #  --Check/prepare Data---
    datatype <- match.arg(datatype)
    if (length(time) != length(sdata)) {
        stop("time and sdata must have the same length")
    }
    in.time <- time
    dTmp <- dataPrep(time, sdata, datatype, rc.data)
    time <- dTmp$time
    sfract <- dTmp$sfract
    x1 <- dTmp$x1
    x2 <- dTmp$x2
    Ni <- dTmp$Ni
    rc.data <- dTmp$rc.data
    if(in.time[1]>0){
      time <- time[-1]
      sfract <- sfract[-1]
      x1 <- c(x1[-c(1,length(x1))], x1[1])
      x2 <- c(x2[-c(1,length(x2))], 0)
      Ni <- Ni[-1]
      rc.data <- rc.data[-1]
    }
    
    
    #  --Produce initial parameter values---
    tt <- time
    sf <- sfract
    if(length(init.params) == 1) {  
        ii <- indexFinder(sfract, 0.5)
        if (ii == -1) {
            warning("ERROR: no survival fraction data below the .5 level.\n
                    Cannot use the initial r s k estimator.  You must supply initial r s k estimates")
            return(-1)
        }
        #else rsk <- c(1/time[ii], 0.01, 0.1) 
        else   	slope <- (sf[ii]-sf[ii-1]) /(tt[ii]-tt[ii-1])
        t50 <- tt[ii] + (0.5-sf[ii])/slope
        nslope <- slope*t50
        
        # Script for setting r.s.slope - the data frame used by function rsk.init
        # to produce initial r,s estimates from an estimated slope of the survival 
        # curve at the inflection point.
        c1<-c(0.000,0.100,0.200,0.300,0.400,0.500,0.600,0.700,0.800,0.900,0.930,0.950,0.960,0.970,
             0.980,0.990,0.991,0.992,0.993,0.994,0.995,0.996,0.997,0.998)
        c2<-c(1.48260200,1.40208800,1.31746800,1.22796200,1.13249500,1.02951700,0.91664920,0.78987300,0.64132740,
             0.45059940,0.37620130,0.31747940,0.28374720,0.24554230,0.20032630,0.14153800,0.13426380,0.12657460,
             0.11839010,0.10959890,0.10004140,0.08947239,0.07747895,0.06325606)
        c3<-c(-0.2143371,-0.2315594,-0.2518273,-0.2761605,-0.3061416,-0.3443959,-0.3956923,-0.4699249,-0.5925325,
             -0.8638228,-1.0422494,-1.2411051,-1.3920767,-1.6126579,-1.9815621,-2.8115970,-2.9646636,-3.1455461,
             -3.3638416,-3.6345703,-3.9827944,-4.4543775,-5.1451825,-6.3036320)
        r.s.slope<-data.frame(c1,c2,c3)
        dimnames(r.s.slope)[[2]]<-c("r","s","slope")
        rm(c1)
        rm(c2)
        rm(c3)
        
        rp<-r.s.slope[,1]
        sp<-r.s.slope[,2]
        slope.p<-r.s.slope[,3]
        
        sze<-length(slope.p)
        if(slope.p[sze] < nslope && nslope < slope.p[1])  { # check if normalized slope (nslope) is on the chart.
          for (i in 2:sze) {       
            if ( (slope.p[i-1] - nslope)*(slope.p[i] - nslope) < 0.0) {
              rri<-( rp[i-1] + (rp[i]-rp[i-1])*(nslope-slope.p[i-1])/(slope.p[i]-slope.p[i-1]) )/t50
              ssi<-( sp[i-1] + (sp[i]-sp[i-1])*(nslope-slope.p[i-1])/(slope.p[i]-slope.p[i-1]) )/sqrt(t50)
              break
            }
          }
        } else {    
          if (nslope <= slope.p[sze]) {
            rri<-rp[sze]/t50
            ssi<-sp[sze]/sqrt(t50)
          }
          else {
            rri<-rp[1]/t50
            ssi<-sp[1]/sqrt(t50)
          }
        }
        
        ssi<-ssi/1.1 # ssi was consistently overestimated above.
        
        #  --estimate initial k---
        
        #use rri,ssi and a data point (tt[ii],sf[ii]) to solve for kki 
        #      ..using the actual survival function.
        ii<-indexFinder(sf,.94)-1
        if (ii <= 1) {     # In case no data points between sruv=1 and surv=.94.
          ii<-2
          warning(message="WARNING: Initial time step may be too long.")
        }
        kki <- -(1/tt[ii])*log(sf[ii]/SurvFn.k(tt[ii],rri,ssi,0))
        if (kki <= 0) {
          kki <- (1.0 - sf[ii])/tt[ii] 
        }
        rsk <- c(rri, ssi, kki)
    } else {                      # use user specified init params
        rsk <- init.params
    }
    if (rsk[1] == -1) {
        stop
    }
    if (silent == FALSE) {
        print(cbind(c("Initial r", "initial s", "initial k"), rsk))
    }
    
    #  --create dataframe for sa---
    dtfm <- data.frame(x1 = x1, x2 = x2, Ni = Ni)
    
    
    #  --run MLE fitting routine---   
    # --conduct Newton-Ralphoson algorithm directly --
    fit.nlm <- nlminb(start = rsk, objective = logLikelihood.k, lower = lower, upper = upper, xx1 = x1, xx2 = x2, NNi = Ni)
    # -- if k<0 run again with k=0 --
#     if(fit.nlm$par[3]<0){
#       #k.final <- 0
#       warning("WARNING: k<0 on initial run.  Trying again with k=0.")
# #       fit.nlm <- nlminb(start = rsk, objective = logLikelihood.k, lower = c(lower[1:2], 0), upper = c(upper[1:2],0), xx1 = x1, xx2 = x2, NNi = Ni)
#     }
    
    # --save final param estimates---
    r.final <- fit.nlm$par[1]
    s.final <- abs(fit.nlm$par[2])
    k.final <- fit.nlm$par[3]
    mlv <- fit.nlm$obj  
    if (silent == FALSE) {print(cbind(c("estimated r", "estimated s", "estimated k", "minimum -loglikelihood value"),
                                      c(r.final, s.final, k.final, mlv)))}
    
    #  == end MLE fitting == = 
    
    #  --compute standard errors---
    if (se != FALSE) {
        s.e. <- stdErr.k(r.final, s.final, k.final, x1, x2, Ni, se)  
        if (silent == FALSE){print(cbind(c("sd for r", "sd for s", "sd for k"), s.e.))}        
    }
    
    #  --plotting and goodness of fit---
if (pplot != F || gfit != F) {
  plotting.k(r.final,s.final,k.final,mlv,time,sfract,x1,x2,Ni,pplot,tlab,lplot,cplot,Iplot,gfit,rc.data)
}
    # ............................................................................................
    
    # --return final param values---
    sigd <- 5   #significant digits of output
    if(se != F){
        params <- c(r.final, s.final, k.final)
        pvalue <- c(1-pnorm(r.final/s.e.[1]), 1-pnorm(s.final/s.e.[2]), 1-pnorm(k.final/s.e.[3]))
        std <- c(s.e.[1], s.e.[2], s.e.[3])
        out <- signif(cbind(params, std, pvalue), sigd)
        return(out)
    } else {
        return(signif(c(r.final, s.final, k.final)))
    }
}


#' Plotting function for 2-process vitality model. 4-param
#' 
#' None.
#' 
#' @param r.final r estimate
#' @param s.final s estimate
#' @param k.final k estimate
#' @param mlv TODO mlv
#' @param time time vector
#' @param sfract survival fraction
#' @param x1 Time 1
#' @param x2 Time 2
#' @param Ni Initial population
#' @param pplot Boolean. Plot cumulative survival fraction?
#' @param Iplot Boolean. Plot incremental survival?
#' @param Mplot Boolean. Plot mortality rate? 
#' @param tlab Character, label for time axis
#' @param rc.data Booolean, right-censored data?
plotting.k <- function(r.final,s.final,k.final,mlv,time,sfract,x1,x2,Ni,pplot,tlab,lplot,cplot,Iplot,gfit,rc.data){
  #  Function to provide plotting and goodness of fit computations
  #
  
  # --plot cumulative survival---
  if (pplot != F) {
    ext<-max(pplot,1)
    par(mfrow=c(1,1))
    len<-length(time)
    tmax <-ext*time[len]
    plot(time,sfract,xlab=tlab,ylab="survival fraction",ylim=c(0,1),xlim=c(0,tmax))
    xxx<-seq(0,tmax,length=200)
    lines(xxx,SurvFn.k(xxx,r.final,s.final,k.final))
    title("Cumulative Survival Data and Vitality Model Fitting")
  }
  
    # --likelihood and likelihood contour plots---
    if(lplot != F) {
      profilePlot <- function(r.f,s.f,k.f,x1,x2,Ni,mlv,cplot){

        SLL <- function(r,s,k,x1,x2,Ni){sum(logLikelihood.k(c(r,s,k),x1,x2,Ni))}
        
        rf<-.2; sf<-.5; kf<-1.0; fp<-40 # rf,sf,kf - set profile plot range (.2 => plot +-20%), 2*fp+1 points
        rseq <-seq((1-rf)*r.f,(1+rf)*r.f, (rf/fp)*r.f)
        sseq <-seq((1-sf)*s.f,(1+sf)*s.f, (sf/fp)*s.f)
        if (k.f > 0) {
          kseq <-seq((1-kf)*k.f,(1+kf)*k.f, (kf/fp)*k.f)
        } else {  #if k=0..
          kseq <-seq(.00000001,.1,length=(2*fp+1))
        }
        
        rl <-length(rseq)
        tmpLLr <-rep(0,rl)
        tmpLLs <-tmpLLr
        tmpLLk <-tmpLLr
        for (i in 1:rl) {
          tmpLLr[i] <-SLL(rseq[i],s.f,k.f,x1,x2,Ni)
          tmpLLs[i] <-SLL(r.f,sseq[i],k.f,x1,x2,Ni)	
          tmpLLk[i] <-SLL(r.f,s.f,kseq[i],x1,x2,Ni)
        }
        
        par(mfrow=c(1,3))
        
        rlim1 <-rseq[1]
        rlim2 <-rseq[rl]
        if (r.f < 0) {  #even though r should not be <0
          rlim2 <-rseq[1]
          rlim1 <-rseq[rl]
        }
        
        plot(r.f,LL<-SLL(r.f,s.f,k.f,x1,x2,Ni),
             xlim=c(rlim1,rlim2), xlab="r",ylab="Likelihood");  
        lines(rseq,tmpLLr)
        legend(x="topright", legend=c("r.final", "Likelihood varying r"), pch=c(1,NA), lty=c(NA, 1))
        plot(s.f,LL,xlim=c(sseq[1],sseq[rl]), xlab="s",ylab="Likelihood");  
        lines(sseq,tmpLLs)
        legend(x="topright", legend=c("s.final", "Likelihood varying s"), pch=c(1,NA), lty=c(NA, 1))
        title("Likelihood Plots")
        plot(k.f,LL,xlim=c(kseq[1],kseq[rl]), 
             ylim=c(1.1*min(tmpLLk)-.1*(mLk<-max(tmpLLk)),mLk),xlab="k",ylab="Likelihood"); 
        lines(kseq,tmpLLk)
        legend(x="topright", legend=c("k.final", "Likelihood varying k"), pch=c(1,NA), lty=c(NA, 1))
        
        
        # --  for contour plotting  ----------------------------------------
        
        if (cplot==T) {
          rl2<-(rl+1)/2; rl4<-20; st<-rl2-rl4; nr<-2*rl4+1
          tmpLLrs <-matrix(rep(0,nr*nr),nrow=nr,ncol=nr)
          for(i in 1:nr) {
            for(j in 1:nr) {
              tmpLLrs[i,j] <-SLL(rseq[i+st-1],sseq[j+st-1],k.f,x1,x2,Ni)
            }
          }
          
          lvv<-seq(mlv,1.02*mlv,length=11)  #99.8%, 99.6% ... 98%
          par(mfrow=c(1,1))
          contour(rseq[st:(rl2+rl4)],sseq[st:(rl2+rl4)],tmpLLrs,levels=lvv,xlab="r",ylab="s")
          title("Likelihood Contour Plot of r and s",  "Outermost ring is likelihhod 98% (of max) level, innermost is 99.8% level.")
          points(r.f,s.f,pch="*",cex=3.0)
          points(c(rseq[st],r.f),c(s.f,sseq[st]),pch="+",cex=1.5)
        }
        
      }
      profilePlot(r.final,s.final,k.final,x1,x2,Ni,mlv,cplot)
    }
    
  # --calculations for goodness of fit---
  if(gfit!=F) { # then gfit must supply the population number
    isp <-survProbInc.k(r.final,s.final,k.final,x1,x2)
    C1.calc<-function(pop, isp, Ni){
        #  Routine to calculate goodness of fit (Pearson's C -type test)
        #  pop - population number of sample
        #  isp  - Modeled: incemental survivor Probability 
        #  Ni   - Data:  incemental survivor fraction (prob)
        #
        # Returns:  a list containing:
        # C1, dof, Chi2    (retrieve each from list as ..$C1   etc.)
        
        if(pop<35){
          if(pop<25){
            warning(paste("WARNING: sample population (",as.character(pop),") is too small for 
                          meaningful goodness of fit measure.  Goodness of fit not being computed"))
            return()
          } else {
            warning(paste("WARNING: sample population (",as.character(pop),") may be too small for 
                          meaningful goodness of fit measure"))
          }
          }
        np <- pop * isp # modeled population at each survival probability level.
        tmpC1 <- 0
        i1<-1; i<-1; cnt<-0
        len <- length(np)
        while(i <= len) {
          idx <- i1:i
          # It is recommended that each np[i] >5 for meningful results.  Where np<5
          #   points are grouped to attain that level.   I have fudged it to 4.5 ... 
          #  (as some leeway is allowed, and exact populations are sometimes unknown).
          if(sum(np[idx]) > 4.5) {
            cnt <-cnt+1
            # Check if enough points remain.  If not, they are glommed onto previous grouping.
            if(i < len && sum(np[(i + 1):len]) < 4.5) {
              idx <- i1:len
              i <- len
            }
            
            sNi <- sum(Ni[idx])
            sisp <- sum(isp[idx])
            tmpC1 <- tmpC1 + (pop * (sNi - sisp)^2)/sisp
            i1 <-i+1			
          }
          
          i <-i+1
        }
        C1 <- tmpC1
        dof <-cnt-1-3    # degrees of freedom (3 is number of parameters).
        if (dof < 1) {
          warning(paste("WARNING: sample population (",as.character(pop),") is too small for 
		meaningful goodness of fit measure (DoF<1).  Goodness of fit not being computed"))
          return()
        }
        
        chi2<-qchisq(.95,dof)
        return(list(C1=C1,dof=dof,chi2=chi2))
        }
    C1dof <- C1.calc(gfit,isp,Ni)
    C1 <-C1dof$C1
    dof <-C1dof$dof
    chi2 <-C1dof$chi2
    
    print(paste("Pearson's C1=",as.character(round(C1,3)),"   chisquared =", as.character(round(chi2,3)),"on",as.character(dof),"degrees of freedom"))
    
    #   Note: The hypothesis being tested is whether the data could reasonably have come 
    #   from the assumed (vitality) model.
    if(C1 > chi2){
      print("C1 > chiSquared;   should reject the hypothesis becasue C1 falls outside the 95% confidence interval.")
    } else {
      print("C1 < chiSquared;   should Not reject the hypothesis because C1 falls inside the 95% confidence interval")		
    }
  }
  
  # --Incremental mortality plot
  if (Iplot != F) {
    par(mfrow=c(1,1))
    if (rc.data != F) {
      ln <-length(Ni)-1
      x1 <-x1[1:ln]
      x2 <-x2[1:ln]
      Ni <-Ni[1:ln]
    }
    ln <-length(Ni)
    #scale<-(x2-x1)[Ni==max(Ni)]
    scale<-max( (x2-x1)[Ni==max(Ni)] )
    
    ext<-max(pplot,1)
    
    npt<-200*ext
    xxx <-seq(x1[1],x2[ln]*ext,length=npt)
    xx1 <-xxx[1:(npt-1)]
    xx2 <-xxx[2:npt]
    sProbI <-survProbInc.k(r.final[1],s.final[1],k.final[1],xx1,xx2)
    
    ytop <-1.1*max( max(sProbI/(xx2-xx1)),Ni/(x2-x1) )*scale		
    plot((x1+x2)/2,Ni*scale/(x2-x1),ylim=c(0,ytop),xlim=c(0,ext*x2[ln]),xlab=tlab,ylab="incremental mortality")
    title("Probability Density Function")
    lines((xx1+xx2)/2,sProbI*scale/(xx2-xx1))
  }
  
  return()	 
}



#' The cumulative survival distribution function for 2-process 3-parameter
#' 
#' None.
#' 
#' @param xx vector of ages
#' @param r r value
#' @param s s value
#' @param k k value
#' @return vector of FF?
SurvFn.k <- function(xx, r, s, k) {
  # pnorm is: cumulative prob for the Normal Dist.
  tmp1 <- (sqrt(1/xx) * (1 - xx * r))/s    #  xx=0 is ok.  pnorm(+-Inf) is defined
  tmp2 <- (sqrt(1/xx) * (1 + xx * r))/s
  
  # --safeguard if exponent gets too large.---
  tmp3 <- 2*r/(s*s)
  
  if(tmp3 >250){   
    q <- tmp3/250 
    
    if(tmp3 >1500){
      q <- tmp3/500
    }
    
    valueFF <-(1.-(pnorm(-tmp1) + (exp(tmp3/q) *pnorm(-tmp2)^(1/q))^(q)))*exp(-k*xx) 
  }		    
  else {
    valueFF <-(1.-(pnorm(-tmp1) + exp(tmp3) *pnorm(-tmp2)))*exp(-k*xx)   #1-G
  }
  if ( all(is.infinite(valueFF)) ) {
    warning(message="Inelegant exit caused by overflow in evaluation of survival function.  
		Check for right-censored data. Try other initial values.")
  }
    
    return(valueFF)	
}


#' Calculates incremental survival probability for 2-process 3-parameter r, s, k
#' 
#' None
#' 
#' @param r r value
#' @param s s value
#' @param k k value
#' @param xx1 xx1 vector
#' @param xx2 xx2 vector
#' @return Incremental survival probabilities.
survProbInc.k <- function(r, s, k, xx1, xx2){
    value.iSP <- -(SurvFn.k(xx2, r, s, k) - SurvFn.k(xx1, r, s, k))
    value.iSP[value.iSP < 1e-18] <- 1e-18   # safeguards against taking Log(0)
    value.iSP
}


#' Gives log likelihood of 2-process 4-parameter model
#' 
#' None
#' 
#' @param par vector of parameter(r, s, lambda, beta)
#' @param xx1 xx1 vector
#' @param xx2 xx2 vector
#' @param NNi survival fractions
#' @return log likelihood
logLikelihood.k <- function(par, xx1, xx2, NNi) {
    # --calculate incremental survival probability--- (safeguraded > 1e-18 to prevent log(0))
    iSP <-  survProbInc.k(par[1], par[2], par[3], xx1, xx2)
    loglklhd <- -NNi*log(iSP)
    if (par[3] < 0) {
      loglklhd<-loglklhd + par[3]*par[3]*1e4
    }
    return(sum(loglklhd)) ## remove sum()? 
}


#' Standard errors for 3-param r, s, k
#' 
#' Note: if k <= 0, can not find std Err for k.
#' 
#' @param r r value
#' @param s s value
#' @param k k value
#' @param xx1 age 1 (corresponding 1:(t-1) and 2:t)
#' @param x2 age 2
#' @param Ni survival fraction
#' @param pop initial population (total population of the study)
#' @return standard error for r, s, k, u.
stdErr.k <- function(r,s,k,x1,x2,Ni,pop){  
  # function to compute standard error for MLE parameters r,s,k in the vitality model
  # Arguments:
  #  r,s,k - final values of MLE parameters
  #  xx1,xx2  time vectors (steps 1:(T-1) and 2:T)
  #  Ni - survival fraction
  #  pop - total population of the study
  #
  # Return:
  #  standard error for r,s,k
  # Note: if k <or= 0, can not find std Err for k.
  #
  LL <-function(a,b,c,r,s,k,x1,x2,Ni){sum(logLikelihood.k(c(r+a,s+b,k+c),x1,x2,Ni))}
  
  #initialize hessian for storage
  if (k > 0) {
    hess <- matrix(0,nrow=3,ncol=3)
  } else {
    hess <- matrix(0,nrow=2,ncol=2)
  } 
  
  #set finite difference intervals
  h <-.001
  hr <-abs(h*r)
  hs <-h*s*.1
  hk <-h*k*.1
  
  
  
  #Compute second derivitives (using 5 point)
  # LLrr
  f0 <-LL(-2*hr,0,0,r,s,k,x1,x2,Ni)
  f1 <-LL(-hr,0,0,r,s,k,x1,x2,Ni)
  f2 <-LL(0,0,0,r,s,k,x1,x2,Ni)
  f3 <-LL(hr,0,0,r,s,k,x1,x2,Ni)
  f4 <-LL(2*hr,0,0,r,s,k,x1,x2,Ni)
  
  fp0 <-(-25*f0 +48*f1 -36*f2 +16*f3 -3*f4)/(12*hr)
  fp1 <-(-3*f0 -10*f1 +18*f2 -6*f3 +f4)/(12*hr)
  fp3 <-(-f0 +6*f1 -18*f2 +10*f3 +3*f4)/(12*hr)
  fp4 <-(3*f0 -16*f1 +36*f2 -48*f3 +25*f4)/(12*hr)
  
  LLrr <-(fp0 -8*fp1 +8*fp3 -fp4)/(12*hr)
  
  # LLss
  f0 <-LL(0,-2*hs,0,r,s,k,x1,x2,Ni)
  f1 <-LL(0,-hs,0,r,s,k,x1,x2,Ni)
  # f2 as above
  f3 <-LL(0,hs,0,r,s,k,x1,x2,Ni)
  f4 <-LL(0,2*hs,0,r,s,k,x1,x2,Ni)
  
  fp0 <-(-25*f0 +48*f1 -36*f2 +16*f3 -3*f4)/(12*hs)
  fp1 <-(-3*f0 -10*f1 +18*f2 -6*f3 +f4)/(12*hs)
  fp3 <-(-f0 +6*f1 -18*f2 +10*f3 +3*f4)/(12*hs)
  fp4 <-(3*f0 -16*f1 +36*f2 -48*f3 +25*f4)/(12*hs)
  
  LLss <-(fp0 -8*fp1 +8*fp3 -fp4)/(12*hs)
  
  # LLkk
  if (k > 0) {
    f0 <-LL(0,0,-2*hk,r,s,k,x1,x2,Ni)
    f1 <-LL(0,0,-hk,r,s,k,x1,x2,Ni)
    # f2 as above
    f3 <-LL(0,0,hk,r,s,k,x1,x2,Ni)
    f4 <-LL(0,0,2*hk,r,s,k,x1,x2,Ni)
    
    fp0 <-(-25*f0 +48*f1 -36*f2 +16*f3 -3*f4)/(12*hk)
    fp1 <-(-3*f0 -10*f1 +18*f2 -6*f3 +f4)/(12*hk)
    fp3 <-(-f0 +6*f1 -18*f2 +10*f3 +3*f4)/(12*hk)
    fp4 <-(3*f0 -16*f1 +36*f2 -48*f3 +25*f4)/(12*hk)
    
    LLkk <-(fp0 -8*fp1 +8*fp3 -fp4)/(12*hk)
  }
  
  
  #-------end second derivs---
  # do mixed partials (4 points)
  # LLrs
  m1 <-LL(hr,hs,0,r,s,k,x1,x2,Ni)
  m2 <-LL(-hr,hs,0,r,s,k,x1,x2,Ni)
  m3 <-LL(-hr,-hs,0,r,s,k,x1,x2,Ni)
  m4 <-LL(hr,-hs,0,r,s,k,x1,x2,Ni)
  
  LLrs <-(m1 -m2 +m3 -m4)/(4*hr*hs)
  
  if (k > 0) {
    # LLrk
    m1 <-LL(hr,0,hk,r,s,k,x1,x2,Ni)
    m2 <-LL(-hr,0,hk,r,s,k,x1,x2,Ni)
    m3 <-LL(-hr,0,-hk,r,s,k,x1,x2,Ni)
    m4 <-LL(hr,0,-hk,r,s,k,x1,x2,Ni)
    
    LLrk <-(m1 -m2 +m3 -m4)/(4*hr*hk)
    
    # LLsk
    m1 <-LL(0,hs,hk,r,s,k,x1,x2,Ni)
    m2 <-LL(0,-hs,hk,r,s,k,x1,x2,Ni)
    m3 <-LL(0,-hs,-hk,r,s,k,x1,x2,Ni)
    m4 <-LL(0,hs,-hk,r,s,k,x1,x2,Ni)
    
    LLsk <-(m1 -m2 +m3 -m4)/(4*hs*hk)
  }
  
  if (k > 0) {
    diag(hess) <-c(LLrr,LLss,LLkk)*pop
    hess[2,1]<-hess[1,2]<-LLrs*pop
    hess[3,1]<-hess[1,3]<-LLrk*pop
    hess[3,2]<-hess[2,3]<-LLrk*pop
  } else {
    diag(hess) <-c(LLrr,LLss)*pop
    hess[2,1]<-hess[1,2]<-LLrs*pop
  }
  #print(hess)
  hessInv <-solve(hess)
  #print(hessInv)
  
  #compute correlation matrix:
  sz <-3
  if (k <= 0) { sz <-2 }
  corr<- matrix(0,nrow=sz,ncol=sz)
  for (i in 1:sz) {
    for (j in 1:sz) {
      corr[i,j] <- hessInv[i,j]/sqrt(abs(hessInv[i,i]*hessInv[j,j]))
    }
  }
  #print(corr)
  if ( abs(corr[2,1]) > .98 ) {
    warning("WARNING: parameters r and s appear to be closely correlated for this data set.  
            s.e. may fail for these parameters.")
  }
  if (  sz == 3 && abs(corr[3,2]) > .98 ) {
    warning("WARNING: parameters s and k appear to be closely correlated for this data set.  
            s.e. may fail for these parameters.")
  }
  if (  sz == 3 && abs(corr[3,1]) > .98 ) {
    warning("WARNING: parameters r and k appear to be closely correlated for this data set.  
            s.e. may fail for these parameters.")
  }
  
  
  se <-sqrt(diag(hessInv))
  
  #  Approximate s.e. for cases where calculation of s.e. failed:
  if( sum( seNA<-is.na(se) ) > 0 ) {
    se12 <-sqrt(diag(solve(hess[c(1,2)	,c(1,2) ])))
    if (k > 0) {
      se13 <-sqrt(diag(solve(hess[c(1,3)	,c(1,3) ])))
      se23 <-sqrt(diag(solve(hess[c(2,3)	,c(2,3) ])))
    }
    
    if(seNA[1] == T) {
      if(!is.na(se[1]<-se12[1]) || (k>0 && !is.na(se[1]<-se13[1])) )
        warning("* s.e. for parameter r is approximate.")
      else warning("* unable to calculate or approximate s.e. for parameter r.")
    }
    if(seNA[2] == T) {
      if(!is.na(se[2]<-se12[2]) || (k>0 && !is.na(se[2]<-se23[1])) )
        warning("* s.e. for parameter s is approximate.")
      else warning("* unable to calculate or approximate s.e. for parameter s.")
    }
    if(k>0 && seNA[3] == T) {
      if(!is.na(se[3]<-se13[2]) ||  !is.na(se[3]<-se23[2]) )
        warning("* s.e. for parameter k is approximate.")
      else warning("* unable to calculate or approximate s.e. for parameter k.")
    }
    
     
  }
  
  
  #######################
  if (k <= 0) {
    se <-c(se,NA)
  }
  return(se)
  }
