#' Animal activity statistics
#'
#' Provides functions to estimate and compare activity parameters from sensor data.
#'  
#' @details Sensors that record active animals (eg camera traps) build up a record of 
#' the distribution of activity over the course of the day. Records are more frequent 
#' when animals are more active, and less frequent or absent when animals are inactive. 
#' The area under the distribution of records thus contains information on the overall 
#' level of activity in a sampled population. This package provides tools for plotting 
#' activity distributions, quantifying the overall level of activity with error, and 
#' statistically comparing distributions through bootstrapping.
#'  
#' The core function is \code{fitact}, which creates an \code{actmod} object containing 
#' the circular kernel PDF, and the activity level estimate derived from this. The 
#' generic plot function for \code{actmod} objects plots the distribution. Functions 
#' starting with \code{compare} make statistical comparisons between distributions or 
#' activity estimates. Note that all time or other circular data should be in radians 
#' (in the range 0 to 2*pi).
#'  
#' @references Rowcliffe, M., Kays, R., Kranstauber, B., Carbone, C., Jansen, P.A. (2014) Quantifying animal activity level using camera trap data. Methods in Ecology and Evolution.
#' @seealso \code{\link{overlap}}
#' @docType package
#' @name activity
NULL

#' Animal record time of day data
#' 
#' Barro Colorado Island 2008 data: times of day at which animal records occured 
#' (\code{time}), together with species (\code{species}).
#' 
#' @format A dataframe with 17820 observations and 2 variables.
#' @source http://dx.doi.org/10.6084/m9.figshare.1160536
#' @name BCItime
#' @docType data
NULL

#' Animal speed data
#' 
#' Barro Colorado Island 2008 data: speeds of animal passages past camera traps 
#' (\code{speed}), together with species (\code{species}) and time of day (\code{time}) 
#' for each record.
#' 
#' @format A dataframe with 2204 observations and 3 variables.
#' @source http://dx.doi.org/10.6084/m9.figshare.1160536
#' @name BCIspeed
#' @docType data
NULL

#' Activity model class.
#' 
#' An S4 class describing activity models fitted to time of observation data.
#' 
#' @slot data Object of class \code{"numeric"}, the input data.
#' @slot wt Object of class \code{"numeric"}, weights applied to the data.
#' @slot bw Object of class \code{"numeric"}, kernel bandwidth.
#' @slot adj Object of class \code{"numeric"}, kernel bandwidth adjustment multiplier.
#' @slot pdf Object of class \code{"matrix"} describing fitted probability density function: 
#'  Column 1: Sequence of radian times at which PDF evaluated (specifically seq(0,2*pi,pi/256)). 
#'  Column 2: Corresponding circular kernel PDF values. 
#' Additionally if errors boostrapped:
#'  Column 3: PDF standard error. 
#'  Column 4: PDF lower confidence interval. Column 5: PDF upper confidence interval.
#' @slot act Object of class \code{"numeric"} giving activity level estimate and, if errors boostrapped, standard error and 95 percent confidence limits. 
#' @method plot \code{signature(x = "actmod")}: Plots PDF, confidence interval if calculated, and optionally data distribution.
#' @export
setClass("actmod", 
         representation(data="numeric", wt="numeric", bw="numeric", adj="numeric", 
                        pdf="matrix", act="numeric"))

#' An S4 class describing linear-circular relationships.
#' 
#' @slot data Object of class \code{"data.frame"}, the input data, with columns 
#' \code{lindat} (linear data) and \code{circdat} (circular data).
#' @slot fit Object of class \code{"data.frame"}, summary of the model fit, with columns:
#'  \code{x}: A regular ascending sequence from 0 to 2*pi at which other columns evaluated;
#'  \code{fit}: The linear fitted values;
#'  \code{p}: The two tailed probability of observing the fitted values under a random (null) circular distribution;
#'  \code{nullLCL}: The lower confidence limit of the null distribution;
#'  \code{nullUCL}: The upper confidence limit of the null distribution.
#' @method plot \code{signature(x = "actmod")}: Plots PDF, confidence interval if calculated, and optionally data distribution.
#' @export
setClass("lincircmod", representation(data="data.frame", fit="data.frame"))

#' Impute empirical circular distribution.
#'
#' Imputes values at given points on an empirical circular distribution.
#' 
#' Note that x is assumed circular, so first and last \code{y} values should be equal. Evaluation points \code{xvals} should also be within the range of \code{x}.
#' 
#' @param xvals Numeric circular values at which to evaluate \code{y}
#' @param x Evenly spaced ascending numeric sequence of circular values
#' @param y Empirical numeric output distribution matched with \code{x}
#' @return A numeric vector of \code{y} values evaluated at \code{xvals}
#' @examples
#' #Abstract example
#' x <- seq(0,2*pi,length.out=11)
#' y <- c(0,1,2,3,4,5,4,3,2,1,0)
#' yfromx(0:6,x,y)
#' 
#' #BCI data example
#' #Weighting ocelot activity pattern to correct for variation in speed
#' data(BCIspeed)
#' data(BCItime)
#' #Fit linear-circular model to log(speed)
#' i <- BCIspeed$species=="ocelot"
#' lcfit <- fitlincirc(BCIspeed$time[i]*2*pi, log(BCIspeed$speed[i]), reps=50)
#' #Fit weighted activity model using yfromx to create weights
#' j <- BCItime$species=="ocelot"
#' tdat <- BCItime$time[j]*2*pi
#' w <- 1/yfromx(tdat, lcfit@@fit$x, exp(lcfit@@fit$fit))
#' mod <- fitact(tdat, wt=w, sample="none")
#' plot(mod)
#' #Oveplot unweighted model for comparison
#' mod2 <- fitact(tdat, sample="none")
#' plot(mod2, lcol=3, add=TRUE)
#' @export
yfromx <- function(xvals, x, y)
{ if(length(x)!=length(y)) stop("x and y lengths are unequal")
  if(max(diff(x))-min(diff(x))>1e-10) stop("x is not an evenly spaced ascending sequence")
  if(abs(y[1]-y[length(y)])>1e-10) stop("first and last y values should be equal")
  if(min(xvals)<min(x) | max(xvals)>max(x)) stop("some xvals are not within the range of x")
  
  unit <- x[2]
  xvals[xvals==0] <- max(x)
  i <- ceiling(xvals/unit)
  y[i] + (y[i+1]-y[i])*(xvals-x[i])/unit
}

#' Calculate circular kernel bandwidth.
#' 
#' Uses an optimisation procedure to calculate the circular kernel bandwidth giving the best fit to the data.
#' 
#' Mainly for internal use.
#' 
#' @param dat Numeric data vector of radian times.
#' @param K Integer number of values of kappa over which to maximise (see references for details).
#' @return Single numeric bandwidth value.
#' @references Ridout, M.S. & Linkie, M. (2009) Estimating overlap of daily activity patterns from camera trap data. Journal of Agricultural Biological and Environmental Statistics, 14, 322-337.
#' @export
bwcalc <- function(dat,K=3)
{  if(!all(dat>=0 & dat<=2*pi)) warning("some dat values are <0 or >2*pi, expecting radian data")
   if(max(dat)<1) warning("max(dat) < 1, expecting radian data")
   
   minfunc <- function(kap,k,dat)
   {	trigmom <- circular::trigonometric.moment(circular::circular(dat),k,center=T)$rho
     (besselI(kap,k)/besselI(kap,0) - trigmom)^2
   }
   kapk.calc <- function(k,dat)
     optimise(minfunc,c(0,100),k,dat)$minimum
   kap <- max(sapply(1:K, kapk.calc, dat))
   ((3*length(dat)*kap^2*besselI(2*kap,2)) / (4*pi^0.5*besselI(kap,0)^2))^(2/5)
}

#'  Circular kernel probability density function.
#'  
#'  Optionally weighted Von Mises kernel probability densities.
#'  
#'  If \code{bw} not provided it is calculated internally using \code{bw.calc}. The \code{adj} argument is used to adjust \code{bw} to facilitate exploration of fit flexibility.
#'  
#'  @param x Numeric vector of radian times at which to evaluate the PDF.
#'  @param dat Numeric vector of radian time data to which the PDF is fitted.
#'  @param wt A numeric vector of weights for each \code{dat} value.
#'  @param bw Numeric value for kernel bandwidth.
#'  @param adj Numeric kernel bandwidth multiplier.
#'  @return Numeric vector of probability densities evaluated at \code{x}.
#'  @seealso \code{\link{density.circular}, \link{bwcalc}}
#'  @examples
#'  #Example with made up input
#'  tt <- runif(100,0,2*pi)
#'  xx <- seq(0,2*pi, pi/256)
#'  pdf <- dvmkern(xx, tt)
#'  plot(xx, pdf, type="l")
#'  @export
dvmkern <- function(x,dat,wt=NULL,bw=NULL,adj=1)
{  if(!all(dat>=0 & dat<=2*pi)) warning("some dat values are <0 or >2*pi, expecting radian data")
   if(max(dat)<1) warning("max(dat) < 1, expecting radian data")
   if(!all(x>=0 & x<=2*pi)) warning("some x values are <0 or >2*pi, expecting radian values")
   if(!is.null(wt) & length(wt)!=length(dat)) stop("dat and wt have different lengths")
   
   if(is.null(bw)) bw <- bwcalc(dat)
   if(is.null(wt)) wt <- rep(1,length(dat))
   dx <- expand.grid(dat,x)
   dif <- abs(dx[,2]-dx[,1])
   i <- dif>pi
   dif[i] <- 2*pi-dif[i]
   prob <- circular::dvonmises(circular::circular(dif),circular::circular(0),bw*adj)
   apply(matrix(prob*wt, nrow=length(dat)),2,sum)/sum(wt)
}

#'  Random circular kernel numbers.
#'  
#'  Random numbers drawn from a fitted Von Mises kernel distribution.
#'  
#'  @details If \code{fit} is a matrix, it should have two columns:
#'   [,1] precisely seq(0,2*pi,pi/256) sequence of radian values at which pdf evaluated; 
#'   [,2] corresponding pdf values.
#'  @param n Integer number of numbers to return.
#'  @param fit An emprical distribution contained in either a matrix or a fitted \code{actmod} object (see Details).
#'  @return A numeric vector of radian numbers.
#'  @examples
#'  #Matrix input
#'  data(BCItime)
#'  tm <- 2*pi*BCItime$time[BCItime$species=="paca"]
#'  xx <- seq(0,2*pi, pi/256)
#'  rn <- rvmkern(100, cbind(xx, dvmkern(xx, tm)))
#'  
#'  #actmod input
#'  fit <- fitact(tm, sample="n")
#'  rvmkern(100,fit)
#'  @export
rvmkern <- function(n,fit)
{ fclass <- class(fit)
  if(fclass=="actmod") fit <- fit@pdf else
    if(fclass=="matrix")
    { if(nrow(fit)!=513 | ncol(fit)!=2) 
      stop("fit matrix should have exactly 513 rows and 2 columns")
      if(max(abs(seq(0,2*pi,pi/256)-fit[,1]))>0)
        stop("first column of fit matrix should be seq(0,2*pi,pi/256)")
      if(abs(1-sum(fit[,2])*2*pi/512)>0.01) 
        warning("second column of fit matrix doesn't look like a pdf and was rescaled to sum to 1")
    } else
      stop("fit should be of either matrix or actmod class")
  cumpdf <- c(0,cumsum(fit[-513,2]*2*pi/512))
  cumpdf <- cumpdf/cumpdf[513]
  rn <- runif(n)
  i <- findInterval(rn,cumpdf)
  fit[i,1] + ((rn-cumpdf[i])/(cumpdf[i+1]-cumpdf[i]))*2*pi/512
}

#' Fit activity model to time-of-day data
#' 
#' Fits a circular kernel density to radian time-of-day data and estimates activity 
#' level from this distribution. Optionally bootstraps the distribution, in which 
#' case SEs and confidence limits are also stored for activity level and PDF.
#' 
#' @details The bandwidth adjustment multiplier \code{adj} is provided to allow 
#' exploration of the effect of adjusting the internally calculated bandwidth on 
#' accuracy of activity level estimates. The alternative bootstrapping methods 
#' defined by \code{sample} are:
#'  \code{data}: sample from the data;
#'  \code{model}: sample from the fitted probability density distribution;
#'  \code{none}: no bootstrapping.
#' Confidence interval coverage seems to be better at large sample size 
#' (greater than 100-200) using \code{"model"}, but better at small sample size 
#' when using \code{"data"}. The reason for this  needs further investigation.
#' @param dat A numeric vector of radian time-of-day data.
#' @param wt A numeric vector of weights for each \code{dat} value.
#' @param reps Number of boostrap iterations to perform. Ignored if sample=="none".
#' @param bw Numeric value for kernel bandwidth. If NULL, calculated internally.
#' @param adj Numeric bandwidth adjustment multiplier.
#' @param sample Character string defining sampling method for bootstrapping errors (see details).
#' @param show Logical whether or not to show a progress bar while bootstrapping.
#' @return An object of type \code{actmod}
#' @examples
#' #Fit without confidence limits
#' data(BCItime)
#' tdat <- 2*pi*BCItime$time[BCItime$species=="ocelot"]
#' mod1 <- fitact(tdat, sample="none")
#' plot(mod1)
#'  
#' #Fit with confidence limits (limited reps to speed up)
#' mod2 <- fitact(tdat, reps=10)
#' plot(mod2)
#'  
#' #Fit weighted function to correct for detection radius 1.21 times higher
#' #by day than by night, assuming day between pi/2 (6am) and pi*2/3 (6pm)
#' weight <- 1/ifelse(tdat>pi/2 & tdat<pi*3/2, 1.2, 1)
#' mod3 <- fitact(tdat, wt=weight, sample="none")
#' plot(mod3)
#' #Overplot unweighted version for comparison
#' plot(mod1, add=TRUE, lcol=3)
#' @export
fitact <- function(dat, wt=NULL, reps=1000, bw=NULL, adj=1, sample=c("data","model","none"), show=TRUE)
{ if(!all(dat>=0 & dat<=2*pi)) warning("some dat values are <0 or >2*pi, expecting radian data")
  if(max(dat)<1) warning("max(dat) < 1, expecting radian data")
  if(!is.null(wt) & length(wt)!=length(dat)) stop("dat and wt have different lengths")
  
  sample <- match.arg(sample)
  if(is.null(bw)) bw <- bwcalc(dat)
  x <- seq(0,2*pi,pi/256)
  pdf <- dvmkern(x, dat, wt, adj, bw)
  act <- 1/(2*pi*max(pdf))  
  if(sample=="none")
    sepdf <- lclpdf <- uclpdf <- seact <- lclact <- uclact <- numeric(0) else
    { if(sample=="model") samp <- matrix(rvmkern(reps*length(dat), cbind(x,pdf)), ncol=reps) else
      if(sample=="data") samp <- matrix(sample(dat, reps*length(dat), replace=TRUE), ncol=reps)
      if(show)
        pdfs <- pbapply::pbapply(samp, 2, function(dat) dvmkern(x,dat,wt,adj=adj)) else
          pdfs <- apply(samp, 2, function(dat) dvmkern(x,dat,wt,adj=adj))
      sepdf <- apply(pdfs,1,sd)
      lclpdf <- apply(pdfs,1,quantile,probs=0.025)
      uclpdf <- apply(pdfs,1,quantile,probs=0.975)
      acts <- 1/(2*pi*apply(pdfs,2,max))
      seact <- sd(acts)
      lclact <- quantile(acts,0.025)
      uclact <- quantile(acts,0.975)    
    }
  if(is.null(wt)) wt <- 1
  new("actmod", data=dat, wt=wt, bw=bw, adj=adj, 
      pdf=cbind(time=x, pdf=pdf, se=sepdf, lcl=lclpdf, ucl=uclpdf),
      act=c(act=act, se=seact, lcl=lclact, ucl=uclact))
}

#' Compare circular distributions.
#' 
#' Test probability that two sets of circular observations come from the same distribution.
#' 
#' Bootstrap test calculates overlap index 2 (see references) for the observed data samples, then generates a null distribution of overlap indices using data sampled randomly with replacement from the combined data. This randomised distribution is then used to estimate the probability that the observed overlap arose by chance.
#' 
#' @param y1,y2 Numeric vectors of radian data.
#' @param reps Number of bootstrap iterations.
#' @return A named 2-element vector: Overlap = observed overlap index; p = probability observed index arose by chance.
#' @references Ridout, M.S. & Linkie, M. (2009) Estimating overlap of daily activity patterns from camera trap data. Journal of Agricultural Biological and Environmental Statistics, 14, 322-337.
#' @seealso \code{\link{overlapEst}}
#' @examples
#' #Example with bootstrap reps limited to speed up
#' data(BCItime)
#' tPaca <- 2*pi*BCItime$time[BCItime$species=="paca"]
#' tRat <- 2*pi*BCItime$time[BCItime$species=="rat"]
#' compareCkern(tPaca,tRat,reps=10)
#' @export
compareCkern <- function(y1,y2,reps=1000)
{  ovl <- overlap::overlapEst(y1,y2)[2]
   names(ovl) <- NULL
   samp <- matrix(sample(c(y1,y2), reps*length(c(y1,y2)), replace=TRUE), ncol=reps)
   res <- pbapply::pbsapply(1:reps, function(i)
     overlap::overlapEst(samp[1:length(y1),i], samp[length(y1)+(1:length(y2)),i])[2])
   fun <- ecdf(res)
   c(Overlap=ovl, p=fun(ovl))
}

#' Compare activity level estimates
#' 
#' Wald test for the statistical difference between two or more activitiy level estimates.
#' 
#' Uses a Wald test to ask whether the difference between estimates a1 and a2 is 
#' significantly different from 0: statistic W = (a1-a2)^2 / (SE1^2+SE2^2) tested 
#' on chi-sq distribution with 1 degree of freedom.
#' 
#' @param fits A list of fitted \code{actmod} objects
#' @return A matrix with 4 columns: 1. differences between estimates; 2. SEs of the differences; 3. Wald statistics; 4. p-values (H0 is no difference between estimates). Matrix rows give all possible pairwise comparisons, numbered in the order in which they entered in the list \code{fits}.
#' @examples
#' #Test whether paca have a sigificantly different activity level from rat.
#' #Bootstrap reps limited to speed up example.
#' data(BCItime)
#' tPaca <- 2*pi*BCItime$time[BCItime$species=="paca"]
#' tRat <- 2*pi*BCItime$time[BCItime$species=="rat"]
#' (fPaca <- fitact(tPaca, reps=10))
#' (fRat <- fitact(tRat, reps=10))
#' compareAct(list(fPaca,fRat))
#' @export
compareAct <- function(fits)
{ if(class(fits)!="list" | !all(unlist(lapply(fits,class))=="actmod")) 
    stop("fits must be a list of actmod objects")
  if(min(unlist(lapply(fits, function(x) length(x@act))))==1)
    stop("all input model fits must be boostrapped")
  
  len <- length(fits)
  i <- rep(1:(len-1), (len-1):1)
  j <- unlist(sapply(2:len, function(i) i:len))
  acts <- unlist(lapply(fits, function(fit) fit@act[1]))
  seacts <- unlist(lapply(fits, function(fit) fit@act[2]))
  dif <- acts[i]-acts[j]
  vardif <- seacts[i]^2 + seacts[j]^2
  W <- dif^2/vardif
  prob <- 1-pchisq(W,1)
  res <- cbind(Difference=dif, SE=sqrt(vardif), W=W, p=prob)
  dimnames(res)[[1]] <- paste(i,j,sep="v")
  res
}

#' Compare activity across between times of day
#' 
#' Uses a Wald test to statistically compare activity levels at given radian times of day for a fitted activity distribution.
#' 
#' Bootrapping the activity model yields standard error estimates for the PDF. This function uses these SEs to compute a Wald statistic for the difference between PDF values (by inference activity levels) at given times of day: statistic W = (a1-a2)^2 / (SE1^2+SE2^2) tested on chi-sq distribution with 1 degree of freedom.
#' 
#' @param fit Fitted \code{actmod} object with errors boostrapped (fit using \code{fitact} with \code{sample} argument != "none").
#' @param times Numeric vector of radian times of day at which to compare activity levels. All pairwise comparisons are made.
#' @return A matrix with 4 columns: 1. differences between PDF values; 2. SEs of the differences; 3. Wald statistics; 4. p-values (H0 is no difference between estimates). Matrix rows give all possible pairwise comparisons, numbered in the order in which they appear in vector \code{times}.
#' @examples
#' data(BCItime)
#' tPaca <- 2*pi*BCItime$time[BCItime$species=="paca"]
#' fPaca <- fitact(tPaca, reps=10)
#' plot(fPaca, hrs=FALSE, frq=FALSE)
#' compareTimes(fPaca, c(5.5,6,0.5,1))
#' @export
compareTimes <- function(fit, times)
{ if(class(fit)!="actmod") stop("fit input must be an actmod object")
  if(!all(times>=0 & times<=2*pi)) stop("some times are <0 or >2*pi, expecting radian data")
  
  len <- length(times)
  i <- rep(1:(len-1), (len-1):1)
  j <- unlist(sapply(2:len, function(i) i:len))
  k <- findInterval(times, fit@pdf[,1])
  p <- (times-fit@pdf[,1][k]) / (fit@pdf[,1][k+1]-fit@pdf[,1][k])
  pdfs1 <- fit@pdf[,2][k]
  pdfs2 <- fit@pdf[,2][k+1]
  pdfs <- pdfs1 + p*(pdfs2-pdfs1)
  sepdfs1 <- fit@pdf[,3][k]
  sepdfs2 <- fit@pdf[,3][k+1]
  sepdfs <- sepdfs1 + p*(sepdfs2-sepdfs1)
  dif <- pdfs[i]-pdfs[j]
  vardif <- sepdfs[i]^2 + sepdfs[j]^2
  W <- dif^2/vardif
  prob <- 1-pchisq(W,1)
  res <- cbind(Difference=dif, SE=sqrt(vardif), W=W, p=prob)
  dimnames(res)[[1]] <- paste(i,j,sep="v")
  res
}
#' Show activity level estimate
#' 
#' Prints the \code{act} slot (activity level estimate) from an \code{actmod} object.
#' 
#' @param object Object of class \code{actmod}.
#' @export
setMethod("show", "actmod", function(object) print(object@act))

#' Plot activity distribution
#' 
#' Plot an activity probability distribution from a fitted \code{actmod} object.
#' 
#' @param x Object of class \code{actmod}.
#' @param hrs Logical defining whether to plot x axis as hours (default) or radians.
#' @param frq Logical defining whether to plot y axis as frequency (default) or probability density.
#' @param dat Data distribution plotting style, one of "h" (histogram, default), "r" (rug), "n" (no data plotting, PDF only).
#' @param add Logical defining whether to create a new plot (default) or add the probability density to an existing plot (in which case no data are plotted).
#' @param dcol Numeric or character defining colour of data lines.
#' @param lcol Numeric or character defining colour of PDF lines.
#' @param ... Additional plotting arguments passed to internal plot call. NB axis limits are calculated internally and cannot be reset.
#' @docType methods
#' @rdname plot-actmod-methods
#' @export
setMethod("plot", "actmod",
          function(x, hrs=TRUE, frq=TRUE, dat=c("histogram","rug","none"), add=FALSE, dcol=1, lcol=2, ...)
          {	dat <- match.arg(dat)
            fit <- x
            x <- fit@pdf[,1]
            y <- fit@pdf[,2]
            if(ncol(fit@pdf)==5)
            { lcl <- fit@pdf[,4]
              ucl <- fit@pdf[,5]
            } else lcl <- ucl <- numeric(0)
            data <- fit@data
            if(hrs)
            {	x <- x*12/pi
              data <- data*12/pi
              maxbrk <- 24
              xaxticks <- c(0,24,4)
            }else
            {	maxbrk <- 2*pi
              xaxticks <- NULL
            }
            h <- hist(data, breaks=seq(0,maxbrk,maxbrk/24), plot=F)
            if(frq) 
            {	y <- y*length(data)*pi/12 
              lcl <- lcl*length(data)*pi/12 
              ucl <- ucl*length(data)*pi/12 
              d <- h$counts
            }else
            {	if(hrs)
            { y <- y*pi/12
              lcl <- lcl*pi/12
              ucl <- ucl*pi/12
            }
            d <- h$density
            }
            if(!add) 
            {	time <- 0
              frequency <- 0
              plot(time,frequency, type="n", xaxp=xaxticks, 
                   xlim=c(0,max(x)), ylim=c(0,max(y,d,ucl)), ...)
              if(dat=="histogram")
                lines(h$breaks, c(d,d[1]), type="s", col=dcol) else
                  if(dat=="rug")
                    for(i in 1:length(data))
                      lines(rep(data[i],2),max(y)*-c(0.05,0.02), lwd=0.1, col=dcol)
            }
            lines(x, y, col=lcol)
            if(length(lcl)>0)
            { lines(x, lcl, col=lcol, lty=2)
              lines(x, ucl, col=lcol, lty=2)
            }
          }
)

#' Linear-circular kernel fit
#'
#' Fits a Von Mises kernel distribution describing a linear variable as a function of a circular predictor.
#'
#' @param x Numeric vector of radian values at which to evaluate the distribution.
#' @param circdat Numeric vector of radian data matched with \code{lindat}.
#' @param lindat Numeric vector of linear data matched with \code{circdat}.
#' @return A numeric vector of fitted \code{lindat} values matched with \code{x}.
#' @references Xu, H., Nichols, K. & Schoenberg, F.P. (2011) Directional kernel regression for wind and fire data. Forest Science, 57, 343-352.
#' @examples
#' data(BCIspeed)
#' i <- BCIspeed$species=="ocelot"
#' sp <- log(BCIspeed$speed[i])
#' tm <- BCIspeed$time[i]*2*pi
#' circseq <- seq(0,2*pi,pi/256)
#' trend <- lincircKern(circseq, tm, sp)
#' plot(circseq, trend, type="l")
#' @export
lincircKern <- function(x,circdat,lindat)
{ if(length(lindat)!=length(circdat)) 
  stop("lindat and circdat lengths are unequal")
  if(min(circdat)<0 | max(circdat)>2*pi) 
    stop("circdat values not between 0 and 2*pi, expecting radian data")
  if(min(x)<0 | max(x)>2*pi) 
    stop("x values not between 0 and 2*pi, expecting radian values")
  
  hs <- 1.06 * min(sd(lindat), (quantile(lindat,0.75)-quantile(lindat,0.25))/1.34) * 
    length(lindat)^-0.2
  bw <- 1/hs^2
  dx <- expand.grid(circdat,x)
  dif <- abs(dx[,2]-dx[,1])
  i <- dif>pi
  dif[i] <- 2*pi-dif[i]
  prob <- matrix(circular::dvonmises(circular::circular(dif),circular::circular(0),bw), nrow=length(circdat))
  apply(prob,2,function(z) mean(z*lindat)/mean(z))
}

#' Linear-circular regression
#'
#' Fits a Von Mises kernel distribution describing a linear variable as a function 
#' of a circular predictor, and boostraps the null distribution in order to evaluate 
#' significance of radial variation in the linear variable.
#'
#' Deviation of \code{lindat} from the null expecation is assessed either visually 
#' by the degree to which the fitted distribution departs from the null confidence 
#' interval (use generic plot function), or quantitatively by column \code{p} of 
#' slot \code{fit} in the resulting \code{lincircmod-class} object.
#'
#' @param circdat Numeric vector of radian data matched with \code{lindat}.
#' @param lindat Numeric vector of linear data matched with \code{circdat}.
#' @param pCI Single numeric value between 0 and 1 defining proportional confidence interval to return.
#' @param reps Integer number of bootstrap repetitions to perform.
#' @param res Resolution of fitted distribution and null confidence interval - specifically a single integer number of points on the circular scale at which to record distributions.
#' @return An object of type \code{\link{lincircmod-class}}
#' @references Xu, H., Nichols, K. & Schoenberg, F.P. (2011) Directional kernel regression for wind and fire data. Forest Science, 57, 343-352.
#' @examples
#' #Example with reps limited to increase speed
#' data(BCIspeed)
#' i <- BCIspeed$species=="ocelot"
#' sp <- log(BCIspeed$speed[i])
#' tm <- BCIspeed$time[i]*2*pi
#' mod <- fitlincirc(tm, sp, reps=50)
#' plot(mod, CircScale=24, xaxp=c(0,24,4), 
#'      xlab="Time", ylab="log(speed m/s)")
#' legend(8,-3, c("Fitted speed", "Null CI"), col=1:2, lty=1:2)
#' @export
fitlincirc <- function(circdat,lindat,pCI=0.95,reps=1000,res=512)
{ if(length(lindat)!=length(circdat)) 
  stop("lindat and circdat lengths are unequal")
  if(min(circdat)<0 | max(circdat)>2*pi) 
    stop("circdat values not between 0 and 2*pi, expecting radian data")
  if(max(circdat)<1) 
    warning("max(circdat) < 1, expecting radian data")
  
  n <- length(circdat)
  x <- seq(0,2*pi,2*pi/res)
  
  bs <- pbapply::pbsapply(1:reps, function(i)
  {  j <- sample(1:n,n,TRUE)
     lincircKern(x,circdat[j],lindat)
  })
  nulllcl <- apply(bs,1,quantile,(1-pCI)/2)
  nullucl <- apply(bs,1,quantile,(1+pCI)/2)
  fit <- lincircKern(x,circdat,lindat)
  p <- sapply(1:(res+1), function(i)
  { f <- ecdf(bs[i,])
    f(fit[i])
  })
  p[p>0.5] <- 1-p[p>0.5]
  p <- 2*p
  
  new("lincircmod", data=data.frame(circdat=circdat, lindat=lindat),
      fit=data.frame(x=x, fit=fit, p=p, nullLCL=nulllcl, nullUCL=nullucl))
}

#' Plot linear-circular relationship
#' 
#' Plot linear against circular data along with the fitted and null confidence limit distributions from a fitted \code{lincircmod} object.
#' 
#' @param x Object of class \code{lincircmod}.
#' @param CircScale Single numeric value defining the plotting maximum of the circular scale.
#' @param tlim Numeric vector with two elements >=0 and <=1 defining the lower and upper limits at which to plot distributions; default plots the full range.
#' @param fcol,flty,ncol,nlty Define line colour (\code{col}) and type (\code{lty}) for fitted (\code{f}) and null (\code{n}) distributions; input types as for \code{col} and \code{lty}, see \code{\link{par}}.
#' @param ... Additional arguments passed to the inital plot construction, affecting axes and data plot symbols.
#' @docType methods
#' @rdname plot-lincircmod-methods
#' @export
setMethod("plot", "lincircmod",
          function(x, CircScale=2*pi, tlim=c(0,1), fcol="black", flty=1, ncol="red", nlty=2, ...)
          { if(min(tlim)<0 | max(tlim)>1 | length(tlim)!=2) stop("tlim should contain two values >=0 and <=1")

            x<-mod
            fit <- x@fit
            dat <- x@data
            xx <- fit$x*CircScale/(2*pi)
            LinearData <- dat$lindat
            CircularData <- dat$circdat*CircScale/(2*pi)
            range <- tlim*CircScale
            plot(CircularData, LinearData, ...)
            if(range[1]<range[2])
            { i <- xx>=range[1] & xx<=range[2]
              lines(xx[i], fit$fit[i], col=fcol, lty=flty)
              lines(xx[i], fit$nullLCL[i], col=ncol, lty=nlty)
              lines(xx[i], fit$nullUCL[i], col=ncol, lty=nlty)
            } else 
            {  i <- xx>=range[1]
               lines(xx[i], fit$fit[i], col=fcol, lty=flty)
               lines(xx[i], fit$nullLCL[i], col=ncol, lty=nlty)
               lines(xx[i], fit$nullUCL[i], col=ncol, lty=nlty)
               i <- xx<=range[2]
               lines(xx[i], fit$fit[i], col=fcol, lty=flty)
               lines(xx[i], fit$nullLCL[i], col=ncol, lty=nlty)
               lines(xx[i], fit$nullUCL[i], col=ncol, lty=nlty)
            }
          }
)
