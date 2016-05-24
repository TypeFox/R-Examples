# robust.filter - robust time series filters based on robust regression techniques,                                        #
#                 with outlier replacement (based on robust scale estimation), shift detection and window width adaptation #
#                                                                                                                          #
# Authors:                                                                                                                 #
#   Prof. Dr. Roland Fried         <fried@statistik.uni-dortmund.de>                                                       #
#   Dipl.-Stat. Karen Schettlinger < schettlinger@statistik.uni-dortmund.de>                                               #

# Bugs:                                                                                                                    #
# Missing values in the time series cannot be handled, yet.                                                                #

                                                
robust.filter <- function(y,                     # input time series data (numeric vector or ts object)
                          width,                 # (minimal) window width
                          trend="RM",            # (regression) method for the level estimation
                          scale="QN",            # method for robust scale estimation
                          outlier="T",           # method for outlier treatment
                          shiftd=2,              # factor regulating the size of shifts to be detected
                          wshift=floor(width/2), # number of observations for shift detection
                          lbound=0.1,            # lower bound for scale estimation
                          p=0.9,                 # fraction of identical values for special rules
                          adapt=0,               # fraction regulating the adaption of the window width
                          max.width=width,       # maximal window width
                          online=FALSE,          # indicator for estimation at the rightmost point or in the centre of each time window
                          extrapolate=TRUE       # indicator for extrapolation of the estimations to the edges
                          ) {

  ## save the name of the input time series
  ts.name <- deparse(substitute(y))

  ## Stopping rules and error messages

  # Errors concerning the input time series
  if(is.ts(y)){y <- as.numeric(y)} # coerce a time series object into a numeric vector
  if (!is.null(dim(y))){
    stop("\n        argument 'y' must be a one dimensional numeric vector or time series object")
  }
  if ( any(is.infinite(y), na.rm=TRUE) ){
    stop("argument 'y' contains Inf or -Inf values")
  }
  if ( any(is.na(y)) ){
    stop("argument 'y' contains missing values")
  }
  
  ## Error for missing width
  if (missing(width)){
    stop("argument 'width' is missing with no default")
  }

  # length of the time series
  N <- length(y)
  # minimal window width (for an adaptive filter)
  min.width <- width

  # Errors concerning the methods for level extraction, scale estimation and outlier treatment
  all.trend   <- c("RM", "MED", "LMS", "LTS")
  all.scale   <- c("MAD", "QN", "SN", "LSH")
  all.outlier <- c("T", "L", "M", "W")
  
  if( !(trend %in% all.trend) ) {
    stop("invalid specification of 'trend': possible values are ", paste(all.trend, collapse=", "))
  }
  if( !(scale %in% all.scale) ) {
    stop("invalid specification of 'scale': possible values are ", paste(all.scale, collapse=", "))
  }
  if( !(outlier %in% all.outlier) ) {
    stop("invalid specification of 'outlier': possible values are ", paste(all.outlier, collapse=", "))
  }
  
  ## Error concerning wshift
  if( wshift==0 ) { stop("invalid specification of 'wshift': value must be positive") } 
  
  ## Errors concerning the window width &
  ## Internal indices, window width and functions for online level estimation
  if( !online ){
    if( (!identical(all.equal(width%%1,0),TRUE)) | (identical(all.equal(width%%2,0),TRUE)) | (width < 3) ) {
      stop("argument 'width' must be an odd positive integer >= 3 for offline estimation ('online=FALSE')")
    }
    m <- (width-1)/2  # Window half-width
    h<-m              # Delay
    ho<-0
    if (adapt>0){
      madapt <- ceiling(m/2)
    }
    window.width <- 2 * m + 1
    ## Helper vector of indices for a symmetric window around some point 't' with
    ## width 'window.width'
    j <- (-m):m
    b <- m
    t <- m + 1
    ## Assign function to estimate level:
    trend.extraction <- match.arg(trend, all.trend)
    reg <- match.trend.extraction(trend.extraction)

  } else {
  ## Internal indices, window width and functions for offline level estimation
    if( (!identical(all.equal(width%%1,0),TRUE)) | (width < 3) ) {
      stop("argument 'width' must be a positive integer >= 3")
    }
    m <- width-1
    h <-0
    ho <-0
    window.width <- width
    if( adapt>0 ){
      madapt=ceiling(m/3)
    }
    j <- (-width+1):0
    b <- width-1
    t <- width
#    if (trend !="RM" && trend != "MED"){
#        stop("Only RM and MED implemented for full online modus.")
#        }
    reg <- switch(trend,
                  MED=trendMED,
                  RM=RMon,
                  LMS=LMSon,
                  LTS=LTSon)
    trend.extraction <- match.arg(trend, all.trend)
  }

  ## Further errors concerning the options
  if (N < width){
    stop("argument 'y' must contain at least 'width' elements")
  }
  # shiftd
#  if( (!identical(all.equal(shiftd%%1,0),TRUE)) | (shiftd < 0) ) {
#    stop("argument 'shiftd' must be a positive integer")
#  }
  if( shiftd < 0 ) {
    stop("argument 'shiftd' must be a positive numeric value")
  }
  # wshift
  if( (!identical(all.equal(wshift%%1,0),TRUE)) | (wshift < 0) ) {
    stop("argument 'wshift' must be a positive integer")
  }
  mh <- floor(width/2)  # Window half-width for checking for too many replacements
  if (mh<wshift) {
    stop("argument 'wshift' cannot be larger than half the 'width'")
  }
  # max.width
  if( (!identical(all.equal(max.width%%1,0),TRUE)) | (max.width < width) ) {
    stop("argument 'max.width' must be a positive integer >= 'width'")
  }
  # lbound
  if( lbound < 0){
    stop("argument 'lbound' must be a positive real value")
  }

  # adapt
  if( (adapt != 0) & ((adapt < 0.6)|(adapt > 1)) ){
    stop("argument 'adapt' must be either 0 or a numeric value in [0.6,1]")
  }
#  if (adapt>1 || adapt<0 ){
#    stop("adapt must be between 0.6 and 1, or 0")
#  }
#  if (0<adapt && adapt<0.6 ){
#    stop("adapt must be between 0.6 and 1, or 0")
#  }

  # p
  # if p is not in [2/3, 1] you cannot say that 100p percent dominate the window
  if ( (p < 2/3) | (p > 1) ) {
    stop("argument 'p' must be a real value in [2/3,1].")
  }

  ## Retrive function for scale estimation.
  scale.estimation <- match.arg(scale, c("QN","SN","MAD","LSH"))
  scale <- match.scale.estimator(scale.estimation)

  ## scale.est - helper function to estimate scale.
  scale.est <- function(x, corr, fallback = NULL) {
      sigma <- fallback
      if (d1 == 0) {
          if (sum(x != 0) > 1) {
              sigma <- scale(res.t[res.t != 0], corr)
          }
      }
      else {
          sigma <- scale(x, corr)
      }
      ## Make sure we don't go below lbound.
      if (is.na(sigma) || sigma < lbound) {
          sigma <- lbound
      }
      return(sigma)
  }

  ## Factors for time correction in the scale estimations
  timecor <- timecorrection[, paste(outlier, scale.estimation, sep = "-")]
  if (N > length(timecor)){
      timecor <- c(timecor, rep(timecor[length(timecor)], N - length(timecor)))
  }

  ## Set constants / factors for outlier treatment
  if (outlier == "T") { # Trimming
      d0 <- 3
      d1 <- 0
  } else if (outlier == "L") { # downweighting large outliers
      d0 <- 3
      d1 <- 1
  } else if (outlier == "M") { # downweighting moderately sized outliers
      d0 <- 2
      d1 <- 1
  } else if (outlier == "W") { # winsorising
      d0 <- 2
      d1 <- 2
  }
  outlier.treatment <- outlier
  
## Preallocate space for the results
  ## Vector of observations with outliers replaced according to 'outlier.treatment'
  y.tilde <- y
  ## Outlier vector with backwards substituted values after level shift
  ol <- numeric(N)
  ## Outlier vector showing positive/negative outliers
  outlier <- integer(N)
  mu <- rep(NA, N)            # signal level
  beta <- rep(NA, N)          # slope
  sigma <- rep(NA, N)         # variability / scale
  level.shift <- rep(NA, N)   # vector indicating the time points of the shifts (and their detection)

  ## internal initial values
  initialise <- 1
  level.shift[1] <- 0
  ma <- m


########################################################################
  ## Starting the main programm loop #######
    while (t <= N - h) {
  # indicator for adaption of window width
        adind=0
        if (adapt>0 && online){
          j=(-ma:0)
          width=ma+1
        }
        if (adapt>0 && !online){
          j=(-ma:ma)
          wshift=ma
          width=2*ma+1
          madapt=ceiling(ma/2)
        }
        mh <- floor(width/2)
        j1 <- (-999)
        level.shift [t] <- 0
####### Step Ia: I N I T I A L  F I T T I N G  ###############################
        if (initialise == 1) {
            tmp <- reg(y.tilde[t + j])
            mu[t] <- tmp[1]
            beta[t] <- tmp[2]
            res.t <- y.tilde[t + j] - (mu[t] + j * beta[t]) # residual vector for this window
            n.replaced <- sum(y[t + j] != y.tilde[t + j])   # number of replaced observations
            sigma[t] <- scale.est(res.t, width - n.replaced, sigma[t - 1])# * timecor[t - ho]
      # Step Ib: Replacing and indicating outliers
            first.res <- y.tilde[t + j] - (mu[t] + j * beta[t])
            for (i in t + j) {
                if (abs(first.res[i - t + ma + 1]) > d0 * sigma[t]) {
                  y.tilde[i] <- mu[t] + beta[t] * (i - t) + d1 *
                    sign(first.res[i - t + ma + 1]) * sigma[t]
                    if ( first.res[i - t + ma + 1] > 0) {
                    outlier[i] <- 1
                    ol[i] <- 1
                  } else {
                    outlier[i] <- (-1)
                    ol[i] <- (-1)
                  }
                }
            }
        }
######## E N D  I N I T I A L I Z A T I O N ########

######## STEP 1: H A V E  T O O  M A N Y  O B S E R V A T I O N S  B E E N  R E P L A C E D? ###########
        if (sum(ol[t + j] > 0) > mh)           # positive outliers
            y.tilde[t + j][ol[t + j] > 0] <- y[t + j][ol[t + j] > 0]
        if (sum(ol[t + j] < 0) > mh)           # negative outliers
            y.tilde[t + j][ol[t + j] < 0] <- y[t + j][ol[t + j] < 0]
        if (sum(ol[t + j] == 0) < floor(mh/3)) # all outliers
            y.tilde[t + j] <- y[t + j]
######## STEP 2: L O C A L  M O D E L  F I T T I N G ###########################################
        tmp <- reg(y.tilde[t + j])
        mu[t] <- tmp[1]
        beta[t] <- tmp[2]

    ## Estimating location in the case of only two or three dominating
    ## values in at least 100p percent of the window
    y.freq <- as.matrix(table(y.tilde[t+j]))
    ## index for the values with a frequency higher than p/3
    i3     <- which(y.freq >= 1/3*p*width)
    ## index for the values with a frequency higher than p/2
    i2     <- which(y.freq >= 1/2*p*width)

    ## The exceptional rules are applied if there are two or three
    ## values which together amount to at least 100p percent of the
    ## window
    if (nrow( y.freq) == 2
        || nrow(y.freq) == 3
        || (any(y.freq >= 1/3*width)
            && length(i3) >= 2)) {
      if (length(i2) == 2 || nrow(y.freq) == 2 || length(i3) == 2) {
        if (length(i2) == 2 || nrow( y.freq) == 2) {
          ## mean as location estimate for two differing values of
          ## y.tilde
          midlevel <- mean(as.numeric(dimnames(y.freq)[[1]][i2]))
        } else {
          ## mean as location estimate for two differing values of
          ## y.tilde
          midlevel <- mean(as.numeric(dimnames(y.freq)[[1]][i3]))
        }
      } else {
        ## medium level as location estimate for three differing
        ## values of y.tilde
        midlevel <- as.numeric(dimnames(y.freq)[[1]][2])
      }
      ## assigning the midlevel as the estimated location and 0 as the
      ## slope at time t
      mu[t] <- midlevel
      beta[t] <- 0
    }

        res.t <- y.tilde[t + j] - (mu[t] + j * beta[t]) # residual vector for this window
        n.replaced <- sum(y[t + j] != y.tilde[t + j])   # number of replacements
        sigma[t] <- scale.est(res.t, width - n.replaced, sigma[t - 1])   #*timecor[t - m]
######### E X T R A P O L A T I O N  A T  I N I T I A L I Z A T I O N  #####################
        if (initialise == 1) {
            tmp <- -b:-1
            mu[t + tmp] <- mu[t] + tmp * beta[t]
            beta[t + tmp] <- beta[t]
            sigma[t + tmp] <- sigma[t]
            level.shift[(t - b+1):(t - 1)] <- 0
        }

########### STEP 3: L E V E L  S H I F T  D E T E C T I O N ##################
    ## Write out the right part of the residual vector at point t
    ## based on the 'true' data vector without replacements for outliers
        res.y.t <- y[t + j] - (mu[t] + j * beta[t])
         right.res.t <- res.y.t[(width-wshift+1):width]

        if (sum(right.res.t > shiftd * sigma[t]) > floor(wshift/2) +   1) {
      ## Positive level shift:
      ## Calculate the index 'j1' for the time point 't + j1' where the level shifted
      ## and set the indicator variable to the time point of detection (for a positive shift)
            index <- 1:wshift
            j1 <- min(index[right.res.t > shiftd * sigma[t]])
            if( online ){
              ls <- t -wshift+ j1
            } else {
              ls <- t +(ma-wshift)+ j1
              level.shift[(t + 1):(ls - 1)] <- 0
            }
            level.shift[ls] <- t
        }  else if (sum(right.res.t < (-shiftd) * sigma[t]) > floor(wshift/2) +   1) {
       ## Negative level shift:
      ## Calculate the index 'j1' for the time point 't + j1' where the level shifted
      ## and set the indicator variable to -1 (for a negative shift)
            index <- 1:wshift
            j1 <- min(index[right.res.t < (-shiftd) * sigma[t]])
            if( online ) {
              ls <- t -wshift+ j1
            } else {
              ls <- t +(ma-wshift)+ j1
              level.shift[(t + 1):(ls - 1)] <- 0
            }
            level.shift[ls] <- -t
        }

      ## STEP 4:  Rules in case of a level shift
        if (j1 != -999) {
         ## Extrapolating the old trend until before the level shift
            if( online ) {
              if (level.shift[ls] > 0) {
                 y.tilde [(ls):(t)][ol[ls:(t)] > 0] <- y[ls:(t)][ol[ls:(t)] > 0]
                 outlier[(ls):(t)][ol[ls:(t)] > 0] <- 0
              } else if (level.shift[ls] < 0) {
                 y.tilde[(ls):(t)][ol[ls:(t)] < 0] <- y[ls:(t)][ol[ls:(t)] < 0]
                 outlier[(ls):(t)][ol[ls:(t)] < 0] <- 0
              }
              level.shift[t+1] <- 0
              ma=m
              b<-t-ls
              t <- ls + ma-wshift +1
              if (extrapolate){
               b<-ma-wshift+1
              } else{
                b<-ma-wshift-b+1
              }
            } else {
              mu[(t + 1):(ls - 1)] <- mu[t] + 1:(ma-wshift+j1 - 1) * beta[t]
              beta[(t + 1):(ls - 1)] <- beta[t]
              sigma[(t + 1):(ls - 1)] <- sigma[t]
      ## Replacing positive (negative) outliers after the level shift
      ## by the original y-values in case of a positive (negative)
      ## level shift and adjust the outlier vector
              if (level.shift[ls] > 0) {
                 y.tilde[(ls):(t + ma)][ol[ls:(t + ma)] > 0] <- y[ls:(t +ma)][ol[ls:(t + ma)] > 0]
                 outlier[(ls):(t + ma)][ol[ls:(t + ma)] > 0] <- 0
              } else if (level.shift[ls] < 0) {
                 y.tilde[(ls):(t + ma)][ol[ls:(t + ma)] < 0] <- y[ls:(t +ma)][ol[ls:(t + ma)] < 0]
                 outlier[(ls):(t + ma)][ol[ls:(t + ma)] < 0] <- 0
              }
              ma=m             # minimize window width
              t <- ls+ma       # move window after the shift
              b <- ma # t-ls   # number of observations to be extrapolated to the left
            }

#           print(c(t,ls,j1,b))
            initialise <- 1
            ol[t + j] <- 0
            }
######### R U L E S  I N  C A S E  O F  N O  S H I F T ################################
    if (j1 == -999) {
          initialise <- 0
######### STEP 5: A D A P T A T I O N   N E C E S S A R Y ? ###########################
            level.shift[t]=0
            if (adapt>0 && width>min.width){
                if( online ) {
                   ad.res=res.y.t[(width-madapt+1):width]
                   if (sum(ad.res<0)>adapt*length(ad.res)){adind=1}
                   if (sum(ad.res>0)>adapt*length(ad.res)){adind=1}
                } else {
                   ad.res=res.y.t[(madapt+1):(width-madapt)]
                   if (sum( ad.res<0)>adapt*length(ad.res)){adind=1}
                   if (sum(ad.res>0)>adapt*length(ad.res)){adind=1}
                }
              }
            if (adind==1 ) {ma=ma-1}
######### STEP 6: O U T L I E R  D E T E C T I O N  (at next time point) ###########################
            if (adind==0){
               if ( online ) {tp=1
                 mu[t]=max(mu[t],min(y[(t-wshift+1):t]))
                 mu[t]=min(mu[t],max(y[(t-wshift+1):t]))
               } else {
                 tp=ma+1
               }
                 tt=t+tp
                 next.res <- y.tilde[tt] - (mu[t] + tp *  beta[t])
                 if ((abs(next.res) > d0 * sigma[t]) && (tt <= N)) {
                   y.tilde[tt] <- mu[t] + tp * beta[t] + d1 * sign(next.res) * sigma[t]
                   if (next.res > 0) {
                     outlier[tt] <- 1
                     ol[tt] <- 1
                   } else {
                     outlier[tt] <- (-1)
                     ol[tt] <- (-1)
                   }
                 }
#                print(c(t,width,max.width))
                 t=t+1
                 if (adapt>0 && width<max.width && ma<(t-1)){
                    ma=ma+1
                 #  if (!online) {t=t-1}
                 }
                 if (!online) {
                   ma=max(min(ma,N-t),m)
                 }
              }
           }
#### E N D  O F  M A I N  P R O G R A M  L O O P ###########
         }
      t = t-1

    if( !online ) {
      tmp <- 1:h
      mu[t + tmp] <- mu[t] + tmp * beta[t]
      beta[t + tmp] <- beta[t]
      sigma[t + tmp] <- sigma[t]
      level.shift[t+tmp] <- 0
    }

# Output:
    return( 
      structure( 
        list(
          level = mu, slope = beta, sigma =  sigma,
          ol = ol, level.shift = level.shift, 
          y = y, width = width, 
          trend = trend.extraction, scale = scale.estimation, outlier = outlier.treatment,
          shiftd = shiftd, wshift = wshift, lbound = lbound, p = p, 
          adapt = adapt, max.width = max.width,
          online=online, extrapolate=extrapolate, ts.name=ts.name
        ), 
      class = "robust.filter")
    )       
}






# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Default output
print.robust.filter <- function(x, ...) {
    # length of the input time series
    N <- length(x$y)
    if(N <= 50){
      print(list(level=round(x$level,6), slope=round(x$slope,6), sigma=round(x$sigma,6)))
    } else {
      cat("$level")
      L <- data.frame("[1]",t(round(x$level[1:5],6)),"...")
      dimnames(L)[[1]] <- c(" ")
      dimnames(L)[[2]] <- c(" ","  ","   ", "    ","     ","      ","       ")
      print(L[1,])
      cat("$slope")
      Sl <- data.frame("[1]",t(round(x$slope[1:5],6)),"...")
      dimnames(Sl)[[1]] <- c(" ")
      dimnames(Sl)[[2]] <- c(" ","  ","   ", "    ","     ","      ","       ")
      print(Sl[1,])
      cat("$sigma")
      S <- data.frame("[1]",t(round(x$sigma[1:5],6)),"...")
      dimnames(S)[[1]] <- c(" ")
      dimnames(S)[[2]] <- c(" ","  ","   ", "    ","     ","      ","       ")
      print(S[1,])
      cat(N-5," observations omitted \n")
    }
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Default plot
plot.robust.filter <- function(x, ...) {
    # Length of the time series
    N <- length(x$y)

    # Setting the y-limits
    ylims <- c(min(x$y,min(x$level,na.rm=TRUE),na.rm=TRUE),max(x$y,max(x$level,na.rm=TRUE),na.rm=TRUE))
    xlims <- c(1,N)

    # Defining the title
    if(x$online){
      ol.text <- "Online "
    } else {
      ol.text <- ""
    }
    if(x$adapt==0){
      adapt.text <- paste(" with Window Width ",x$width,sep="")
    } else {
      adapt.text <- " with Adaptive Window Width"
    }
    titel <- paste(ol.text,x$trend," ",x$outlier,"-",x$scale," Filter",adapt.text,sep="")

    # Plot
    plot(x$y, type="l", main=titel, ylim=ylims, xlim=xlims,xlab="Time",ylab=x$ts.name)
    lines(x$level, col="red", lty=1, lwd=2)
    # Legend
    legend(x="topright", bty="n",legend=c("Time Series", "Filtered Signal"),lty=c(1, 1),col=c("black", "red"),lwd=c(1, 2))
}





#####
## Helper functions for online regression / level estimation

RMon <- function(x)   {
   l <- 1:length(x)
   medquotient.t <- c()
   for (k in l){ medquotient.t <- c(medquotient.t, median((x[k] - x[l!=k])/(k -
   l[l!=k])))}
   beta.RM <- median(medquotient.t)
   mu.RM   <- median(x + beta.RM*(length(x)-l))
   c(mu.RM, beta.RM)}

LMSon <- function(y) {
  n <- length(y)
  x<- (-n+1):0
  mdl <- lqs(x, y, method="lms")
  return(coef(mdl))
}

LTSon <- function(y) {
  n <- length(y)
  x <- (-n+1):0
  mdl <- lqs(x, y, method="lts")
  return(coef(mdl))
}
