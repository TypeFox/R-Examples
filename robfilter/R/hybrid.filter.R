# hybrid.filter - robust time series filters based on robust location and /or regression estimators                        #
#                                                                                                                          #
# Authors:                                                                                                                 #
#   Prof. Dr. Roland Fried         <fried@statistik.uni-dortmund.de>                                                       #
#   Dipl.-Stat. Karen Schettlinger <schettlinger@statistik.uni-dortmund.de>                                                #


hybrid.filter <- function(y,                # input time series data (numeric vector or ts object)               
                          width,            # window width                                                       
                          method="all",     # filtering methods                                                  
                          minNonNAs=3,      # required minimal number of non-missing values in each window (half)
                          extrapolate=TRUE  # indicator for extrapolation of the estimations to the edges        
                          ){
### FIXME:
sub.minNonNAs <- minNonNAs # minimal number required in each window half -> ALS ZUSÄTZLICHE OPTION???
###

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

  # Errors concerning the filter methods
  all.methods <- c("MED", "RM", "MEAN", 
                   "FMH", "PFMH", "CFMH",
                   "MH", "PRMH", "CRMH", 
                   "MMH", "PRMMH", "CRMMH")       # character string of all possible methods 
  RMH.methods <- c("PRMH","CRMH","PRMMH","CRMMH") # methods using half window RM subfilters
  # character string of chosen methods (method.names):
  if (any(method == "all")){
    method.names <- all.methods
  } else {
    method.names <- method
  }
  n.met  <- length(method.names)
  if (!all(method.names %in% all.methods)){
    stop("invalid specification of 'method': possible values are ", paste(all.methods, collapse=", "))
  }

  # Errors concerning the window widths and minimum number of non-missing values
  if (missing(width)){
    stop("argument 'width' is missing with no default")
  }
  if ( (identical(all.equal(width%%2,0),TRUE)) | (!identical(all.equal(width%%1,0),TRUE)) | (width < 3) ){
    stop("argument 'width' must be an odd positive integer >= 3")
  }
  if ( (!identical(all.equal(minNonNAs%%1,0),TRUE)) | (minNonNAs < 0) ){
    stop("argument 'minNonNAs' must be an integer >= 0")
  }
  if (width < minNonNAs)  { 
    stop(paste("argument 'width' cannot be smaller than 'minNonNAs' (=",minNonNAs,")",sep=""))
  }
  if (length(y) < width){
    stop("argument 'y' must contain at least 'width' elements")
  }
  
  ## Internal indices
  N <- length(y)
  m <- floor(width/2)
  i <- -m:m
  i.F <- (-m):(-1)
  i.B <- 1:m
  ## weights for weighted average = LS fit
#  h <- (4*m - 6*c(1:m) + 2) / (m*(m-1))                                                    #### RAUSNEHMEN??? JETZT: lm()

  ## Information message:
  if( any(is.na(y)) ){ 
   cat(paste( sum(is.na(y)), "out of", N, "time series values in", ts.name, "are missing. \n"))
  }
    
  ## Preallocate space for level estimations (filled with NAs)
    results <- vector("list",length(method.names))
    for(mn in seq(along=method.names)){ 
        results[[mn]] <- rep(NA,N)
    }
    names(results) <- method.names

  ## Preallocate space for slopes
  RM.b   <- rep(NA,N)
  RM.F.b <- rep(NA,N)
  RM.B.b <- rep(NA,N)
  LS.F.b <- rep(NA,N)
  LS.B.b <- rep(NA,N)

# # # # #
# main programme loop:
for (t in (m+1):(N-m)) {
# # # # #

  ## Data in the current window around time 't'
  yw <- y[t+i]
  ## Values left of the window centre
  ywF <- yw[1:m]
  ## Values right of the window centre
  ywB <- yw[(m+2):width]    
  ## y.C:   Center value
  y.C <- y[t]

  ## 'Enough' observations in the whole window?
  enough.yw  <- (sum(!is.na(yw)) >= minNonNAs) & (sum(!is.na(yw)) > 0)
  ## 'Enough' observations in the left (F) and right (B) window half?
  enough.ywF <- (sum(!is.na(ywF)) >= sub.minNonNAs) & (sum(!is.na(ywF)) > 0)
  enough.ywB <- (sum(!is.na(ywB)) >= sub.minNonNAs) & (sum(!is.na(ywB)) > 0)
    
  ########################################################################
  ## Non-robust (sub-)filters
  ########################################################################
  if( any(method.names %in% c("MEAN","FMH","PFMH","CFMH")) ) {
# # # ## Running mean
    if( ("MEAN" %in% method.names) & enough.yw ) {
      results$MEAN[t] <- mean(yw, na.rm=TRUE)
    }

    if( any(method.names %in% c("FMH","PFMH","CFMH")) ) {
      ## Half-window means (for 'enough' obs. in each window half)
      if( enough.ywF ){
          mean.F <- mean(ywF, na.rm=TRUE)
      } else {
          mean.F <- NA
      }
      if( enough.ywB ) {
          mean.B <- mean(ywB, na.rm=TRUE)
      } else {
          mean.B <- NA
      }

# # # ## FIR median hybrid filter
      if( "FMH" %in% method.names ) {
        cFMH <- c(mean.F, y.C, mean.B)
        if( !all(is.na(cFMH)) ){
          results$FMH[t] <- median( cFMH, na.rm=TRUE)
        }
      }

      if( any(method.names %in% c("PFMH","CFMH")) ) {
      ## Half-window weighted means / LS prediction for centre
        if( enough.ywF ){
#          wm.F <- sum(na.omit(rev(h)*ywF))
          lmF       <- lm(ywF ~ i.F)
          wm.F      <- lmF$coef[1]
          LS.F.b[t] <- lmF$coef[2]
        } else {
          wm.F <- NA
        }
        if( enough.ywB ) {
#          wm.B <- sum(na.omit(h*ywB))
          lmB       <- lm(ywB ~ i.B)
          wm.B      <- lmB$coef[1]
          LS.B.b[t] <- lmB$coef[2]
        } else {
          wm.B <- NA
        }

# # #   ## Predictive FIR median hybrid filter
        if( "PFMH" %in% method.names ) {
          cPFMH <- c(wm.F,y.C,wm.B)
          if( !all(is.na(cPFMH)) ){
            results$PFMH[t] <- median( cPFMH, na.rm=TRUE)
          }
        }

# # #   ## Combined FIR median hybrid filter
        if( "CFMH" %in% method.names ) {
          cCFMH <- c(wm.F,mean.F,y.C,mean.B,wm.B)
          if( !all(is.na(cCFMH)) ){
            results$CFMH[t] <- median( cCFMH, na.rm=TRUE)
          }
        }

      } # end of PFMH, CFMH
    } # end of FMH, PFMH, CFMH
  } # end of non-robust filters
  ########################################################################


  ########################################################################
  ## Robust (sub-)filters
  ########################################################################
# # ## Simple RM filter 
  if( ("RM" %in% method.names) & (sum(!is.na(yw)) >= minNonNAs) ) {
    ## RM: Centre intercept estimation by RM regression
    index  <- seq(along=i)
    slopes <- sapply(index, function(k) median((yw[k] - yw[-k])/(index[k]-index[-k]), na.rm=TRUE))
    if( all(is.na(slopes)) ) {
      RM.slope <- NA
    } else {
      RM.slope <- median(slopes, na.rm=TRUE)
    }
    RM.b[t] <- RM.slope
    results$RM[t] <- median(yw - RM.slope*i, na.rm=TRUE)
  }
 
  ## Filters with whole-window median subfilter
  if( any(method.names %in% c("MED","MMH","PRMMH","CRMMH")) ){
    ## Whole-window median  
    if( enough.yw ){
      M.C <- median(yw, na.rm=TRUE)
    } else {
      M.C <- NA
    }
# # ## Running median 
    if ("MED" %in% method.names ){
      results$MED[t] <- M.C
    }
  }
 
  # # # # # # # # # # # # # # # # #
  ## Half window median subfilters 
  if( any(method.names %in% c("MH","CRMH","MMH","CRMMH")) ){
    ## Forward prediction by the median (if 'enough' values in left window half)
    if( enough.ywF ){
      M.F <- median(ywF,na.rm=TRUE)
    } else {
      M.F <- NA
    }
 
    ## Backward prediction by the median (if 'enough' values in right window half)
    if( enough.ywB ) {
      M.B <- median(ywB, na.rm=TRUE) 
    } else {
     M.B <- NA
   }

# # ## Median hybrid filter
   if( "MH" %in% method.names ) {
     cMH <- c(M.F, y.C, M.B)
     if( !all(is.na(cMH))){
       results$MH[t] <- median( cMH, na.rm=TRUE )
     }
   }

# # ## Median/median hybrid filter
   if( "MMH" %in% method.names ) {
     cMMH <- c(M.F, M.C, M.B)
     if( !all(is.na(cMMH)) ){
       results$MMH[t] <- median( cMMH, na.rm=TRUE )
     }
   }

 } # end of median subfilters     
 # # # # # # # # # # # # # # # # #


 # # # # # # # # # # # # # # # # #
 ## Half window RM subfilters     
 if( any(method.names %in% RMH.methods) ){
   ## Forward prediction by RM regression (if 'enough' values in left window half)
   if( sum(!is.na(ywF)) >= sub.minNonNAs ){
     index <- seq(along=ywF)
     slopes <- sapply(index, function(k) median((ywF[k] - ywF[-k])/(index[k]-index[-k]), na.rm=TRUE))
     if(all(is.na(slopes))){
       RM.F.slope <- NA
     } else{
       RM.F.slope <- median(slopes, na.rm=TRUE)
     }
     RM.F.b[t] <- RM.F.slope
     RM.F <- median(ywF - RM.F.slope*i.F, na.rm=TRUE)
   } else {
     RM.F <- NA
   }

   ## Backward prediction by RM regression (if 'enough' values in right window half)
   if( sum(!is.na(ywB)) >= sub.minNonNAs ) {
     index  <- seq(along=ywB)
     slopes <- sapply(index, function(k) median((ywB[k]-ywB[-k])/(index[k]-index[-k]),na.rm=TRUE))
     if (all(is.na(slopes))){
         RM.B.slope <- NA
     } else {
         RM.B.slope <- median(slopes, na.rm=TRUE)
     }
     RM.B.b[t] <- RM.B.slope
     RM.B <- median( ywB - RM.B.slope*i.B, na.rm=TRUE)
   } else {
     RM.B <- NA
   }
   
# # # ## Predictive RM hybrid filter
    if( "PRMH" %in% method.names ) {
      cPRMH <- c(RM.F, y.C, RM.B)
      if( !all(is.na(cPRMH)) ){
        results$PRMH[t] <- median( cPRMH, na.rm=TRUE )
      }
    }

# # # ## Predictive RM/median hybrid filter
    if( "PRMMH" %in% method.names ) {
      cPRMMH <- c(RM.F, M.C, RM.B)
      if( !all(is.na(cPRMMH)) ){
        results$PRMMH[t] <- median( cPRMMH, na.rm=TRUE )
      }
    }

# # # ## Combined RM hybrid filter
    if( "CRMH" %in% method.names ) {
      cCRMH <- c(RM.F, M.F, y.C, M.B, RM.B)
      if( !all(is.na(cCRMH)) ){
        results$CRMH[t] <- median( cCRMH, na.rm=TRUE )
      }
    }

# # # ## Combined RM/median hybrid filter
   if( "CRMMH" %in% method.names ) {
     cCRMMH <- c(RM.F, M.F, M.C, M.B, RM.B)
     if( !all(is.na(cCRMMH)) ){
       results$CRMMH[t] <- median( cCRMMH, na.rm=TRUE)
     }
   }

 } # end of RM subfilters         
 # # # # # # # # # # # # # # # # #

# # # # #
} # end of for-loop over all time points
# # # # #


  # Storing all (half window) LS and RM slopes in one data frame 
  # and extrapolating the level estimations to the edges of the time series based on the slope estimations (for regression-based subfilters)
  slope   <- NULL
  s.names <- NULL
  
  # Slopes for extrapolation of the trend in the first window to the left time series edge / in the last window to the right edge of time
  bl <- rep(NA,n.met)
  br <- rep(NA,n.met)

  if( "RM" %in% method.names ) { 
    slope   <- cbind(slope,RM.b) 
    s.names <- c(s.names,"RM")
    if(extrapolate){
      w.RM  <- which(method.names == "RM")
      bl[w.RM]  <- RM.b[m+1]
      br[w.RM]  <- RM.b[N-m]
    }
  }
 
  if( any(method.names %in% RMH.methods) ) {
    slope   <- cbind(slope, RM.F.b, RM.B.b) 
    s.names <- c(s.names, "RM.left", "RM.right")
    if(extrapolate){
      w.RMH <- which(method.names %in% RMH.methods)
      bl[w.RMH] <- RM.F.b[m+1] 
      br[w.RMH] <- RM.B.b[N-m]
    }
  } 

  if( any(method.names %in% c("PFMH","CFMH")) ) { 
    slope   <- cbind(slope, LS.F.b, LS.B.b) 
    s.names <- c(s.names, "LS.left", "LS.right")
    if(extrapolate){
      w.LSH <- which(method.names %in% c("PFMH","CFMH"))
      bl[w.LSH] <- LS.F.b[m+1] 
      br[w.LSH] <- LS.B.b[N-m]
    }
  }

  if( any(method.names %in% c("MED","MEAN","FMH","MH","MMH")) & extrapolate) { 
    w.loc <- which(method.names %in% c("MED","MEAN","FMH","MH","MMH"))
    bl[w.loc] <- 0
    br[w.loc] <- 0
  } 
  
  dimnames(slope)[[2]] <- s.names
  if( extrapolate ) {
    # Extrapolation for the slope values
    if(!is.null(slope)){
      slope <- apply(slope,2,function(x) {x[1:m] <- x[m+1]; x[(m+1):(N-m)] <- x[(m+1):(N-m)]; x[(N-m+1):N] <- x[N-m]; return(x)})
    }
    # Extrapolation for the level
    for(i in 1:n.met){
        results[[i]][1:m]       <- results[[i]][m+1] - (m:1)*bl[i]
        results[[i]][(N-m+1):N] <- results[[i]][N-m] + (1:m)*br[i]
    }
  }

  return( structure( list( level=as.data.frame(results), slope=as.data.frame(slope),
                           y=y, width=width, method=method, extrapolate=extrapolate, ts.name=ts.name),
          class="hybrid.filter")
        )
}



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Default output
print.hybrid.filter <- function(x, ...) {

    # Number of filter methods
    n.met   <- dim(x$level)[2]

    S <- summary(x)
    N <- length(x$y)
    start <- (x$width + 1)/2
    if(all(is.na(x$level[start+c(0:2),]))){ 
        start <- min(which(apply(x$level, 1, function(L){any(!is.na(L))}) )) 
    }

    # Print level output
    level.output <- rbind(rep("...",n.met),round(x$level[start+c(0:2),],6),rep("...",n.met))
    rownames(level.output)[c(1,5)] <- c(" ","  ")
    cat("$level \n") # print level
    print(level.output)
    cat("('",S["level","Class"],"' with ",N," obs. of", S["level","Length"], " variables)\n \n",sep="")
    # Print slope output if not NULL
    num.s <- as.numeric(S["slope","Length"])
    if( num.s > 0 ) {
      if(num.s ==1) {
        slope.output <- as.data.frame(c("...",round(x$slope[start+c(0:2),],6),"..."))
        dimnames(slope.output)[[2]] <- names(x$slope)
        dimnames(slope.output)[[1]] <- c(" ",start+c(0:2),"  ")
        vars <- "variable"
      } else {
        slope.output <- rbind(rep("...",num.s),round(x$slope[start+c(0:2),],6),rep("...",num.s))
        rownames(slope.output)[c(1,5)] <- c(" ","  ")
        vars <- "variables"
      }
        cat("$slope \n")
        print(slope.output)
        cat("('",S["slope","Class"],"' with ",N," obs. of ", num.s, " ", vars,")\n \n",sep="")
    }
}



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Default plot
plot.hybrid.filter <- function(x, ...) {
    # Length of the time series
    N <- length(x$y)
    
    # Names & number of filter methods
    m.names <- names(x$level)
    n.met   <- length(m.names)
    
    # Setting the y-limits
    ylims <- c(min(x$y,min(x$level,na.rm=TRUE),na.rm=TRUE),max(x$y,max(x$level,na.rm=TRUE),na.rm=TRUE))
    xlims <- c(1,N)
    
    # Defining the title
    if(n.met > 1) { 
        titel <- "Hybrid Filters"
    } else {
        titel <- "Hybrid Filter"
    }

    # Possible colors
    mcols <- c("red","green3","blue","skyblue2","orange","yellow","grey","purple2","darkgreen","cyan","brown4","pink")

    # Plot
    par(mar=c(4,4,4,7),oma=rep(0,4),mgp=c(2.5,1,0))
    plot(x$y,xlim=xlims,ylim=ylims,type="l",xlab="Time", ylab=x$ts.name, main=titel)
    
    for(i in seq(along=m.names)){
        lines(x$level[[i]],col=mcols[i],lwd=2)
    }

    # Legend
    par(xpd=TRUE,cex=0.8)
#    legend(par("usr")[2],mean(c(par("usr")[3],par("usr")[4])),c(x$ts.name,m.names),xjust=0,yjust=0.5,lty=rep(1,n.met+1),lwd=c(1,rep(2,n.met)),col=c("black",mcols[1:n.met]),bty="n")
    legend(par("usr")[2],mean(c(par("usr")[3],par("usr")[4])),c("Time Series",m.names),xjust=0,yjust=0.5,lty=rep(1,n.met+1),lwd=c(1,rep(2,n.met)),col=c("black",mcols[1:n.met]),bty="n")
    par(xpd=FALSE,cex=1)
}
