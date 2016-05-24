# dw.filter - robust time series filters based on robust regression techniques in nested time windows                      #
#             with outlier detection based on robust scale estimation                                                      #
#                                                                                                                          #
# Authors:                                                                                                                 #
#   Prof. Dr. Roland Fried         <fried@statistik.uni-dortmund.de>                                                       #
#   Dipl.-Stat. Karen Schettlinger <schettlinger@statistik.uni-dortmund.de>                                                #

### FIXME: Output of sigma?

dw.filter <- function(y,                        # input time series data (numeric vector or ts object)                                
                      outer.width, inner.width, # window widths                                                                       
                      method="all",             # filtering methods                                                                   
                      scale="MAD",              # scale estimation method                                                             
                      d=2,                      # trimming constant (for a d*sigma rule)                                              
                      minNonNAs=5,              # required minimal number of non-missing values in each window                        
                      online=FALSE,             # indicator for estimation at the rightmost point or in the centre of each time window
                      extrapolate=TRUE          # indicator for extrapolation of the estimations to the edges                         
                      ) {
                      
  ## save the name of the input time series
  ts.name <- deparse(substitute(y))

  ## Stopping rules and error messages

  # Errors concerning the input time series
  if(is.ts(y)){y <- as.numeric(y)} # coerce a time series object into a numeric vector
  if (!is.null(dim(y))){
    stop("\n        argument 'y' must be a one dimensional numeric vector or time series object")
  }
  if ( any(is.infinite(y), na.rm=T) ){
    stop("argument 'y' contains Inf or -Inf values")
  }

  # Errors concerning the filter methods
  all.methods <- c("MED", "RM", "MTM", "TRM", "MRM", "DWMTM", "DWTRM", "DWMRM", "DWRM")  # character string of all possible methods
  # character string of chosen methods (method.names)
  if (any(method == "all")){
    method.names <- all.methods
  } else {
    method.names <- method
  }
  n.met  <- length(method.names)
  if (!all(method.names %in% all.methods)){
    stop("\n        invalid specification of 'method': possible values are ", paste(all.methods, collapse=", "))
  }

  # Errors concerning the window widths
  if(missing(outer.width)){
    if(missing(inner.width) & any(substr(method.names,1,2)=="DW")) { 
      stop("arguments 'outer.width' and 'inner.width' are missing with no default") 
    } else {
      stop("\n        argument 'outer.width' is missing with no default")
    }
  } else {
    if( (!identical(all.equal(outer.width%%1,0),TRUE)) | (outer.width < 0) ){
        stop("argument 'outer.width' must be a positive integer") 
    }
    if(missing(inner.width)){ 
        if(any(substr(method.names,1,2)=="DW")){
            stop("argument 'inner.width' is missing with no default")
        } else {
            inner.width <- outer.width
        }
    } else {
        if( (!identical(all.equal(inner.width%%1,0),TRUE)) | (inner.width < 0) ){
            stop("argument 'inner.width' must be a positive integer") 
        }
    }
  }
 
  if( (!online) & (identical(all.equal(outer.width%%2,0),TRUE)) ){ 
    stop("argument 'outer.width' must be an odd integer for offline estimation ('online=FALSE')") 
  }
  if( (!online) & (identical(all.equal(inner.width%%2,0),TRUE)) ){ 
    stop("argument 'inner.width' must be an odd integer for offline estimation ('online=FALSE')") 
  }

  if (outer.width < inner.width){ 
    stop("argument 'outer.width' must be larger than or equal to 'inner.width'")
  }
  if (inner.width < 3)          { 
    stop("argument 'inner.width' must be at least 3")
  }
  if (inner.width < minNonNAs)  { 
    stop(paste("argument 'inner.width' cannot be smaller than 'minNonNAs' (=",minNonNAs,")",sep=""))
  }
  if (length(y) < outer.width)  { 
    stop("argument 'y' must contain at least 'outer.width' elements")
  }

  # Errors concerning the methods for scale estimation, outlier treatment and treatment of missing values
  if ( !(scale %in% c("QN","SN","MAD")) ){ 
    stop("invalid specification of 'scale': possible values are MAD, QN and SN")
  }
  if (d < 0){ 
    stop("argument 'd' must be a non-negative number") 
  }

  # Errors concerning the minimum number of observations within one window
  if(minNonNAs < 0) {
     stop("argument 'minNonNAs' must be a positive integer")
  } else {
    if( identical(all.equal(minNonNAs%%1,0),TRUE) ){
        minNonNAs <- floor(minNonNAs)   # defining it as 'true' integer
    } else {
        stop("argument 'minNonNAs' must be a positive integer")
    }
  }


  # Internal indices within the time windows
  N <- length(y)
  if (online) {
    outer.i <- c((-outer.width+1):0)
    inner.i <- c((-inner.width+1):0)
  } else {
    m <- floor(outer.width/2)
    outer.i <- c(-m:m)
    l <- floor(inner.width/2)
    inner.i <- c(-l:l)
  }
  
  ## Information messages
  if( any(is.na(y)) ){ 
   cat(paste( sum(is.na(y)), "out of", N, "time series values in", ts.name, "are missing. \n"))
  }
  if ("MED" %in% method.names){
   cat("For the MED filter the window width 'outer.width' =",paste(outer.width), "was used. \n")
  }
  if ("RM" %in% method.names){
   cat("For the RM filter the window width 'outer.width' =",paste(outer.width), "was used. \n")
  }
  if ("MTM" %in% method.names){
   cat("For the MTM filter the window width 'outer.width' =",paste(outer.width), "was used. \n")
  }
  if ("TRM" %in% method.names){
   cat("For the TRM filter the window width 'outer.width' =",paste(outer.width), "was used. \n")
  }
  if ("MRM" %in% method.names){
   cat("For the MRM filter the window width 'outer.width' =",paste(outer.width), "was used. \n")
  }
  
  ## Preallocate space for slopes and results (filled with NAs)
    results <- vector("list",length(method.names))
    slopes  <- vector("list",length(method.names))
    for(mn in seq(along=method.names)){ 
        results[[mn]] <- rep(NA,N)
        slopes[[mn]]  <- rep(NA,N)
    }
    names(results) <- method.names
    names(slopes)  <- method.names

  ## Preallocate space for scales (filled with NAs)
  outer.loc.sigma <- rep(NA,N)
  outer.reg.sigma <- rep(NA,N)
  inner.loc.sigma <- rep(NA,N)
  inner.reg.sigma <- rep(NA,N)

  if (online){ 
    tpoints <- outer.width:N 
  } else {  
    tpoints <- (m+1):(N-m) 
  }

  ## Definition of the estimate used for scale estimation
    if (any(method.names %in% c("MTM","TRM","MRM","DWMTM","DWTRM","DWMRM"))) {
      if (scale == "MAD"){
###        scale.est <- rf.MAD
        scale.est <- mad
      }
      if (scale == "SN"){
###        scale.est <- rf.SN
        scale.est <- Sn
      }
      if (scale == "QN"){
###        scale.est <- rf.QN
        scale.est <- Qn
      }  
    }

# # # # #
# main programme loop:
  for (t in tpoints) {
# # # # #

  ## Observations in the current outer and inner window
  outer.y <- y[t+outer.i]
  inner.y <- y[t+inner.i]

  ## The following estimations only take place if there are a sufficient number of observations within the outer window
  if(sum(!is.na(outer.y)) >= minNonNAs){                                            ### Mehr als minNonNAs obs. in outer window (sonst NA)

# # #
    ## Location based filters
    if (any(method.names %in% c("MED","MTM"))) {
      ## Location estimation in the outer window by the median
      outer.med <- median(outer.y, na.rm=TRUE)
      
      ## (MED) Simple median filter / running median
      if ("MED" %in% method.names){ 
        results$MED[t] <- outer.med
        slopes$MED[t]  <- 0
      }
      ## (MTM) Modified Trimmed Mean filter
      if ("MTM" %in% method.names) {
        # scale estimation
        outer.loc.s        <- scale.est(na.omit(outer.y))
        outer.loc.sigma[t] <- outer.loc.s
        # outlier identification
        ol.in <- which(abs(outer.y-outer.med) > d*outer.loc.s) 
        if (length(ol.in)> 0) {
          y.mtm <- outer.y[-ol.in]
        } else {
          y.mtm <- outer.y
        }
        results$MTM[t] <- mean(y.mtm, na.rm=TRUE)
        slopes$MTM[t]  <- 0
      }
    } # end of MED MTM
# # #

# # #
    ## Level estimation by repeated median regression (based on the outer slope)
    if (any(method.names %in% c("RM","MRM","TRM"))) {
      ## Estimating the RM-slope from the outer window
      outer.slopes <- c()
      outer.slopes <- sapply(1:outer.width, function(k)
                             median((outer.y[k]-outer.y[-k])/(outer.i[k]-outer.i[-k]),
                                    na.rm=TRUE)
                             )
      if (all(is.na(outer.slopes))){                                   
        outer.slope.t <- NA                                            
      } else {
        outer.slope.t <- median(outer.slopes, na.rm=TRUE)
      }

      ## Estimating the RM intercept from the outer window based on the outer slope
      outer.mu.t <- median(outer.y - outer.slope.t*outer.i, na.rm=TRUE)
      
      ## (RM) Simple RM filter
      if ("RM" %in% method.names) {
        slopes$RM[t] <- outer.slope.t
        results$RM[t]   <- outer.mu.t
      }

      if(any(method.names %in% c("MRM","TRM"))) { 
        ## Scale estimation in the outer window based on the RM residuals
        outer.res.t <- outer.y - (outer.mu.t + outer.i*outer.slope.t)
        outer.reg.s <- scale.est(na.omit(outer.res.t))
        outer.reg.sigma[t] <- outer.reg.s
        ## Indices for outlying values (based on the calculations in the outer window)
        ol.out <- which(abs(outer.res.t) > d*outer.reg.s) 
        if (length(ol.out) > 0) {
          trimmed.i <- outer.i[-ol.out]
          trimmed.y <- outer.y[-ol.out]
        } else {
          trimmed.i <- outer.i
          trimmed.y <- outer.y
        }
      
        ## (MRM) Final RM Fit based on the trimmed values
        if ("MRM" %in% method.names) {
            tslopes <- sapply(1:length(trimmed.y), function(k) 
                            median((trimmed.y[k] - trimmed.y[-k])/(trimmed.i[k]-trimmed.i[-k]),
                                    na.rm=TRUE))
            if (all(is.na(tslopes))) { 
                slope.t <- NA
            } else {
                slope.t <- median(tslopes, na.rm=TRUE)
            }
            slopes$MRM[t]  <- slope.t
            results$MRM[t] <- median(trimmed.y - slope.t*trimmed.i, na.rm=TRUE)
        } 

        ## (TRM) Final Least Squares Fit based on the trimmed values
        if ("TRM" %in% method.names) {
            LS <- coef(lm(trimmed.y ~ trimmed.i))
            results$TRM[t] <- LS[1]
            slopes$TRM[t] <- LS[2]
        }
      } # end of MRM, TRM
    } # end of RM, MRM TRM
# # #

# # #
  ## DW-Methods: Estimations based on smaller inner window widths
  if( any(method.names %in% c("DWMTM","DWTRM","DWMRM","DWRM")) ) {
  if( sum(!is.na(inner.y)) >= minNonNAs ){                                          ### Mehr als minNonNAs obs. in inner window (sonst NA)
    
      ## (DWMTM) Modified Trimmed Mean with smaller inner window width 
      if ("DWMTM" %in% method.names) {
        # scale estimation
        inner.s   <- scale.est(na.omit(inner.y))
        inner.loc.sigma[t] <- inner.s  
        # level and slope estimation
        inner.med <- median(inner.y, na.rm=TRUE)
        ol.in <- which(abs(outer.y - inner.med) > d*inner.s)
        if(length(ol.in)> 0){
          y.dwmtm <- outer.y[-ol.in]
        } else {
          y.dwmtm <- outer.y
        }
        results$DWMTM[t] <- mean(y.dwmtm, na.rm=TRUE)
        slopes$DWMTM[t]  <- 0
      } # end of DWMTM

    if (any(method.names %in% c("DWRM","DWMRM","DWTRM"))) {
      ## Estimating the RM-slope from the inner window
      inner.slopes <- sapply(1:inner.width, function(k) {
        median((inner.y[k]-inner.y[-k])/(inner.i[k]-inner.i[-k]), na.rm=TRUE)})      
      if (all(is.na(inner.slopes))){
        inner.slope.t <- NA
      } else {
        inner.slope.t <- median(inner.slopes, na.rm=TRUE)
      }

      ## (DWRM) Double Window Repeated Median
      if ("DWRM" %in% method.names) {
        slopes$DWRM[t]  <- inner.slope.t
        results$DWRM[t] <- median(outer.y - inner.slope.t*outer.i, na.rm=TRUE)
      } # end of DWRM

      ## DWTRM, DWMRM
      if (any(method.names %in% c("DWMRM","DWTRM"))) {
        # intercept for RM regression in inner window
        inner.mu.t   <- median(inner.y - inner.slope.t*inner.i, na.rm=TRUE)
        # scale estimation based on RM residuals in inner window
        inner.res.t  <- inner.y - (inner.mu.t + inner.i*inner.slope.t)
        inner.reg.s  <- scale.est(na.omit(inner.res.t))
        inner.reg.sigma[t] <- inner.reg.s
        ## Trimming based on the RM fit in inner window
        ol.inner <- which(abs(outer.y - (inner.mu.t + outer.i*inner.slope.t) ) > d*inner.reg.s)
        if (length(ol.inner) > 0) {
            trimmed.i <- outer.i[-ol.inner]
            trimmed.y <- outer.y[-ol.inner]
        } else {
            trimmed.i <- outer.i
            trimmed.y <- outer.y
        }
    
        ## (DWMRM) Final RM Fit based on the trimmed values
        if ("DWMRM" %in% method.names) {
            tslopes <- sapply(1:length(trimmed.y),
                            function(k) {
                            median((trimmed.y[k] - trimmed.y[-k])/(trimmed.i[k]-trimmed.i[-k]),
                                    na.rm=TRUE)})
            if (all(is.na(tslopes))) {
                slope.t <- NA
            } else {
                slope.t <- median(tslopes, na.rm=TRUE)
            }
            
            slopes$DWMRM[t]   <- slope.t
            results$DWMRM[t]  <- median(trimmed.y - slope.t*trimmed.i, na.rm=TRUE)
        }
    
        ## (DWTRM) Final Least Squares Fit based on the trimmed values
        if ("DWTRM" %in% method.names) {
            LS <- lm(trimmed.y ~ trimmed.i)$coef            
            results$DWTRM[t] <- LS[1]                             
            slopes$DWTRM[t]  <- LS[2]                             
        }
    } # end of DWTRM, DWMRM

    } # end of DWTRM, DWMRM, DWRM
# # #
  } # end of 'enough' (i.e. more than minNonNAs) non-missing observations in inner window   ###
  } # end of DWMTM, DWTRM, DWMRM, DWRM
# # #
  } # end of 'enough' (i.e. more than minNonNAs) non-missing observations in outer window   ###


# # # # #
  } # end of for-loop over all time points
# # # # #

  ## Extrapolation of the end values for the signal estimates
  if (extrapolate) {
    for(m.name in method.names){
      if (m.name %in% names(slopes)) {
        temp.level <- results[[m.name]]
        temp.slope <-  slopes[[m.name]]
        if (online){
          temp.level[1:(outer.width-1)] <- temp.level[outer.width]-((outer.width-1):1)*temp.slope[outer.width]
          temp.slope[1:(outer.width-1)]  <- temp.slope[outer.width]
        } else {      
          temp.level[1:m] <- temp.level[m+1]-(m:1)*temp.slope[m+1]
          temp.slope[1:m]  <- temp.slope[m+1]
          temp.level[(N-m+1):N] <- temp.level[N-m]+(1:m)*temp.slope[N-m]
          temp.slope[(N-m+1):N]  <- temp.slope[N-m]
        }
        results[[m.name]] <- temp.level
        slopes[[m.name]]  <- temp.slope
      }
    }
    if( !all(method.names %in% c("MED","RM")) ){
        if(online){
            outer.loc.sigma[1:(outer.width-1)] <- outer.loc.sigma[outer.width]
            outer.reg.sigma[1:(outer.width-1)] <- outer.reg.sigma[outer.width]
            if(any(method.names %in% c("DWMTM","DWMRM","DWTRM"))){ 
               inner.loc.sigma[1:(outer.width-1)] <- inner.loc.sigma[outer.width]
               inner.reg.sigma[1:(outer.width-1)] <- inner.reg.sigma[outer.width]
            }
        } else {
           outer.loc.sigma[1:m]       <- outer.loc.sigma[m+1]                   
           outer.reg.sigma[1:m]       <- outer.reg.sigma[m+1]                   
           outer.loc.sigma[(N-m+1):N] <- outer.loc.sigma[N-m]                   
           outer.reg.sigma[(N-m+1):N] <- outer.reg.sigma[N-m]                   
            if(any(method.names %in% c("DWMTM","DWMRM","DWTRM"))){ 
                inner.loc.sigma[1:m]       <- inner.loc.sigma[m+1]
                inner.reg.sigma[1:m]       <- inner.reg.sigma[m+1]
                inner.loc.sigma[(N-m+1):N] <- inner.loc.sigma[N-m]
                inner.reg.sigma[(N-m+1):N] <- inner.reg.sigma[N-m]
            }
        }
    }
  }

  ## data frame with scale estimates from inner and outer time window
  sigma <- data.frame(inner.loc.sigma=inner.loc.sigma, inner.reg.sigma=inner.reg.sigma, outer.loc.sigma=outer.loc.sigma, outer.reg.sigma=outer.reg.sigma)

  return( structure( list( level=as.data.frame(results), slope=as.data.frame(slopes), sigma=sigma,
                           y=y, outer.width=outer.width, inner.width=inner.width, method=method, 
                           scale=scale, d=d, minNonNAs=minNonNAs,
                           online=online, extrapolate=extrapolate, ts.name=ts.name),
          class="dw.filter")
        )
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Default output
print.dw.filter <- function(x, ...) {

    # Number of filter methods
    n.met   <- dim(x$level)[2]

    S <- summary(x)
    N <- length(x$y)
    start <- ifelse(x$online, x$outer.width, (x$outer.width+1)/2)
    if(all(is.na(x$level[start+c(0:2),]))){ 
        start <- min(which(apply(x$level, 1, function(L){any(!is.na(L))}) )) 
    }

    # Define level and slope output
    if (n.met == 1){
      level.output <- rbind(rep("...",n.met),matrix(round(x$level[start+c(0:2),],6), ncol=1),rep("...",n.met))
      slope.output <- rbind(rep("...",n.met),matrix(round(x$slope[start+c(0:2),],6), ncol=1),rep("...",n.met))
    } else {
      level.output <- rbind(rep("...",n.met),round(x$level[start+c(0:2),],6),rep("...",n.met))
      slope.output <- rbind(rep("...",n.met),round(x$slope[start+c(0:2),],6),rep("...",n.met))  
    }
    rownames(level.output) <- c(" ", start+c(0:2), "  ")
    rownames(slope.output) <- c(" ", start+c(0:2), "  ")
    # Print level and slope output
    cat("$level \n")
    print(level.output)
    cat("('",S["level","Class"],"' with ",N," obs. of", S["level","Length"], " variables)\n \n",sep="")
    cat("$slope \n")
    print(slope.output)
    cat("('",S["slope","Class"],"' with ",N," obs. of", S["slope","Length"], " variables)\n \n",sep="")
    # Output sigma if any inner or outer scale was evaluated
    if( any(!is.na(x$sigma)) ){
        s.est <- apply(x$sigma,2,function(a) !all(is.na(a)))
        sigma.output <- rbind(rep("...",length(which(s.est))),round(x$sigma[start+c(0:2),s.est],6),rep("...",length(which(s.est))))
#        sigma.output <- rbind(rep("...",n.met),round(x$sigma[start+c(0:2),],6),rep("...",n.met))
        rownames(sigma.output)[c(1,5)] <- c(" ","  ")
        cat("$sigma \n")
        print(sigma.output)
        cat("('",S["sigma","Class"],"' with ",N," obs. of", S["sigma","Length"], " variables: ",length(which(s.est))," non-missing variables displayed)\n \n",sep="")
    }

}



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Default plot
plot.dw.filter <- function(x, ...) {
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
        t1 <- "Double Window Filters"
    } else {
        t1 <- "Double Window Filter"
    }
    titel <- ifelse(x$online,paste("Online ",t1,sep=""),t1)

    # Possible colors
    mcols <- c("red","green3","blue","skyblue2","orange","yellow","grey","purple2","darkgreen")#"olivedrab4")

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
