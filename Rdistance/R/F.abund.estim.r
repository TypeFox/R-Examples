F.abund.estim <- function(dfunc, detection.data, transect.data,
                          area=1, ci=0.95, R=500, by.id=FALSE,
                          plot.bs=FALSE){
  
  # Stop and print error if key columns of detection.data or transect.data contain NAs
  if(any(is.na(detection.data$dist))) stop("Please remove rows for which detection.data$dist is NA.")
  if(any(is.na(detection.data$siteID))) stop("Please remove rows for which detection.data$siteID is NA.")
  if(any(is.na(detection.data$groupsize))) stop("Please remove rows for which detection.data$groupsize is NA.")
  
  if(any(is.na(transect.data$siteID))) stop("Please remove NA's from transect.data$siteID.")
  if(any(is.na(transect.data$length))) stop("Please remove NA's from transect.data$length.")
  
  
  # Plotting
  f.plot.bs <- function(x, xscl, yscl, ...) {
    x.seq <- seq(x$w.lo, x$w.hi, length = 200)
    g.at.x0 <- x$g.x.scl
    x0 <- x$x.scl
    y <- like(x$parameters, x.seq - x$w.lo, series = x$series, 
              expansions = x$expansions, w.lo = x$w.lo, w.hi = x$w.hi)
    f.at.x0 <- like(x$parameters, x0 - x$w.lo, series = x$series, 
                    expansions = x$expansions, w.lo = x$w.lo, w.hi = x$w.hi)
    yscl <- g.at.x0/f.at.x0
    lines(x.seq, y * yscl, ...)
  }
  if (plot.bs) {
    tmp <- plot(dfunc) 
    x.scl.plot <- tmp$xscl.plot
    y.scl.plot <- tmp$yscl
    like <- match.fun(paste(dfunc$like.form, ".like", sep = ""))
  }
  
  # Apply truncation specified in dfunc object (including dist equal to w.lo and w.hi)
  (detection.data <- detection.data[detection.data$dist >= dfunc$w.lo & detection.data$dist <= dfunc$w.hi, ])

  # sample size (number of detections, NOT individuals)
  (n <- nrow(detection.data))
  
  # group sizes
  (avg.group.size <- mean(detection.data$groupsize))

  # total transect length and ESW
  (tot.trans.len <- sum(transect.data$length))
  (esw <- ESW(dfunc))  # get effective strip width

  # estimate abundance
  (n.hat <- avg.group.size * n * area/(2 * esw * tot.trans.len))
  
  # store output
  ans <- dfunc
  ans$n.hat <- n.hat
  ans$n <- n
  ans$area <- area
  ans$esw <- esw
  ans$tran.len <- tot.trans.len
  ans$avg.group.size <- avg.group.size



  if (!is.null(ci)) {
    # Compute bootstrap CI by resampling transects

      g.x.scl.orig <- dfunc$call.g.x.scl  # g(0) or g(x) estimate
      
      n.hat.bs <- rep(NA, R)  # preallocate space for bootstrap replicates of nhat
      
      # Turn on progress bar (if utils is installed)
      if ("utils" %in% installed.packages()[, "Package"]) {
        pb <- txtProgressBar(1, R)
        show.progress = TRUE
      } else show.progress = FALSE
      
      
      # Bootstrap
      cat("Computing bootstrap confidence interval on N...\n")
      for(i in 1:R){
        # sample rows, with replacement, from transect data
        new.transect.data <- transect.data[sample(nrow(transect.data), nrow(transect.data), replace=TRUE), ]
        
        new.trans <- as.character(new.transect.data$siteID)  # which transects were sampled?
        trans.freq <- data.frame(table(new.trans))  # how many times was each represented in the new sample?
        
        # subset distance data from these transects
        if( class(new.transect.data$siteID) == "factor" ){
          new.trans <- unique(droplevels(new.transect.data$siteID))
        } else {
          new.trans <- unique(new.transect.data$siteID)
        }
        new.detection.data <- detection.data[detection.data$siteID %in% new.trans, ]  # this is incomplete, since some transects were represented > once
        
        # replicate according to freqency in new sample
        # merge to add Freq column to indicate how many times to repeat each row
        red <- merge(new.detection.data, trans.freq, by.x="siteID", by.y="new.trans")
        # expand this reduced set my replicating rows
        new.detection.data <- red[rep(seq.int(1, nrow(red)), red$Freq), -ncol(red)]
        
        # Extract distances
        new.x <- new.detection.data$dist
        
        #update g(0) or g(x) estimate.
        if (is.data.frame(g.x.scl.orig)) {
          g.x.scl.bs <- g.x.scl.orig[sample(1:nrow(g.x.scl.orig), 
                                            replace = TRUE), ]
        } else g.x.scl.bs <- g.x.scl.orig
        
        
        # estimate distance function
        dfunc.bs <- F.dfunc.estim(new.x, likelihood = dfunc$like.form, 
                                  w.lo = dfunc$w.lo, w.hi = dfunc$w.hi, expansions = dfunc$expansions, 
                                  series = dfunc$series, x.scl = dfunc$call.x.scl, 
                                  g.x.scl = g.x.scl.bs, observer = dfunc$call.observer, 
                                  warn = FALSE)
        
        # Store ESW if it converged
        if (dfunc.bs$convergence == 0) {
          esw.bs <- ESW(dfunc.bs)
          if (esw.bs <= dfunc$w.hi) {
            
            # Calculate observed metrics
            # sample size
            n.bs <- nrow(new.detection.data)
            
            # group sizes
            avg.group.size.bs <- mean(new.detection.data$groupsize)
            
            # Store observed metrics
            #esw <- ESW(dfunc.bs)  #get effective strip width
            tot.trans.len.bs <- sum(new.transect.data$length)
            n.hat.bs[i] <- avg.group.size.bs * n.bs * area/(2 * esw.bs * tot.trans.len.bs)  # area stays same as original?   
            
          }  # end if esw.bs <= w.hi
          if (plot.bs) 
            f.plot.bs(dfunc.bs, x.scl.plot, y.scl.plot, col = "blue", lwd = 0.5)
          
          if (show.progress) setTxtProgressBar(pb, i)
        }  # end if dfunc.bs converged
      }  # end bootstrap
    
    
    # close progress bar  
    if (show.progress) close(pb)
    
    # plot red line of original fit again (over bs lines)
    if (plot.bs) f.plot.bs(dfunc, x.scl.plot, y.scl.plot, col = "red", lwd = 3)
    
    
    # Calculate CI from bootstrap replicates using bias-corrected bootstrap method in Manly text
    p <- mean(n.hat.bs > n.hat, na.rm = TRUE)
    z.0 <- qnorm(1 - p)
    z.alpha <- qnorm(1 - ((1 - ci)/2))
    p.L <- pnorm(2 * z.0 - z.alpha)
    p.H <- pnorm(2 * z.0 + z.alpha)
    ans$ci <- quantile(n.hat.bs[!is.na(n.hat.bs)], p = c(p.L, p.H))
    ans$B <- n.hat.bs
    if (any(is.na(n.hat.bs))) cat(paste(sum(is.na(n.hat.bs)), "of", R, "iterations did not converge.\n"))
    
  }  else {
    # Don't compute CI if ci is null
    ans$B <- NA
    ans$ci <- c(NA, NA)
    }  # end else

  
  
  # Compute transect-level densities
  if (by.id) {
    
    # Starting df
    nhat.df <- transect.data[, c("siteID", "length")]
    
    # Summarize raw count (truncated observations excluded previously) by transect
    rawcount <- data.frame(rawcount = tapply(detection.data$groupsize, detection.data$siteID, sum))
    rawcount <- cbind(siteID = rownames(rawcount), rawcount)
    
    # Merge and replace NA with 0 for 0-count transects
    nhat.df <- merge(nhat.df, rawcount, by="siteID", all.x=TRUE)
    nhat.df$rawcount[is.na(nhat.df$rawcount)] <- 0
    
    # Calculate transect-level abundance (density)
    nhat.df$nhat <- (nhat.df$rawcount * area) / (2 * esw * nhat.df$length)   

    # Check that transect-level abundances match total abundance
    #mean(nhat.df$nhat)
    #ans$n.hat
    
    # Remove the length column
    nhat.df <- nhat.df[, -2]
    
    # Save in output list
    ans$nhat.df <- nhat.df
  }  # end if by.id
  
  
  
  # Output
  ans$alpha <- ci
  class(ans) <- c("abund", class(dfunc))
  ans

  
}  # end function
