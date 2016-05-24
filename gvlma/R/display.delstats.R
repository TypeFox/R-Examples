"display.delstats" <-
function(deletedStatvals, deletedpvals, nsd = 3, TukeyStyle = TRUE,
         statname = "G",
         pointlabels)
  {
    n <- length(deletedStatvals)
    
    if (length(deletedpvals) != n)
      stop("deleted test statistic values must be same length as p values")
    if (missing(pointlabels)) pointlabels <- as.character(1:n)
    grange <- range(deletedStatvals) 
    prange <- range(deletedpvals) # maximum range is (0, 1)
    glim <- max(abs(grange)) * c(-1,1) * 1.1
    plim <- max(abs(prange)) * c(0,1) * 1.1
    if (missing(statname))
      {
        argname <- deparse(substitute(deletedStatvals))
        if (length(grep("G", argname))>0) statname <- "G"
        else if (length(grep("1", argname))>0) statname <- "S1"
        else if (length(grep("2", argname))>0) statname <- "S2"
        else if (length(grep("3", argname))>0) statname <- "S3"
        else if (length(grep("4", argname))>0) statname <- "S4"
        else statname <- "xx"
      }
    xlabnm = switch(statname,
      "G" = expression(paste("Deleted ", {G[4]}^2, " statistic (% change)")),
      "S1" = expression(paste("Deleted ", {S[1]}^2, " statistic (% change)")),
      "S2" = expression(paste("Deleted ", {S[2]}^2, " statistic (% change)")),
      "S3" = expression(paste("Deleted ", {S[3]}^2, " statistic (% change)")),
      "S4" = expression(paste("Deleted ", {S[4]}^2, " statistic (% change)")),
      "xx" = paste("Deleted statistics")
      )
    plot(glim, plim, type = "n",
         xlab = xlabnm,
         ylab = expression(paste("Deleted p value")),
	 main="Outlying and Influential Observtions")
    points(deletedStatvals, deletedpvals, pch = 1)
    meanG <- mean(deletedStatvals)
    sdG <- sqrt(var(deletedStatvals))
    meanp <- mean(deletedpvals)
    sdp <- sqrt(var(deletedpvals))
### If TukeyStyle, then use nsd*(inter-quartile range) to delimit
### unusual observations.  Otherwise use nsd*(standard deviation) as
### the limit.
    if (TukeyStyle)
      {
	quantsG <- quantile(deletedStatvals,probs=c(.25,.50,.75))
	iqrG <- quantsG[3]-quantsG[1]
	lfG <- quantsG[1] - nsd*iqrG
	ufG <- quantsG[3] + nsd*iqrG
	quantsp <- quantile(deletedpvals,probs=c(.25,.50,.75))
	iqrp <- quantsp[3]-quantsp[1]
	lfp <- quantsp[1] - nsd*iqrp
	ufp <- quantsp[3] + nsd*iqrp
      }
    else
      {
        ufG <- meanG + nsd * sdG
        lfG <- meanG - nsd * sdG
        ufp <- meanp + nsd * sdp
        lfp <- meanp - nsd * sdp
      }
    lfp <- max(lfp, 0)
    ufp <- min(ufp, 1)
    flagG <- ifelse((deletedStatvals > ufG) | (deletedStatvals < lfG),
                    TRUE, FALSE)
    flagp <- ifelse((deletedpvals > ufp) | (deletedpvals < lfp),
                    TRUE, FALSE)
    flag <- flagG | flagp
    if (any(flag))
      {
        points(deletedStatvals[flag], deletedpvals[flag],
               pch = 4, col= 2, lwd = 2) 
        text(deletedStatvals[flag], deletedpvals[flag],
             labels = pointlabels[flag],
             pos = 3, offset = .25)
      }
    abline(v = ufG, lty = 2, lwd = .5)
    abline(v = lfG, lty = 2, lwd = .5)
    abline(h = ufp, lty = 2, lwd = .5)
    abline(h = lfp, lty = 2, lwd = .5)
    z <- data.frame(deletedStatvals, deletedpvals)[flag,]
    rownames(z) <- pointlabels[flag]
    invisible(z)
  }

