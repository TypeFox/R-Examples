rankhazardplot <- function(...) UseMethod("rankhazardplot")

rankhazardplot.default <- function(
  x, coefs = NULL, xp = NULL, refvalues = NULL, 
  refpoints = NULL, confinterval = NULL, select = 1, 
  legendtext = NULL, axistext = NULL, legendlocation = "top", 
  axistextposition = -0.1, reftick = TRUE, refline = FALSE, 
  col.refline = 1, lwd.refline = 1, lty.refline = 2, 
  ylab = NULL, ylim = NULL, yticks = NULL, yvalues = NULL,
  xtext =TRUE,  plottype = "hazard",axes = TRUE, na.rm = TRUE,
  col = NULL, lwd = 1, lty = 1, pch = NULL, cex = 1, 
  bg = "transparent", pt.lwd = 1, draw.confint = NULL,
  col.CI = col, lty.CI = lty +1, lwd.CI = lwd, add = FALSE, 
  graphsbefore = 0, args.legend = NULL, ...)				
{
  if (!identical(plottype, "hazard") && !identical(plottype, "loghazard")) 		
    stop("'plottype' must be  'hazard' or 'loghazard'.")
  if (!is.numeric(graphsbefore) || graphsbefore < 0) warning("'graphsbefore' must be a non-negative integer.")
  if (add && graphsbefore == 0) 
    warning("When 'add = TRUE' the amount of already drawn graphs must be given by 'graphsbefore'.") 
  if(!is.null(draw.confint)) warning("'draw.confint' can't be used without model object.")
  
  if(!is.null(confinterval)){
    x <- confinterval$x
    if (na.rm) x <- na.omit(x)
    x <- confinterval$x[select]
    xp <- confinterval$xp
    if (na.rm) xp <- na.omit(xp)
    xp <- confinterval$xp[select]
    refvalues <- confinterval$refvalues[select]
  }

  if (!is.data.frame(x)) stop("'x' must be a data frame.") 
  if (na.rm) x <- na.omit(x)
  
  n <- dim(x)[1]	# number of observations
  m <- dim(x)[2]	# number of covariates
  
  if (!is.null(coefs) && m != length(coefs)) 
    stop ("The length of the vector 'coefs' must be the same as the number of columns in 'x'.")
  
  if (!is.null(refpoints) && !is.null(refvalues)) 
    stop("Only one of the arguments 'refpoints' and 'refvalues' can be in use at the same time.")
  
  if (!is.null(xp)){
    if (!is.null(coefs)) stop("Only one of the arguments 'coefs' and 'xp' can be in use at the same time.")
    if (!is.data.frame(xp)) stop("'xp' must be a data frame.") 
    if (is.null(refvalues)) stop("When 'xp' is given, also 'refvalues' are required.")	
    if (na.rm) xp <- na.omit(xp)  
    if (!identical(dim(xp),dim(x))) stop("The dimensions of 'xp' and 'x' must be the same.") 
    if (m != length(refvalues)) stop("The length of the vector 'refvalues' must be the same as the number of columns in 'xp'.")
  }else{
    if (is.null(coefs)) stop("Either 'coefs' or 'xp' must be provided.")
    if (!is.null(refvalues)) stop("'refvalues' can only be used with 'xp'.")
    xp <- as.data.frame(t(coefs * t(x))) 
  }

  if(is.null(refvalues)){	
    if (is.null(refpoints)) refpoints <- apply(x, 2, median, na.rm = TRUE)	
    if (m != length(refpoints)) stop("The length of the vector 'refpoints' must be the same as the number of columns in 'x'.")
    refvalues <- coefs*refpoints
  }

  if (!is.numeric(lwd) || lwd < 0) warning("'lwd' must be a positive number.")
  lwd <- rep(lwd, length.out = m)
  lty <- rep(lty, length.out = m)
  cex <- rep(cex, length.out = m)
  bg <- rep(bg, length.out = m)
  if (!is.numeric(pt.lwd)|| pt.lwd < 0) warning("'pt.lwd' must be a positive number.")
  pt.lwd <- rep(pt.lwd, length.out = m)
  
  if (!add && graphsbefore != 0)  warning("'graphsbefore' is not zero even though a new plot is drawn.")

  if (is.null(pch)){pch <- seq(0, m - 1) + graphsbefore} 		
  else{pch <- rep(pch, length.out = m)}							
  if (is.null(col)) {col <- 1:m + graphsbefore }				
  else{ col <- rep(col, length.out = m)	}	
  
  if(is.null(col.CI)) col.CI <- col
    
  if (is.null(legendtext) && !is.null(axistext)) 
    legendtext <- axistext	
  if (!is.null(legendtext) && is.null(axistext)) 
    axistext	<- legendtext	 
  if (is.null(legendtext) && is.null(axistext) && !is.null(names(xp))) 
    legendtext <- names(xp)			
  if (is.null(legendtext) && is.null(axistext) && !is.null(names(coefs))) 
    legendtext <- names(coefs)
  if (is.null(axistext) && !is.null(colnames(x))) 
    axistext <- colnames(x)			
  if (is.null(axistext)) 
    axistext <- legendtext
  
  if (length(axistext) != m) warning("The length of the vector 'axistext' differs from the number of covariates to be drawn.")
  if (length(legendtext) != m) warning("The length of the vector 'legendtext' differs from the number of covariates to be drawn.")
  
  ones <- matrix(1, nrow = n, ncol = 1)	
  y <- xp - ones %*% refvalues

  if(!is.null(confinterval)){
    upp_ci <- confinterval$upp - ones %*% confinterval$upprefvalues
    upp_ci <- upp_ci[select]
    low_ci <- confinterval$low - ones %*% confinterval$lowrefvalues
    low_ci <- low_ci[select]
  }

  if (identical(plottype, "hazard")){
    y <- exp(y)
    if (!is.null(confinterval)){
      low_ci <- exp(low_ci) 
      upp_ci <- exp(upp_ci)  
    }
  }

  yrange <- y    # makes sure that confidence intervals fit to the screen
  if(!is.null(confinterval))
    yrange <- as.data.frame(c(y, low_ci, upp_ci))	 
    
  if (length(ylim)!= 2){ 
    if (!is.null(ylim)) warning("The length of 'ylim' differs from two, and the default is used.")
    maxy <- max(yrange, na.rm = TRUE)
    miny <- min(yrange, na.rm = TRUE) 
  }else{
    maxy <- ylim[2]
    miny <- ylim[1]
  }

  if (identical(plottype, "hazard")) {	
    if(is.null(ylab)) ylab <- "relative hazard"					
    if (is.null(yticks))
      yticks <- c(pretty(c(miny, 1)), pretty(c(1, maxy)))	
    reftickvalue <- 1
    logvar = "y"
  }    											
  if (identical(plottype, "loghazard")) {
    if(is.null(ylab)) ylab <- "logarithm of the relative hazard"										
    if (is.null(yticks))
      yticks <- c(pretty(c(miny, 0)), pretty(c(0, maxy)))
      reftickvalue <- 0	
      logvar = ""
  }

  if (is.null(yvalues)) yvalues <- yticks
    
  quantiles <- c(0, 0.25, 0.5, 0.75, 1)	
  orders <- apply(x, 2, order) 
  scaleranks <- x
  y_ord <- y
  rank_quantile <- matrix(ncol=m, nrow=5)
  y_points <- matrix(ncol=m, nrow=5)
  
  na_sum <- colSums(is.na(x))
 
  for(i in 1:m){
    scaleranks[i] <- c(seq(0, 1, length = n - na_sum[i]), rep(NA, na_sum[i]))
    y_ord[i] <- y[orders[,i],i]
    rank_quantile[,i] <-quantile(1:(n - na_sum[i]), probs = quantiles)
    y_points[,i] <- y_ord[,i][rank_quantile[,i]]
  }
  
  if (!add){
    matplot(1, 1, log=logvar, ylim=c(miny, maxy), xlim=c(0,1), ylab=ylab, xlab="", xaxt="n", yaxt = "n",type = "n",  axes = axes, ...)
  
    if (refline){# draws the reference line
      if (!is.numeric(lwd.refline)) warning("'lwd.refline' must be numeric.")
      abline(h = reftickvalue,  col = col.refline,  lty = lty.refline, lwd = lwd.refline)
    }  
  }
  matplot(x = scaleranks, y = y_ord, type = 'l', col=col, lty=lty, lwd=lwd, add = TRUE, ylab="", xlab="", xaxt="n", yaxt = "n", ...)
  matpoints(quantiles, y_points, pch=pch, col=col, cex=cex, bg=bg, lwd=pt.lwd)
    
  if (!is.null(confinterval)){
    low_ci_ord <- low_ci
    upp_ci_ord <- upp_ci
    low_ci_points <- matrix(ncol=m, nrow=5)
    upp_ci_points <- matrix(ncol=m, nrow=5)
      
    for(i in 1:m){
      low_ci_ord[i] <- low_ci[orders[ , i], i]
      upp_ci_ord[i] <- upp_ci[orders[ , i], i]
      low_ci_points[ , i] <- low_ci_ord[ , i][rank_quantile[ , i]]
      upp_ci_points[ , i] <- upp_ci_ord[ , i][rank_quantile[ , i]]
    }
      
    if (!is.numeric(lwd.CI) || lwd.CI < 0) warning("'lwd.CI' must be a positive number.")
     
    matlines(x = scaleranks, y = low_ci_ord, type = "l", col = col.CI, lty = lty.CI, lwd = lwd.CI, ...) 
    matlines(x = scaleranks, y = upp_ci_ord, type = "l", col = col.CI, lty = lty.CI, lwd = lwd.CI, ...) 
    matpoints(quantiles, low_ci_points, pch = pch, col = col.CI, cex = cex, bg = bg, lwd = pt.lwd)
    matpoints(quantiles, upp_ci_points, pch = pch, col = col.CI, cex = cex, bg = bg, lwd = pt.lwd)
  }
    
  if (xtext)
    for(i in 1:m){
      xlabels <- x[orders[rank_quantile[,i],i], i]    # quantiles for covariate i
      if (is.numeric(xlabels)) xlabels <- signif(xlabels, 3) #rounds numeric labels
        mtext(side = 1, at = c(axistextposition, quantiles),    
              adj = c(1,rep(0.5, length(quantiles))), text = c(axistext[i], as.character(xlabels)), line = i + graphsbefore)
    }
  
  if (!add){
    if (axes){
      axis(1, at = quantiles, labels = FALSE)    # marks ticks on x-axis
      axis(2, at = yticks, labels = FALSE)    # marks ticks on y-axis
      axis(2, at = yvalues, labels = as.character(yvalues))    # marks values on y-axis
      
      if (reftick)    # eboldens the reference tick
        axis(2, at = reftickvalue, labels = FALSE, lwd.ticks = 2)
    }

    if (is.null(args.legend)){
      legend(x = legendlocation, legend = legendtext, col = col, lwd = lwd, 
             pch = pch, lty = lty, bty = "n", pt.cex = cex, pt.lwd = pt.lwd, pt.bg = bg)
    }
    else {
      args.legend1 <- list(x = legendlocation, legend = legendtext, col = col, lwd = lwd, 
                           pch = pch, lty = lty, bty = "n", pt.cex = cex, pt.lwd = pt.lwd, pt.bg = bg)
      if (!identical(class(args.legend), "list")) warning("The class of 'args.legend' should be list.")
      args.legend1[names(args.legend)] <- args.legend
      do.call("legend", args.legend1)
    }
  }

  ### Output to console ####
  quantile_na.rm <- function(...) quantile(..., na.rm = TRUE, probs = quantiles)
    
  rankhazard_quantiles <- t(apply(y, 2, quantile_na.rm))
  colnames(rankhazard_quantiles) <- c("Min.", "1st Qu.", "Median" , "3rd Qu.", "Max.")
  rownames(rankhazard_quantiles) <- legendtext[1:m]

  cat("Y-axis range: ", signif(c(miny, maxy), 3), "\n", "\n")
  if (identical(plottype, "hazard")) cat("Relative hazards for each covarite:", "\n")
  if (identical(plottype, "loghazard")) cat("Logarithm of the relative hazards for each covarite:", "\n")
  print(signif(rankhazard_quantiles, 3))

  if (!is.null(confinterval)){
      
    cat("\n")
    if (identical(plottype, "hazard")) 
      cat("Relative hazards for the confidence intervals of each covariate:", "\n")
    if (identical(plottype, "loghazard")) 
      cat("Logarithm of the relative hazards for the confidence intervals of each covariate:", "\n")
      
    rankhazard_CI_quantiles <- matrix(0, 2 * m, 5)
    colnames(rankhazard_CI_quantiles) <- c("Min.", "1st Qu.", "Median" , "3rd Qu.", "Max.")
    low_legend <- paste("Low", legendtext, sep = "_")
    upp_legend <- paste("Upp", legendtext, sep = "_")
    rownames(rankhazard_CI_quantiles)[2 * 1:m] <- upp_legend
    rownames(rankhazard_CI_quantiles)[2 * 1:m - 1] <- low_legend

    rankhazard_CI_quantiles[2 * 1:m - 1, ] <- t(apply(low_ci, 2, quantile_na.rm))
    rankhazard_CI_quantiles[2 * 1:m, ] <- t(apply(upp_ci, 2, quantile_na.rm))
      
    print(signif(rankhazard_CI_quantiles, 3))
  }		
}