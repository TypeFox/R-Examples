bplot<-function (y, x, bar.col = "gray", loc.meas = mean, sort = FALSE, order = NULL, int = "SE", conf = 0.95, 
    uiw = NULL, liw = NULL, sfrac = 0.1, slty = 1, scol = 1, slwd = 1, exp.fact = 1.5, simlett = FALSE, lett.side = 3,  
    lett = NULL, cex.lett = 1, names.arg = NULL, ylim = NULL, horiz = FALSE, xpd = FALSE,...) 
{



    SE <- tapply(y, x, function(x) {
    ci.mu.t(x[!is.na(x)])$SE
    })
    CI <- tapply(y, x, function(x) {
    ci.mu.t(x[!is.na(x)], conf = conf)$margin
    })
    iqr <- tapply(y, x, IQR)
    n <- as.numeric(summary(as.factor(x)))
    iqr.ci <- 1.58 * iqr/sqrt(n)
    MAD <- tapply(y, x, function(x) {
    mad(x[!is.na(x)])})
    
	if(int == "bootSE"){
    	     lvl <- levels(factor(x))
		bootSE <- 1:length(lvl)
		for(i in 1: length(lvl)){
		bootSE[i] <- bootstrap(y[x==lvl[i]], loc.meas)$res[4]
		}
	}	
     	
    loc.vec <- tapply(y, x, function(x) {
        loc.meas(x[!is.na(x)])})
    
      if (sort == TRUE) {
      o <- order(loc.vec)
      loc.vec <- loc.vec[o]
      SE <- SE[o]
      CI <- CI[o]
      iqr <- iqr[o]
      iqr.ci <- iqr.ci[o]
      MAD <- MAD[o]
      names.arg <- names.arg[o]
	if(!is.null(lett))lett <- lett[o]	
	}
    
      if(!is.null(order)){
	 if(length(order) != nlevels(factor(x))) stop("order must be a vector whose length is equal to the number of factor levels")   
	 loc.vec <- loc.vec[order] 
	 SE <- SE[order]
      CI <- CI[order]
      iqr <- iqr[order]
      iqr.ci <- iqr.ci[order]
      MAD <- MAD[order]
      names.arg <- names.arg[order]
	if(!is.null(lett))lett <- lett[order]	
	}
	 
    if(int == "CI") margin <- CI
    if(int == "SE") margin <- SE
    if(int == "IQR") margin <- iqr
    if(int == "IQR.CI") margin <- iqr.ci
    if(int == "MAD") margin <- MAD
    if(int == "bootSE") margin <- bootSE 
    if(is.null(uiw)) uiw <- loc.vec + margin; liw <- loc.vec - margin
    
                if (simlett == TRUE & is.null(ylim)){ 
                  ylim = c(min(c(0, loc.vec - (margin * exp.fact))), 
                  max(c(0, loc.vec + (margin * exp.fact))))}
                if(simlett == FALSE & is.null(ylim)){
                  ylim <- c(min(0, loc.vec - (margin)), max(0, loc.vec + (margin)))}
                                     
                if(horiz == FALSE) b <- barplot(loc.vec, ylim = ylim , col = bar.col, names.arg = names.arg, xpd = xpd, ...)
                if(horiz == TRUE) b <- barplot(loc.vec, xlim = ylim , col = bar.col, horiz = TRUE, names.arg = names.arg, xpd = xpd, ...)
                if(horiz == FALSE){
                    arrows(b, liw, b, uiw, angle = 90, col = scol, lty = slty, lwd = slwd, length = sfrac)
                    arrows(b, liw, b, uiw, code = 1, angle = 90, col = scol, lty = slty, lwd = slwd, length = sfrac)}
                if(horiz == TRUE){
                    arrows(liw, b, uiw, b, angle = 90, col = scol, lty = slty, lwd = slwd, length = sfrac)
                    arrows(liw, b, uiw, b, code = 1, angle = 90, col = scol, lty = slty, lwd = slwd, length = sfrac)}
			if(simlett == TRUE){
                    mtext(lett, side = lett.side, cex = cex.lett, at = b, line = ifelse(lett.side == 3, 0.5, -0.5), las = ifelse(horiz == TRUE, 2, 1))}
			if(int != "CI"){
			    cat(paste("\n","Bars are ", deparse(substitute(loc.meas)), "s.  Errors are ", int, "s.", "\n\n", sep = ""))}
			if(int == "CI"){
			    cat(paste("\n", "Bars are ", deparse(substitute(loc.meas)), "s.  Errors are ", conf * 100, "% confidence intervals for the true mean.", "\n\n", sep = ""))}
}				