plot.ntest <- function(x, ...) {
  
  ## definitions
  ## -----------
  breaks <- seq(0, 1, 0.025)
  
  if (x$method == "niche identity test"){
    
    # divide screen
    par(mfrow = c(1, 2))
    
    ## Schoener's D
    ## ------------
    d <- x$statistic[1]
    par(mar = c(5, 3, 2, 1))
    f <- hist(x$null.distribution[, 1], breaks = breaks, xlim = c(0, 1), col = "grey80", main = NULL, 		xlab = "D")
    yy <- 0.5 * max(f$c)
    lines(rep(d, 2), c(0, yy), col = "red")
    text(d, yy * 1.1, round(d, digits = 3), adj = 0.5)
    
    ## modified Hellinger distance I
    ## -----------------------------
    i <- x$statistic[2]
    par(mar = c(5, 1, 2, 1))
    f <- hist(x$null.distribution[, 2], breaks = breaks, xlim = c(0, 1), col = "grey80", main = NULL, xlab = "I")
    yy <- 0.5 * max(f$c)
    lines(rep(i, 2), c(0, yy), col = "red")
    text(i, yy * 1.05, round(i, digits = 3), adj = 0.5)
    # add title
    par(mfrow = c(1, 1))
    title(paste("Niche identity test:", paste(x$spec, 		collapse = " - ")))  
  }
  if (x$method == "background similarity test"){
    
    breaks <- seq(0, 1, 0.025)
    
    par(mfcol = c(2, 2))
    
    ## histograms for D
    ## ----------------
    d <- x$statistic[1]
    
    ## X versus random Y
    par(mar = c(1, 5, 5, 0))
    f <- hist(x$nd.x.randomY[, 1], breaks = breaks, xlim = c(0, 1), col = "grey80", 
              main = NULL, ylab = paste(x$spec[1], "vs. bg of", x$spec[2]), ...)
    yy <- 0.5 * max(f$c)
    lines(rep(d, 2), c(0, yy), col = "red")
    text(d, yy * 1.1, round(as.numeric(d), digits = 3), adj = 0.5)
    
    ## Y versus random X
    par(mar = c(5, 5, 2, 0))
    f <- hist(x$nd.y.randomX[, 1], breaks = breaks, xlim = c(0, 1), col = "grey80", 		
              main = NULL, xlab = "D", ylab = paste(x$spec[2], "vs. bg of", x$spec[1]), ...)
    yy <- 0.5 * max(f$c)
    lines(rep(d, 2), c(0, yy), col = "red")
    text(d, yy * 1.1, round(as.numeric(d), digits = 3), adj = 0.5)
    
    ## histograms for I
    ## ----------------
    i <- x$statistic[2]

    ## X versus random Y
    par(mar = c(1, 2, 5, 2))
    f <- hist(x$nd.x.randomY[, 2], breaks = breaks, xlim = c(0, 1), col = "grey80", 		
              main = NULL, ...)
    yy <- 0.5 * max(f$c)
    lines(rep(i, 2), c(0, yy), col = "red")
    text(i, yy * 1.1, round(as.numeric(i), digits = 3), adj = 0.5)
    
    ## Y versus random X
    par(mar = c(5, 2, 1, 2))
    f <- hist(x$nd.y.randomX[, 2], breaks = breaks, xlim = c(0, 1), col = "grey80", 		
              main = NULL, xlab = "I", ...)
    yy <- 0.5 * max(f$c)
    lines(rep(i, 2), c(0, yy), col = "red")
    text(i, yy * 1.1, round(as.numeric(i), digits = 3), adj = 0.5)
    
    
    par(mfrow = c(1, 1))
    par(mar = c(0, 0, 4, 0))
    title(paste("Background similarity test:", paste(x$spec, 		
                                      collapse = " - "))) 
  }
}