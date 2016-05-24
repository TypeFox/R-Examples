plotCurves <- function(x, y, cyc = 1, fluo = 2:ncol(x), one.plot = FALSE, 
                       nrow = ceiling(sqrt(ncol(y))), 
                       CPP = FALSE, ...) {
  
  if(is.null(y)) {
    if (!(class(x) %in% c("matrix", "data.frame")))
      stop("If 'y' is NULL, x must be matrix or data.frame")
    y <- x[, fluo]
    x <- x[, cyc]  
  }
  
  
  testxy(x, y, length = FALSE)
  
  if(one.plot) {
    matplot(x, 
            y,
            xlab = "Cycle number",
            ylab = "Fluorescence",
            type = "l")
  } else {
    if(CPP) {
      cpp.res <- apply(y, 2, function(i) CPP(x, i)[["y.norm"]])
      #y <- apply(y, 2, normalizer)
    }
    
    if(ncol(y) %% nrow != 0) {
      new.columns <- nrow - (ncol(y) %% nrow)
      y <- cbind(y, matrix(rep(NA, new.columns*nrow(y)), 
                           ncol = new.columns))
      colnames(y)[(ncol(y) - new.columns + 1):ncol(y)] <- rep("Empty", new.columns)
      if(CPP)
        cpp.res <- cbind(cpp.res, matrix(rep(NA, new.columns * nrow(y)), 
                                         ncol = new.columns))
    }
    
    #save parameters before invoking layout
    old.oma <- par("oma")
    old.mar <- par("mar")
    old.fig <- par("fig")
    
    
    lay.matrix <- matrix(1L:(2*ncol(y)), nrow = nrow*2)
    
    layout(lay.matrix, heights = rep(c(0.2, 1), nrow*2))
    
    #divide by 2, because of fake plot with column name
    lefts <- lay.matrix[seq(2, nrow(lay.matrix), by = 2)/2, 1]
    
    bottoms <- lay.matrix[nrow(lay.matrix), ]/2
    
    par(oma = c(5, 4, 4, 2))
    par(mar = c(0, 0, 0, 0))
    
    x.lim = range(x)
    y.lim = range(y, na.rm = TRUE)
    if(CPP)
      y.lim <- range(cbind(y, cpp.res), na.rm = TRUE)
    
    
    curve.names <- colnames(y)
    sapply(1L:ncol(y), function(i) {
      plot(0, 0, xlab = "", ylab = "", type = "n", xaxt = "n", yaxt = "n")
      res.NA <- which(is.na(y[, i]))
      bg.color <- ifelse(sum(res.NA) > 0, bg <- 2, bg <- 3)
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
           col = bg.color)
      text(0, 0, curve.names[i])
      
      plot(x, y[, i], xlim = x.lim, ylim = y.lim,  
           xaxt = "n", yaxt = "n", ylab = "", xlab = "", ...)
      rug(res.NA, col = 2, lwd = 1.5, quiet = TRUE)
      
      if(i %in% lefts)
        axis(2)
      if(i %in% bottoms)
        axis(1)
      if(CPP)
        lines(x, cpp.res[, i], col = "red")
    })
    
    par(oma = old.oma)
    par(mar = old.mar)
    par(fig = old.fig)
    invisible()
  }
}