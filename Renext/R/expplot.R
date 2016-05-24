##==============================================================================
## Author Y. Deville
##
## This function plots the sample in 'x' on an exponential plot with a little
## more material if needed.
##==============================================================================

expplot <- function(x,
                    plot.pos = "exp",
                    rate = NULL,
                    labels = NULL,
                    mono = TRUE,
                    ...) {

  if (mono) {
    col1 <- "darkgray"
  } else {
    col <- col2rgb("SteelBlue3")/256
    col1 <- rgb(col[1], col[2], col[3], 0.3)
  }
  
  xs <- sort(x)
  n <- length(x)


  if (!(plot.pos %in% c("exp", "med")))
    stop("plot.pas must be either \"exp\" or \"med\"")
      
  if (plot.pos == "medrank") F <- ((1:n) - 0.3)/(n + 0.4)
  else  F <- (1:n) / (n + 1)
  
  transF <- -log(1-F)
  
  ## Pass to log scale
  par(xlog = TRUE)
  
  plot(x = xs,
       y = transF,
       pch = 16,
       col = col1,
       yaxt = "n",
       ylab = "prob",
       ...)

  ## axis labels
  probs <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.60, 0.70, 0.80, 0.90)
  
  abline(h = -log(probs),
         lty = "dotted",
         col = "gray")
  
  abline(v = axTicks(side = 1),
         lty = "dotted",
         col = "gray")
  
  axis(side = 2,
       at =  -log(probs),
       labels = formatC(1-probs, format = "f", digits = 2))

  if (!is.null(rate)) {
    
    x.g <- seq(from = min(x), to = max(x) , length = 50)
    nW <- length(rate)

    if (mono) {
      cols <- c("black", "darkgray")
      ltys <- c("solid", "dashed", "dotted")
    } else {
      ltys <- "solid"
      cols <- c("orangered", "DarkOliveGreen3", "purple", "pink")
    }

    ltys <- rep(ltys, length.out = nW)
    cols <- rep(cols, length.out = nW)
    
    for (i in 1:nW) {

      if (!mono){
        col1 <- col2rgb(cols[i])/256
        cols[i] <- rgb(col1[1], col1[2], col1[3], 0.9)
      }
      
      transF.g <- -log(1 - pexp(x.g, rate = rate[i]))
      
      lines(x = x.g,
            y = transF.g ,
            col = cols[i],
            lty = ltys[i])
    }
    
    ## candidate position for legend.
    coords <- par()$usr
    
    if (is.null(labels)) labels <- paste("rate = ",
                                         formatC(rate, format = "f", digits = 2))
    else labels <- rep(labels, length.out = nW)
    
    legend(x = range(x)[2]*0.9,
           y = 0.60*coords[3]+0.40*coords[4],
           xjust = 1,
           yjust = 1,
           lty = ltys,
           lwd = rep(2, nW),
           col = cols,
           horiz = FALSE,
           legend = labels)
    
  }
  
}

  
