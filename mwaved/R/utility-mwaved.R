#' @name summary.mWaveD
#' @title Summary Output for the mWaveD object
#'
#' @param object A mWaveD object which is a list containing all the information for a multichannel 
#' deconvolution analysis produced by the \code{\link{multiWaveD}} function.
#' @param ... Arguments to be passed to methods.
#' @description Gives some numerical summaries of a \code{mWaveD} object.
#' @return Text output giving summary information of the input and output analysis including, \itemize{
#' \item Degree of Meyer wavelet used in the analysis.
#' \item Number of observations, within each channel and number of channels present.
#' \item Resolution levels used (j0 to j1)
#' \item Blur type assumed in the analysis (direct, smooth or box.car)
#' \item Matrix summarising the noise levels in each channel (and Fourier decay information for the smooth case)
#' \item Summaries of the severity of the thresholding applied amongst the resolutions.
#' }
#' @seealso \code{\link{multiWaveD}}
#' 
#' @examples
#' library(mwaved)
#' # Simulate the multichannel doppler signal.
#' m <- 3
#' n <- 2^10
#' t <- (1:n)/n
#' signal <- makeDoppler(n)
#' # Create multichannel version with smooth blur
#' shape <- seq(from = 0.5, to = 1, length = m)
#' scale <- rep(0.25, m)
#' G <- gammaBlur(n, shape, scale)
#' X <- blurSignal(signal, G)
#' # Add noise with custom signal to noise ratio
#' SNR <- c(10,15,20)
#' E <- multiNoise(n, sigma = sigmaSNR(X, SNR), alpha = c(0.5, 0.75, 1))
#' # Create noisy & blurred multichannel signal
#' Y <- X + E
#' mWaveDObject <- multiWaveD(Y, G)
#' summary(mWaveDObject)
#' @export
summary.mWaveD <- function(object, ...){
  n <- length(object$estimate)
  m <- dim(object$signal)[2]
  cat("Degree of Meyer wavelet =", object$degree, "  , Coarse resolution level j0 =", object$j0)
  cat("\n")
  cat("Sample size per channel = ", n, ", Maximum possible resolution level = ", log2(n) - 1, ".", sep = '')
  cat("\n\n")
  cat("Number of channels: m =", m,"\n")
  cat('Detected Blur Type:',detectBlur(object$G), '\n\n')
  cat('Resolution selection method: ',object$resolution,'\n\n')
  cat("Estimated Channel information:\n\n")
  
  if (object$blurDetected == "direct" && object$resolution == "smooth") {
    mat <- cbind(round(object$sigma, 3), round(object$alpha, 3), object$blurInfo$freq, rep(object$j1, m))
    colnames(mat) <- c("Sigma.hat", "Alpha", "Fourier number cutoff", "Highest resolution")
    rownames(mat) <- paste("Channel ", 1:m,':', sep='')
    print(mat, ...)
  } else {
    # If Smooth blur is used, display the matrix of values
    if (object$resolution == "smooth"){
      mat <- cbind(round(object$sigma, 3), round(object$alpha, 3), object$blurInfo$freq, object$blurInfo$maxLevels)
      colnames(mat) <- c("Sigma.hat", "Alpha", "Fourier number cutoff", "Highest resolution")
      rownames(mat) <- paste("Channel ", 1:m,':', sep='')
      print(mat, ...)
      cat("\n")
      cat("Estimated best channel = Channel", object$blurInfo$bestChannel)      
    } else {
      if (object$resolution == "block"){
        mat <- cbind(round(object$sigma, 3), round(object$alpha, 3))
        colnames(mat) <- c("Sigma.hat", "Alpha")
        rownames(mat) <- paste("Channel ", 1:m,':', sep='')
        print(mat, ...)
      } else {
        warning('Unrecognised resolution selection method.')
      }
    } 
  }
  cat("\n\n")
  cat("mWaveD optimal finest resolution level j1 =", object$j1)

  cat("\n\n")
  cat("Thresholding method:", object$shrinkType, "   Tuning parameter: eta =", object$eta,'\n\n')
  
  threshMatrix <- cbind(round(object$levelMax,4), round(object$thresh,4), object$percent)
  rownames(threshMatrix) <- paste("Level", object$j0:object$j1,":" )
  colnames(threshMatrix) <- c("Max|w|", "Threshold", "% Shrinkage" )
  print(threshMatrix, ...)
} 

#' @name plot.waveletCoef
#' @title Multi-Resolution Analysis plot of wavelet coefficients
#'
#' @description Plots the wavelet coefficient object in the multiresolution analysis
#' 
#' @param x A list of class waveletCoef.
#' @param y An optional numeric vector of trimmed wavelet coefficients to be overlayed on top of the plot for comparison with the \code{x} wavelet coefficients. 
#' @param labels Optional character vector with two elements to give name labels to \code{x} and \code{y} respectively.
#' @param ... Arguments to be passed to methods.
#' @param lowest Specifies the coarsest resolution to display in the Multi-resolution plot.
#' @param highest Specifies the finest resolution to display in the Multi-resolution plot.
#' @param scaling A numeric value that acts as a graphical scaling parameter to rescale the wavelet coefficients in the plot. A larger scaling value will reduce the size of the coefficients in the plot.
#' @param ggplot A logical value to specify if the user wants to use base graphics (FALSE) or ggplot2  graphics (TRUE).
#' 
#' @seealso \code{\link{multiCoef}} for generating a list of class `waveletCoef`
#' 
#' @export
plot.waveletCoef <- function(x, y = NULL, labels = NULL,  ..., lowest = NULL, highest = NULL, scaling = 1, ggplot = TRUE){
  stopifnot(class(x) == "waveletCoef")
  if (!is.null(y) && class(y) != "waveletCoef") {
    stop('y must be a waveletCoef object')
  }
  n <- length(x$coef)
  J <- floor(log2(n))
  fine <- ceiling(J) - 1
  # Check resolution ranges
  if (is.null(lowest)) {
    lowest <- x$j0
  } else if (lowest < x$j0 ) {
      warning("lowest level shouldn't be smaller than j0 specified in wavelet coefficient object.")
  }
  
  # Check resolution levels aren't empty.
  ind <- which.max(rev(x$coef) != 0)
  lastres <- floor(log2(n - ind + 1)) - 1
  if (is.null(highest)) {
    highest <- lastres
  } else if (highest > fine) {
      warning(paste('highest level too high. Resetting highest level to the maximum at j1 = ', fine))
  } else if (highest < lowest) {
      warning('highest level must be higher than the lowest level.')
      highest <- lowest
  }
  
  js <- rep(lowest:highest, 2^(lowest:highest))
  ks <- unlist(lapply(lowest:highest, function(i) 0:(2^i-1)/2^i))
  wi <- (2^lowest + 1):2^(highest + 1)
  nw <- length(wi)
  w  <- x$coef[wi]*scaling
  wf <- 2.05 * max(abs(w))/scaling
  w  <- w/wf
  ws <- w + js
  
  # check shrink input is a waveletCoef object
  if (!is.null(y)) {
    if (length(x$coef) != length(y$coef)) {
      stop('length of y coefficients is different to the length of x coefficients')
    }
    j0Trim <- y$j0
    if (j0Trim != x$j0) {
      warning('y object has a different coarse resolution j0 than x object, interpret lower resolution levels with caution')
    }
    wShrink <- y$coef[wi]/wf
    survived <- which(wShrink != 0)
    kss <- ks[survived]
    jss <- js[survived]
    wss <- wShrink[survived] + jss
    ns <- 2
    nss <- length(survived)
  } else {
    ns <- 1
  }
  
  mraTitle <- 'MRA'
  mraLabels <- c("Location", "Resolution Level")
  if (!is.null(labels)) {
    if (length(labels) != ns) {
      warning('length of labels might not be long enough.')
      labels = c('x', 'y')
    }
  } else {
    labels = c('x', 'y')
  }
  if (ggplot) {
    ggAvailable <- requireNamespace("ggplot2", quietly = TRUE)
    if (!ggAvailable) {
      ggplot = FALSE
    }
  }
  
  if (ggplot) {
    if (ns == 2) {
      nData <- data.frame(w = c(ws,wss), js = c(js, jss), ks = c(ks, kss), col = rep(labels, c(nw, nss)))
      mraPlot <- ggplot2::ggplot(nData) + ggplot2::geom_segment(ggplot2::aes(x = ks, xend = ks, y = js, yend = w, colour = col, size = col)) + ggplot2::labs(x = mraLabels[1], y = mraLabels[2]) + ggplot2::scale_size_discrete(range = c(1, 2))
      # Fix legend
      mraPlot <- mraPlot + ggplot2::theme(legend.position = "top", axis.text.y = ggplot2::element_text(angle = 90)) + ggplot2::guides(colour = ggplot2::guide_legend(title = mraTitle), size = ggplot2::guide_legend(title = mraTitle))
    } else {
      nData <- data.frame(w = ws, js = js, ks = ks)
      mraPlot <- ggplot2::ggplot(nData) + ggplot2::geom_segment(ggplot2::aes(x = ks, xend = ks, y = js, yend = w), colour = 'red') + ggplot2::ggtitle(mraTitle)
    }
    mraPlot + ggplot2::labs(x = mraLabels[1], y = mraLabels[2]) + ggplot2::scale_y_continuous(breaks = lowest:highest) 
  } else {
    buf <- 0.5
    plot(0, type = "n", xlim = c(0,1), ylim = c(lowest - buf, highest + buf), yaxt = 'n', xlab = mraLabels[1], ylab = mraLabels[2], main = mraTitle)
    axis(2, at = lowest:highest)
    if (!is.null(y)) {
      col <- 'red'
    } else {
      col <- 1
    }
    segments(ks, js, ks, ws, col = col)
    abline(h = lowest:highest, v = axTicks(1), col="gray", lty=3)
    if (!is.null(y)) {
      segments(kss, jss, kss, wss, lwd = 2, col = 'blue')
    }
  }
}

#' @name plot.mWaveD
#' @title Plot Output for the mWaveD object
#' 
#' @description Creates plot output that summarises the \code{mWaveD} object produced by the \code{\link{multiWaveD}} function. 
#'  
#' @param x A mWaveD object to be plotted (list created by \code{\link{multiWaveD}})
#' @param ... Arguments to be passed to methods.
#' @param which A numeric vector that specifies which plots to output. Default value is \code{1:4} which  specifies that all four plots are to be displayed.
#' @param ask A logical value that specifies whether the user is \emph{ask}ed before each plot output.
#' @param singlePlot A logical value that controls whether all plots should appear on a single window. The plot window is resized depending on the value of \code{which}.
#' @param ggplot A logical value to specify if the user wants to use base graphics (FALSE) or ggplot2 graphics (TRUE).
#' 
#' @details Four plots are output that summarise the multichannel input, a visualisation of the characteristics of the channels and the output estimate and a multi-resolution analysis plot.\itemize{
#' \item Plot 1: Multichannel input signal overlayed.
#' \item Plot 2: Estimated output signal using the mWaveD approach.
#' \item Plot 3: Plot of the log decay of Fourier coefficients against the log bounds (direct and smooth case) or the blockwise resolution levels against their limit (box car case)
#' \item Plot 4: Multi-resolution plot of the raw wavelet coefficients and the trimmed wavelet coefficients}
#' @references
#' Kulik, R., Sapatinas, T. and Wishart, J.R. (2014) \emph{Multichannel wavelet deconvolution with long range dependence. Upper bounds on the L_p risk}  Appl. Comput. Harmon. Anal. (to appear in).
#' \url{http://dx.doi.org/10.1016/j.acha.2014.04.004}
#' 
#' Wishart, J.R. (2014) \emph{Data-driven wavelet resolution choice in multichannel box-car deconvolution with long memory}, Proceedings of COMPSTAT 2014, Geneva Switzerland, Physica Verlag, Heidelberg (to appear)
#' @seealso \code{\link{multiWaveD}}
#' @importFrom grid grid.newpage pushViewport viewport grid.layout
#' @export
plot.mWaveD <- function(x, ..., which = 1L:4L, singlePlot = TRUE, ask = !singlePlot, ggplot = TRUE){
  # Check if ggplot is available if requested
  if (ggplot) {
    ggAvailable <- requireNamespace("ggplot2", quietly = TRUE)
    # Check optional dependency
    if (ggAvailable) {
      hsize <- 1
      lsize <- 0.5
      asize <- 0.5
      # Initialise list
      ggList <- list(NULL)
      i <- 1
    } 
  } else {
    ggAvailable <- FALSE
  }
  
  show <- rep(FALSE, 4)
  # Make sure which argument is numeric
  if (!is.numeric(which)) {
    stop('`which` argument must be a vector containing numerical elements')
  }
  # Only consider the integer bits between 1:4
  which <- which[which %in% 1L:4L]
  # Complain if there are no values in 1:4
  if (length(which) == 0) {
    stop('`which` argument must be a vector containing elements: 1, 2, 3 or 4')
  }
  show[which] <- TRUE
  
  n  <- length(x$estimate)
  n2 <- n/2
  m  <- dim(x$signal)[2]
  t  <- (1:n)/n
  
  blurInfo <- x$blurInfo
  resolution <- x$resolution
  j0 <- x$j0
  j1 <- x$j1
  
  estimateTitle <- 'mWaveD estimate'
  signalTitle <- 'Input Signal'
  fourierLabel <- 'Fourier number'
  fourierTitle <- 'Kernel decay'
  blockTitle <- 'Block wise resolution selection'
  mraLabels <- c('raw',paste(x$shrinkType, ' thresholded', sep = ''))
  mraTitle <- 'Multiresolution Analysis'
  
  if (show[3L]) {
    if (resolution != "block") {
      xw = fourierWindow(n)
      blur <- mirrorSpec(blurInfo$decay)
      cut <- mirrorSpec(blurInfo$cutoff)
      ylim <- c(min(blurInfo$cutoff[2, ]), 0)
      if (x$blurDetected == 'direct') {
        xlim = c(-n2, n2)
      } else {
        xbest <- max(blurInfo$freqCutoffs) - 1
        ybest <- cut[n/2 + xbest, blurInfo$bestChannel]
        xlim <- min(2*max(blurInfo$freqCutoffs), n/2)
        xlim <- c(-xlim, xlim)
      }
    } else {
      J    <- floor(log2(n)) - 1
      j    <- j0:min(c(J, 2 * j1))
      blkV <- blurInfo$blockVar[1:length(j)]
      blkc <- blurInfo$blockCutoff[1:length(j)]
      ylim <- range(c(blurInfo$blockVar, blurInfo$blockCutoff))
    }
  }
  
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  
  if (ggplot) {
    if (show[1L]) {
      signalData <- data.frame(Y = as.vector(x$signal), x = rep(t, m), Channel = rep(LETTERS[1:m], each = n))
      signalPlot <- ggplot2::ggplot(signalData, ggplot2::aes_string(x = 'x', y = 'Y', colour = 'Channel')) + ggplot2::geom_line(size = lsize, alpha = asize) + ggplot2::ggtitle(signalTitle) + ggplot2::labs(x = '', y = '')
      ggList[[i]] <- signalPlot
      i <- i + 1
    }
    if (show[2L]) {
      estimateData <- data.frame(Y = as.vector(x$estimate), x = t)
      estimatePlot <- ggplot2::ggplot(estimateData, ggplot2::aes_string(x = 'x', y = 'Y')) + ggplot2::geom_line(size = lsize, alpha = asize) + ggplot2::ggtitle(estimateTitle) + ggplot2::labs(x = '', y = '')
      ggList[[i]] <- estimatePlot
      i <- i + 1
    }
    if (show[3L]) {
      if (resolution != 'block') {
        fourierData <- data.frame(Y = as.vector(blur), x = rep(xw,m), Ycut = as.vector(cut), Channel=rep(LETTERS[1:m],each=n), m = m)
        resolutionPlot <- ggplot2::ggplot(fourierData) + ggplot2::geom_line(ggplot2::aes_string(x = 'x', y = 'Y', colour = 'Channel', group = 'Channel'),size = 1) + ggplot2::geom_line(ggplot2::aes_string(x = 'x', y = 'Ycut', colour = 'Channel'), linetype='dashed', size = 1) + ggplot2::ggtitle(fourierTitle) + ggplot2::labs(x = fourierLabel, y = '') + ggplot2::coord_cartesian(xlim = xlim)
        if (resolution == 'smooth' && x$blurDetected != 'direct') {
          rightLine <- ggplot2::geom_line(ggplot2::aes_string(x = 'x', y = 'y'), linetype = 'dotted', data = data.frame(x = rep(xbest,2), y = c(ybest, -Inf)))
          leftLine <- ggplot2::geom_line(ggplot2::aes_string(x = 'x', y = 'y'), linetype = 'dotted', data = data.frame(x = rep(-xbest,2), y = c(ybest, -Inf)))
          pointDots <- ggplot2::geom_point(ggplot2::aes_string(x = 'xbest', y = 'ybest'), shape = 1, size = 4, data = data.frame(xbest = c(-xbest, xbest), ybest = rep(ybest, 2)))
          resolutionPlot <- resolutionPlot + leftLine + rightLine + pointDots
        }
      } else {
        resolutionData <- data.frame(Y = c(blkV, blkc), x = rep(j,2), colour = rep(c("Resolution var.",'Resolution bounds'), each = length(j)) , Ycut = blkc)
        bestV <- blkV[j == j1]
        highlightData <- data.frame(x = c(j1, j1), y = c(ylim[1], bestV))
        pointData <- data.frame(j1 = j1, bestV = bestV)
        resolutionPlot <- ggplot2::ggplot(resolutionData) + ggplot2::geom_line(ggplot2::aes_string(x = 'x', y = 'Y', colour = 'colour', linetype = 'colour'), size = hsize) +  ggplot2::geom_line(ggplot2::aes_string(x = 'x', y = 'y'), linetype = 'dotted', data = highlightData) + ggplot2::labs(x = 'j', y = '') + ggplot2::geom_point( ggplot2::aes_string(x = 'j1', y = 'bestV'), size = 4, shape = 1, data = pointData)  + ggplot2::scale_color_discrete(labels= c('Resolution bounds', 'Resolution var.'), guide=ggplot2::guide_legend(title.position='left',title.theme = ggplot2::element_text(size=15,angle=0))) + ggplot2::scale_size(guide='none') + ggplot2::guides(colour = ggplot2::guide_legend( title='Blockwise resolution decay')) + ggplot2::theme(legend.position="top", legend.key = ggplot2::element_rect(fill = NA), axis.text.y = ggplot2::element_text(angle = 90)) + ggplot2::scale_linetype_manual(values=c(1,2), name="Blockwise resolution decay", labels=c('Resolution bounds', 'Resolution var.')) + ggplot2::scale_x_continuous(breaks = j)
      }
      ggList[[i]] <- resolutionPlot
      i <- i + 1
    }
    
    if (show[4L]) {
        mraPlot <- plot(x$coef, x$shrinkCoef, highest = j1, labels = c('Raw', paste('Thresholded (', x$shrinkType, ')', sep = '')), ggplot = TRUE)
        ggList[[i]] <- mraPlot
    }
    # Plot them
    if (singlePlot) {
      nPlots <- sum(show)
      if (nPlots < 4) {
        cols <- 1
      } else {
        cols <- 2
      }
      layout <- matrix(seq(1, cols * ceiling(nPlots/cols)),
                       ncol = cols, nrow = ceiling(nPlots/cols),
                       byrow = TRUE)
      if (nPlots == 1) {
        print(ggList[[1]])
      } else {
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        for (i in 1:nPlots) {
          matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
          print(ggList[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                           layout.pos.col = matchidx$col))
        }
      }
    } else {
      if (show[1L]) {
        print(signalPlot)
      }
      if (show[2L]) {
        print(estimatePlot)
      }
      if (show[3L]) {
        print(resolutionPlot)
      }
      if (show[4L]) {
        print(mraPlot)
      }
    }
  } else {
    if (singlePlot) {
      plotDims <- switch(sum(show), c(1, 1), c(2, 1), c(3, 1), c(2, 2))
      par(mfrow = plotDims)
    } else {
      par(mfrow = c(1,1))
    }
    
    if (show[1L]) {
      matplot(t, x$signal, type = 'l', main = signalTitle, ylab = '', xlab = '', lty = 1, cex = 0.8)
      grid()
    }
    
    if (show[2L]) {
      # Plot mWaveD estimate
      plot(t, x$estimate, type = 'l', main = estimateTitle, ylab = '', xlab = '', ...)
      grid()
    }
    
    if (show[3L]) {
      # Plot resolution analysis
      if (resolution != 'block') {
        iw = fourierWindow(n)
        matplot(iw, blur, type = 'l', lty = 1, xlim = xlim, ylim = ylim, main = fourierTitle, xlab = fourierLabel, ylab = "")
        matlines(iw, cut, lty = 2)
        grid()      
        if (resolution == 'smooth' && x$blurDetected != "direct") {
          points(xbest, ybest, col='blue')
          points(-xbest, ybest, col = 'blue')
          xbest <- rep(xbest, 2)
          ybest <- c(ylim[1], ybest)
          lines(xbest, ybest, lty = 'dotted')
          lines(-xbest, ybest, lty = 'dotted')
        }
      } else {
        rang = range(as.vector(c(blkV, blkc)))
        buf = 0.1 * diff(rang)
        ylims = c(rang[1] - buf, rang[2] + buf)
        plot(j, blkV, type = 'b', xlab = 'j', ylab = '', main = blockTitle, ylim = ylims)
        lines(j, blkc, col = 2)
        points(j1, blurInfo$blockVar[j == j1], col='blue')
        lines(c(j1, j1), c(ylims[1], blurInfo$blockVar[j == j1]), lty = 'dashed')
        grid()
      }
    }
    if (show[4L]) {
      plot(x$coef, x$shrinkCoef, highest = j1, ..., ggplot = FALSE)
    }
  }
}

# Auxiliary function to obtain domain for plot as a func of Fourier freq
mirrorSpec <- function(x) {
  if (is.matrix(x)){
    n <- dim(x)[1]
    x <- rbind(as.matrix(x[(n-1):2, ]), x)
  } else {
    n <- length(x)
    x <- c(x[(n-1):2], x)
  }
  x
}

# Auxiliary function to obtain domain for plot as a func of Fourier freq
fourierWindow <- function(n) {
  n2 <- floor(n/2)
  iw = -(n2 - 1):n2
  iw
}

#' @name mWaveDDemo
#' @title Interactive Demonstration
#' @description Interactive Demonstration
#' @importFrom shiny runApp
#' @export
mWaveDDemo <- function (){
  runApp(system.file('mWaveDDemo', package = 'mwaved'))
}