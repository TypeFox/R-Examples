plotPseudospectrum = function(
  ##title<< Plot and calculate the pseudospectrum of spectrally decomposed SSA eigentriples
  ssa.object            ##<< SSA object: the results of a run of ssa().
  ,calc.raw.SSA = TRUE  ##<< logical: Whether to additionally compute the whole spectrum for
                        ##   all un-grouped eigentriples (my slow the process in case of
                        ##   long time series.
  ,plot.spectrum = TRUE ##<< logical: whether to plot the pseudospectrum.
  ,plot.fourier = TRUE  ##<< logical: Whether to plot the Fourier spectrum in the background.
  ,series.orig = c()    ##<< numeric vector: original, non decomposed time series
                        ##   (used to calculate Fourier spectrum). If not supplied,
                        ##   an object with the saved in ssa.object is searched for in
                        ##   all active environments.
  ,pch = 16             ##<< integer: graphic parameter passed to plot() (?par)
  ,col = 'red'          ##<< character: graphic parameter passed to plot() (?par)
  ,show.harmonies = TRUE##<< logical: whether to mark the positions of the harmonies of the oscillation
                        ##   with the highest variance.
  ,label.points = TRUE  ##<< logical: whether to label the points with period values
  ,label.tresh = 0.0005 ##<< numeric: threshold used to label points
  ,call.freq = quote(DetermineFreq(series.t)) ##<< function to use to determine the frequencies
                        ##   of the ssa eigentriples.
  ,print.stat = TRUE    ##<< logical: whether to print status information during the
                        ## calculations.
  , ...                 ## <<  further arguments passed to the plot() call (?par)
  )
  ##description<<
  ## This function plots the pseudospectrum of the results from a SSA run, e.g. it
  ## plots the variance of the individual eigentriples vs. their frequencies. It can also
  ## be used to compute the frequency, variance and period of all SSA eigentriples.
  ##seealso<<
  ## \code{\link[Rssa]{ssa}}
{
  ## check input
  if (!inherits(ssa.object, 'ssa'))
    stop('ssa.object is not of class ssa! Did you run ssa.object <- ssa() to create it?')
  if (plot.spectrum & plot.fourier & length(series.orig) == 0) {
    if (exists(ssa.object$series)) {
      series.orig = get(ssa.object$series)
    } else {
      stop('series.orig as saved in ssa.object not found. Supply values to plot Fourier Spectrum!')
    }
  }

  ## calculate ssa 
  if(calc.raw.SSA) {
    if (print.stat)
      printStatus('Determining raw spectrum (i.e. all eigentriples)')
    rec.raw               <- reconstruct(ssa.object)
    rec.matrix.raw        <- matrix(unlist(rec.raw), byrow = TRUE,nrow = length(rec.raw))
    results.raw           <- matrix(NA, ncol = 3, nrow = length(rec.raw))
    colnames(results.raw) <- c('frequency', 'period', 'variance')
    for (h in 1:length(rec.raw)) {
      if (print.stat)
        if (h %% floor(length(rec.raw) / 10) == 0)
          printStatus(paste('....completed ', round(h / length(rec.raw) * 100, digits = 0), '%', sep = ''))
      series.t                   <- rec.matrix.raw[h, ]
      results.raw[h, 'frequency'] <- eval(call.freq)
      results.raw[h, 'variance']  <- var(series.t)
    }
    results.raw[, 'period'] <- 1/results.raw[, 'frequency']
  } else {
    results.raw = 'not computed'
  }

  ## group eigentriples
  if (print.stat)
    printStatus('Grouping eigentriples.')
  pairs.ssa          <- grouping.auto(ssa.object)

  ## reconstruct grouped triples
  if (print.stat)
    printStatus('Reconstructing grouped eigentriples.')
  rec.pairs          <- reconstruct(ssa.object, groups = pairs.ssa)

  ## compute spectrum
  if (print.stat)
    printStatus('Determining final spectrum (grouped eigentriples)')
  n.pairs            <- length(pairs.ssa)
  n.points           <- length(rec.pairs[[1]])
  rec.matrix.pairs   <- matrix(unlist(rec.pairs), byrow = TRUE, nrow = n.pairs)
  results.pairs      <- matrix(NA, ncol = 3, nrow = n.pairs)
  colnames(results.pairs) <- c('frequency','period','variance')
  for (i in 1:n.pairs) {
    if (print.stat)
      if (i %% floor(n.pairs / 5) == 0)
        printStatus(paste('....completed ', round(i / n.pairs * 100, digits = 0), '%', sep = ''))
    series.t     <- rec.matrix.pairs[i, ]
    results.pairs[i,'frequency'] <- eval(call.freq)
    results.pairs[i,'variance']  <- var(series.t)
  }
  results.pairs[, 'period'] <- 1 / results.pairs[, 'frequency']

  ## do plots
  if (plot.spectrum) {
    ## calculate and plot Forurier Spectrum
    if (print.stat)
      printStatus('Plotting results')
    new.par <- list()
    new.par$mar <- par()$mar+c(0,0,2,1)
    old.par <- par(new.par)
    try(par(new.par), silent = TRUE)
    xlim= c(1 / (2 * n.points), 0.5)

    if (plot.fourier) {
      spec.fourier <- spectrum(series.orig, plot = FALSE, span = 3)
      plot(spec.fourier$freq, spec.fourier$spec, type = 'l', log = 'xy', ylab = '',
           yaxt = 'n', xaxt = 'n', xlab = '', col = 'gray', xlim = xlim, ...)
      axis(4)
      mtext(side = 4, 'raw periodogram', line = 2)
      axis(1)
      par(new = TRUE)
    }

    ## plot SSA pseudospectrum
    if(calc.raw.SSA) {
      n.values <-n.pairs
      pch.vals <-15:18
      col.vals <-1:6
      cols.vals<-2:6
      n.cols<-length(cols.vals)
      n.pch <-ceiling(n.values / n.cols)
      pch.vec<-rep(pch.vals[1:n.pch], times = n.cols)[1:n.values]
      col.vec<-rep(col.vals[1:n.cols], each = n.pch)[1:n.values]
      col.pairs <- col.vec
      pch.pairs <- pch.vec
      triples.all    <- unlist(pairs.ssa)
      length.groups  <- sapply(pairs.ssa, length)
      group.ind.all  <- rep(as.integer(names(length.groups)), times = length.groups)
      group.ind.all.ord<-group.ind.all[order(triples.all)]
      pch.raw        <- pch.vec[group.ind.all.ord]
      pch.raw[pch.raw == 16]<-21
      pch.raw[pch.raw == 15]<-22
      pch.raw[pch.raw == 17]<-24
      pch.raw[pch.raw == 18]<-23
      col.raw        <- col.vec[group.ind.all.ord]
    } else {
      col.pairs<-rep('red', times = n.pairs)
      pch.pairs<-rep(16, times = n.pairs)
    }
    ylim<-range(c(results.raw[, 'variance'], results.pairs[, 'variance']), na.rm = TRUE)
    plot(results.pairs[, 'frequency'], results.pairs[, 'variance'], 
         ylab = 'log(variance)(SSA)', xlab = 'frequency', cex = 1.5, ylim = ylim, 
         log = 'xy', xaxt = 'n', xlim = xlim, pch = pch.pairs, col = col.pairs, yaxt = 'n', ...)
    if (calc.raw.SSA)
      points(results.raw[, 'frequency'], results.raw[, 'variance'], col = col.raw, pch = pch.raw, cex = 2)    
    axis(2)
    xaxis.values <- axTicks(side = 1)
    side.new     <- 3
    axis(side = side.new, at = xaxis.values, labels = round(1/xaxis.values, digits = 2))
    mtext(side = 3, 'period', line = 2)
    if (label.points) {
      tresh.var <- var(series.orig)*label.tresh
      x         <- results.pairs[results.pairs[, 'variance']>tresh.var, 'frequency']
      y         <- results.pairs[results.pairs[, 'variance']>tresh.var, 'variance']
      period.lab<- results.pairs[results.pairs[, 'variance']>tresh.var, 'period']
      if (!(length(y) == 0))
        text(x = x, y = y, labels = round(period.lab, digits = 1), pos = (1:length(x))%%4+1)
    }
    if (show.harmonies) {
      freq.harmonies = results.pairs[which.max(results.pairs[, 'variance']), 'frequency']*c(1:6)
      abline(v = freq.harmonies, lty = 2)
    }
    mtext(side = 3, text = 'SSA pseudospectrum', cex = 1.5, line = 4)
    legend.text<- c('paired eigentr.')
    legend.col <- c('red')
    legend.pch <- 17
    legend.lty <- NA
    if (plot.fourier) {
      legend.text<- c(legend.text, 'periodogram')
      legend.col <- c(legend.col, 'gray')
      legend.pch <- c(legend.pch, NA)
      legend.lty <- c(legend.lty, 1)
    }
    if (calc.raw.SSA) {
      legend.text<- c(legend.text, 'unpaired eigentr.')
      legend.col <- c(legend.col, 'red')
      legend.pch <- c(legend.pch, 24)
      legend.lty <- c(legend.lty, NA)
    }
    legend('topright', merge = TRUE, col = legend.col, lty = legend.lty, pch = legend.pch, 
           legend = legend.text)
    par(old.par)
  }
  
  ## prepare output
  ##value<< list with values
  results <- list(paired = results.pairs  ##<< matrix: frequency, period and variance of the paired reconstruction
                  , raw = results.raw     ##<< matrix: frequency, period and variance of the unpaired reconstruction
                  )
  invisible(results)
}
