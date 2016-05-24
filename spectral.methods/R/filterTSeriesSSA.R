filterTSeriesSSA <- structure(function(
##title<< Decompose a vector (i.e. time series) into spectral bands
series                      ##<< numeric vector: Input time series (no gaps!)
, borders.wl                ##<< list of numeric vectors: Borders of the different periodicity
                            ## bands to extract. Units are the sampling frequency of the series
                            ## (needs one  vector per step (see details)).
, M = rep(floor(length(series) / 3), times = n.steps)  ##<< integer (vector): Window length
                            ## or embedding dimension (see details and ?ssa) (in ssa() this parameter
                            ## is called L).
, n.comp = rep(40, times = n.steps)##<< integer (vector): Amount of SSA components to compute.
                            ## See the help of ssa (package Rssa) for details.
, harmonics=rep(0, times = n.steps)##<< integer (vector): How many harmonics to include in each
                            ## component (see details). No harmonics are used if set to 0 (default).
, tolerance.harmonics = 0.05##<< numeric fraction (0-1): Tolerance to use to determine the
                            ## width of the bands the harmonics are looked for in. The actual
                            ## width is calculated by multiplying the frequency of the "main"
                            ## oscillation with this ratio. Use higher values for oscillations
                            ## with few repetitions (and, hence, wider peaks in a spectrum) and lower
                            ## ones with those with many repetitions (and thus sharper peaks).
, var.thresh.ratio = 0.005  ##<< numeric fraction(0-1): Variance threshold below which eigentriples are
                            ## treated as "noise" and will not be included in the groups. The actual
                            ## threshold is calculated by multiplying the total variance of
                            ## the time series with this fraction.
, grouping = c('grouping.auto', 'groupSSANearestNeighbour')[1]##<< character string: Method to use for grouping
                            ##  the individual SSA eigentriples. 'grouping.auto' uses the function
                            ## of that name in package Rssa, 'groupSSANearestNeighbour' employs a rather crude scheme
                            ## based on finding pairs (or larger groups) in an euclidian
                            ## distance matrix of the reconstructions of all extracted SSA eigentriples.
                            ## See ?grouping.auto or ?groupSSANearestNeighbour for details.
, groupingMethod = 'wcor'   ## character string: Method to use for grouping with the grouping.auto
                            ## function.
, repeat.extr = rep(1, times = length(borders.wl))##<< integer value/vector: How often to repeat the
                            ## extraction. If the respective value is > 1 than the result of the extraction
                            ## is again subject to spectral decomposition/filtering for n times and only
                            ## the (filtered) result is treated as the actual band (see details). 
, recstr.type = 'subtraction' ##<< string: How to determine the high frequency residuals.
                            ## If == 'subtraction', the high frequency signal is computed by subtracting
                            ## all other signals from the original series (see details). For all other
                            ## values only extracted eigentriples with high frequencies are grouped
                            ## in this band.
, pad.series = c(0, 0)      ##<< integer vector (length 2): Length of the part of the series to use for
                            ## padding at the start (first value) and at the end of the series. Values
                            ## of zero (default) cause no padding.
, SSA.methods = c('nutrlan', 'propack', 'eigen', 'svd')  ##<< character vector: Methods to use for the
                            ## SSA computation. First the first method is tried, when convergence
                            ## fails, the second is used and so on. See the help of ssa() in
                            ## package Rssa for details on the methods. It is preferable to use more
                            ## than one method as some methods (especially nutrlan) in some cases do not
                            ## converge. The last two methods are relatively slow.
, center.series = TRUE      ##<< logical: Whether to center the series around zero prior to the computation.
                            ## The (subtracted) mean will be added to the long term 'trend' component,
                            ## e.g. to the step containing an Inf value in  borders.wl. Not centering
                            ## of the series may cause erroneous trend extraction in some cases!
, call.freq = quote(calcFrequency(series.t)) ##<< 'quoted' function call : call to a function to compute
                            ## the frequency of the 'major' oscillation present in some time series.
                            ## This is used to compute the frequency of the (grouped) SSA eigentriples.
                            ## See the help for 'calcFrequency' for details of the default mechanism.
, n.steps = switch(class(borders.wl),list=length(borders.wl), dim(borders.wl)[2]) ##<< integer: 
                            ##  Amount of steps in the process. This argument is only kept
                            ## for backwards compatibility. Do not supply or change any values!
, plot.spectra = FALSE      ##<< logical: Whether to plot pseudo spectra for the different steps.
, second.axis = TRUE        ##<< logical: Whether to plot a second axis with period units
, open.plot = TRUE          ##<< logical: Whether to open a new plotting window for the plots. Set this to
                            ##   FALSE and open and set up a device prior to running the function to specify
                            ##   or change the device options.
, print.stat = TRUE         ##<< logical: whether to print status information during the calculations.
, xlim = c()                ##<< numeric vector: x-limits for the plotted spectra. If not supplied
                            ## it goes from 1 / n....0,5.
, debugging = FALSE         ##<< logical: if TRUE, workspaces are saved that can be used for debugging
                            ##   non convergence cases that do not cause R errors but may indicate
                            ##   a possible error in the settings, data or code.
, ...                       ##<< miscellaneous: further arguments passed to the plotting routines.
)
                              
##description<<
## This function decomposes (or filters) a time series into a set of orthogonal (i.e. additive)
## components with variance on different and distinct timescales (i.e. within different bands).
## It uses the fast and optimized  Singular Spectrum Analysis (SSA) method of Korobeynikov (2010).

##details<<
## Purpose
## The function is based on "singular spectrum analysis" (SSA) (Golyandina et al. (2001))
## based on the Rssa package (Korobeynikov (2010), see Golyandina et al. (2013) for a basic introduction).
##
## Definition of the period borders (borders.wl):
## borders.wl contains the borders of the different periodicity bands to
## extract. Units are sampling frequency of the series. borders.wl has to be a
## list with one element per desired decomposition step (see also 'stepwise
## extraction').  In the case of a single band to be extracted, this vector has
## to consist of two values (e.g. c(<lower border>,<upper.border>)).  In the
## case of the extraction of several bands per step, the upper border of each
## band is automatically the lower border of the next band. A value of c(0, 10,
## 100) would, hence, result in the extraction of two bands, one consisting of
## eigentriples with periods between 0 and 10 timesteps and one with periods
## between 10 and 100 timesteps. As a result the vector has to consist of one
## more value than desired groups. If the vector starts with 0 (and recstr.type
## = 'subtraction') all high frequency oscillations are combined this group. Use
## the value INF as the upper band for the long term trend (if desired)  to savely
## include all low frequency parts here.
##
## Window length (M):
## For optimal detection of spectral components the window length
## or embedding dimension M should be an integer multiple of the period of the
## oscillation that is to be extracted. See also ?ssa (here the same parameter is
## called L)
##
## Padding (pad.series):
## For padding, the series should start and end exactly at the start and end of a
## major oscillation (e.g. a yearly cycle). Additionally the length to use for padding
## should be a integer multiple of this period length. The padding is solved internally
## by adding the indicated part of the series at the start and at the end of the
## series. This padded series is only used internally and only the part of the
## series with original data is returned in the results.
##
## High frequency part (recstr.type):
## Two different ways to compute the high frequency part (e.g. the part with a
## frequency between 0 and the lowest border.wl value) are implemented. The
## standard way (if recstr.type == 'subtraction') is to sum all eigentriples with
## lower frequencies and subtract this sum from the original series to compute
## this as the high frequency residual. Otherwise (if (recstr.type !=
## 'subtraction')) only the eigentriples with such high frequencies that were
## actually extracted are used to build this spectral band.
##
## Stepwise extraction:
## In general the whole algorithm can be run stepwise. This means that first a
## certain spectral band is computed (with all possible harmonies etc. ...) and
## subtracted from the original series. Subsequently this process can be
## repeated with the residual as often as wanted. This allows for the adaptation
## of, for example, M, harmonics, ... to the particular oscillation to be
## extracted. Additionally it often this often leads to a clearer signal
## separation. To implement several steps, each of M, n.comp and harmonics needs
## to be a vector of the length of the amount of steps. For each step the
## corresponding element of this vector is used. borders.wl needs to contain one
## list entry per step which needs to be a vector containing all borders of the
## bands to extract in this step (see 'Definition of the period borders').
##
## Repeated extraction (repeat.extr):
## Especially for the trend it may be advisable to repeat the filtering step for
## this particular band. In this case the result of the first filtering
## (e.g. the sum of all eigentriples within this band) is filtered again n times
## with the same period borders. Finally only the final filtered components are
## summed and treated as the actual spectral band. In many cases this is helpful
## to exclude high frequency parts from the trend or low frequency components.
## It is only possible to repeat the extraction for steps where single bands are
## extracted.
## 
## Visualize results (plot.spectra):
## In the case that diagnostic plots should be plotted (plot.spectra == TRUE)
## one (or more) pseudospectra are plotted. Each point in these represents one
## group of eigentriples determined. For each step a separate plot is
## produced. Colored regions represent the specified spectral bands. grey lines in
## the background represent a Fourier Spectrum of the original series.

##references<<
## Golyandina, N.; Nekrutkin, V.; Nekrutkin, V. & Zhigljavsky, A. (2001),
## Analysis of time series structure: SSA and related techniques, CRC Press    
## Korobeynikov, A. (2010), Computation- and space-efficient implementation of SSA.
## Statistics and Its Interface, Vol. 3, No. 3, Pp. 257-268
## Golyandina, N. & Korobeynikov, A. (2013), 'Basic Singular Spectrum Analysis
## and forecasting with R', Computational Statistics & Data Analysis.

                              
##keyword<<
## SSA, time series, spectral analysis, singular spectrum analysis, spectral decomposition, filter

##seealso<<
##\code{\link{ssa}},\code{\link{calcFrequency}}
{ 
  ## transpose  borders.wl for backwards compatibility
  if (class(borders.wl) == 'data.frame'){
    borders.wl.new = list()
    for (b in 1: dim(borders.wl)[2])
      borders.wl.new[[b]] <- borders.wl[!is.na(borders.wl[,b]), b]
    borders.wl = borders.wl.new
  }

  ## check for NA margins at start and end and clean them
  na.margins <- c()
  if (sum(is.na(series)) > 0) {
    if (isSeriesContinuous(series)) {
      na.margins <- which(is.na(series))
      series     <- series[-na.margins]
    } else {
      stop('Series contains gaps. Remove them before decomposition!')
    }
  }
  
  ## check input and give warnings if necessary
  if (sum(diff(series) == 0) == length(series) - 1) {
    if (print.stat)
      warning('Series is apparently constant. Spectral decomposition makes no sense!')
  }
  if (missing(series) | missing(borders.wl))
    stop('Time series and/or vector with borders for spectral bands needs to be supplied!')
  if (sum(length(series)<(M-1))>0)
    stop('Window length exceeds series length!')
  if (sum(M / n.comp<2)>0 && is.element('nutrlan', SSA.methods)) {
    print(paste('n.comp exceeds M / 2! As using nutrlan in such cases is not appropriate,  ',
            'only the remaining methods are used!', sep = ''))
    SSA.methods <- SSA.methods[-match('nutrlan', SSA.methods)]
    if (length(SSA.methods) == 0)
      stop('No SSA methods left. Specify other methods or adjust n.comp to M!')
  }
  if (sum(n.comp>M)>0)
    stop('Amount of eigentriples to extract exceeds window length!')
  if (sum(length(series)<(n.comp + M))>0)
    stop('M + n.comp exceeds series length!')
  if (sum(!(c(M, n.comp, harmonics)-trunc(c(M, n.comp, harmonics)) == 0))>0 | sum(c(M, n.comp, harmonics)<0)>0 )
    stop('Invalid parameters! M, n.comp and harmonics need to be positive integers!')
  if (sum(!diff(c(length(borders.wl), length(harmonics), length(n.comp), length(M))) == 0)>0)
    stop('M, harmonics and n.comp have to be all vectors of the same length as dim(borders.wl)[1]!')
  if(!is.numeric(pad.series) || !length(pad.series) == 2 || sum(pad.series<0)>0)
    stop('pad.series has to be a positive integer vector of length 2!')
  if(sum(!rapply(borders.wl, function(x)sum(diff(x)<0) == 0))>0)
    stop('Each element of borders.wl has to be a monotonously rising (!) vector!')
  if(sum(pad.series>length(series))>0)
    stop('pad.series must not exceed the length of the timeries!')
  test.vec <- c(harmonics, unlist(borders.wl), M, n.comp, harmonics, pad.series)
  test.vec <- test.vec[is.finite(test.vec)]
  if(sum(!(abs(test.vec-round(test.vec))<.Machine$double.eps^0.5))>0)
    stop('All values of harmonics, borders.wl, M, n.comp and pad.series have to be integer values!')
  if (sum(unlist(lapply(borders.wl, length))>2 & harmonics>0)>0 )
    stop(paste('Do not try to extract harmonics during steps in which',
            'you extract more than one period band.',
            'This is possible theoretically but not implemented.', sep = ''))
  if (!length(borders.wl) == length(repeat.extr) | !class(repeat.extr) == 'numeric')
    stop('Lengths of borders.wl and repeat.extr do not correspond!')
  if (sum(sapply(borders.wl[repeat.extr>1], length)>2))
    stop(paste('If you want to perform multiple, subsequent repeated extractions (repeat.extr>1), ',
            'extract only one band during this step!', sep = ''))
  if (class(series) != 'numeric')
    stop('Only numeric objects are allowed for series (e.g. no ts, zoo, etc ... objects )!')

  if (print.stat)
    printStatus('Spectral decomposition: Preparing calculations.')

  ## pad series
  if (sum(pad.series>0)>0) {
    ind.padded     <- c(0:pad.series[1], seq(  from = (pad.series[1] + length(series) + 1),
            length.out = pad.series[2]  ) )
    ind.padded     <- ind.padded[!ind.padded == 0]
    series         <- c(series[0:pad.series[1]], series,
        series[-(1:(length(series)-pad.series[2]))])
  }

  ## center series
  if (center.series){
    mean.series    <- mean(series)
    series.work    <- series - mean.series
  } else {
    series.work    <- series
  }

  ## prepare stuff
  n.steps <- length(borders.wl)
  if (open.plot & plot.spectra){
    layout(matrix(1:n.steps, ncol = 1))
    par(tcl = 0.2, mgp = c(1, 0, 0), mar = c(0, 0, 0, 0), oma = c(2, 2, 2, 2), ps = 10, cex = 1)
  }
  n.datapts      <- length(series)
  var.total      <- var(series)
  n.series.res   <- length(unlist(borders.wl)) - length(borders.wl)
  series.perstep <- sapply(borders.wl, length) - 1
  series.results <- matrix(NA, ncol = n.datapts, nrow = n.series.res)
  borders.results<- matrix(NA, ncol = 2, nrow = n.series.res)
  borders.results[, 1] <- unlist(sapply(borders.wl, function(x)x[-length(x)]))
  borders.results[, 2] <- unlist(sapply(borders.wl, function(x)x[-1]))
  row.fill       <- 1
  ylims.SSA      <- c()
  ylims.FFT      <- c()
  if (length(xlim) == 0)
    xlim = c(1 / n.datapts, 0.5)
  break.h.loop = FALSE
  if (!(sum(diff(series) == 0) == length(series)-1)) {   ##if series not constant

    ## run calculation for all steps
    for (h in 1:n.steps) {
      if (print.stat)
        printStatus(paste('Spectral decomposition: Starting step ', h, ' of ',
                n.steps, '.', sep = ''))
      n.filtered.series        <- series.perstep[h]
      rows.results.step.single <- (sum(series.perstep[0:(h-1)])) + c(1, series.perstep[h])
      rows.results.step        <- rows.results.step.single[1]:rows.results.step.single[2]
      if (print.stat)
        printStatus('Spectral decomposition: Running SSA.')
      for (j in 1:n.filtered.series) {
        if (sum(series.work == 0) == length(series)) {
          series.results[(-1:rows.results.step[j]-1), ] <- 0
          break.h.loop  =  TRUE
          break
        }
        for (n in 1: repeat.extr[[h]]) {
          if (n == 1) {
            series.input = series.work
          } else {                                          
            series.input = series.results[rows.results.step[j], ]
            if (sum(series.input == 0) == length(series))
              break
          }
          
          ## run SSA
          res.run          <- .calcSSAAllMethods(series.in = series.input, M = M[h], n.comp = n.comp[h],
                                                 SSA.methods = SSA.methods, kind = '1d-ssa',
                                                 GroupEigTrpls = grouping, groupingMethod = groupingMethod,
                                                 debugging = debugging)
          groups.ssa       <- res.run$ssa.groups.t
          n.pairs          <- length(groups.ssa)
          r.final          <- res.run$recstr.res
          r.matrix.final   <- matrix(unlist(r.final), nrow = n.pairs, byrow = TRUE)
          method.used      <- res.run$method.used

          ## determine frequencies/periods of SSA components
          if (print.stat)
            printStatus('Spectral decomposition: Determining frequencies.')
          frequency        <- numeric(length = n.pairs)
          variance         <- numeric(length = n.pairs)
          for (i in 1:n.pairs) {
            series.t     <- r.matrix.final[i, ]
            frequency[i] <- eval(call.freq)
            variance[i]  <- var(series.t)
          }
          period           <- 1 / frequency

          ## determine to which bandwith group each component belongs
          if (print.stat)
            printStatus('Spectral decomposition: Summing groups in wavlength bands.')
          pairs.valid.step        <- c()
          pairs.inband.step       <- c()

          ## determine eigentriples in spectral band
          pairs.valid.band    <- c()
          pairs.inband.band   <- which(period > borders.wl[[h]][j] & period < borders.wl[[h]][j + 1])
          freq.main           <- frequency[pairs.inband.band[which.max(variance[pairs.inband.band])]]

          ## determine eigentriples in harmonic bands
          if (harmonics[h] > 0) {
            if (length(pairs.inband.band) > 0) {
              for (k in 1:harmonics[h]) {
                frequency.harmonic      <- freq.main * (k + 1)
                pairs.inband.t          <- which(frequency < ((1 + tolerance.harmonics) *
                                                              frequency.harmonic  ) &
                                                 frequency > ((1-tolerance.harmonics) *
                                                              frequency.harmonic  ))
                pairs.inband.band       <- union(pairs.inband.band, pairs.inband.t)
              }
            }
            
            ## remove pairs with variance lower than treshhold
            pairs.valid.band   <- pairs.inband.band[which(variance[pairs.inband.band] >  (var.total * var.thresh.ratio))]
          } else if(harmonics[h] == 0) {
            pairs.valid.band   <- pairs.inband.band
          }

          ## merge pairs for bands into those per step
          pairs.valid.step                       <- union(pairs.valid.step, pairs.valid.band)
          pairs.inband.step                      <- union(pairs.inband.step, pairs.inband.band)

          ## sum respective bands and fill to results matrix
          if (length(pairs.valid.band) > 0) {                                  ## if valid eigentriples in band
            series.results[rows.results.step[j], ] <- colSums(matrix(r.matrix.final[pairs.valid.band, ],
                    ncol = n.datapts))
          } else {
            series.results[rows.results.step[j], ] <- 0
            break
          }
        }
        
        

        ## subtract identified eigentriples from series (to use residuals for next iteration)
        series.work <-  series.work - colSums(matrix(series.results[rows.results.step[j], ], ncol = n.datapts))
    }
      ## plot pseudospectrum
      if (plot.spectra) {
          if (print.stat)
              printStatus('Spectral decomposition: Plotting pseudospectrum.')

          ## calculate and plot Fourier Spectrum
          spec.fourier <- spectrum(series.work, plot = FALSE, span = 3)
          plot(spec.fourier$freq, spec.fourier$spec, type = 'l', log = 'xy', ylab = '',
               yaxt = 'n', xaxt = 'n', xlab = '', col = 'gray', xlim = xlim, ylim = ylims.FFT)#, ...)
          if (h == 1)
              ylims.FFT <- 10^par()$usr[3:4]
          axis(4)
          par(new = TRUE)

          ## plot SSA pseudospectrum
          plot(frequency, variance, pch = 1, cex = 0.8, ylab = '', xlab = '', log = 'xy',
               xlim = xlim, xaxt = 'n', ylim = ylims.SSA, ...)
          if (h == 1)
              ylims.SSA <- 10^par()$usr[3:4]

          ## add coloured polygons to mark band widths
          y.extremes <- 10^(par()$usr[c(3, 4)])
          x.extremes <- 10^(par()$usr[c(1, 2)])
          colors <-rainbow(n.filtered.series, alpha = 0.4)
          for (l in 1:n.filtered.series) {
              for (m in 0:(harmonics[h])) {
                  if (m>0) {
                      x.poly <- c(freq.main * (m + 1) * (c(-1, 1, 1, -1) * tolerance.harmonics + 1))
                  } else {
                      x.poly <- 1 / c(borders.wl[[h]][l:(l + 1)], rev(borders.wl[[h]][l:(l + 1)]))
                  }
                  x.poly[x.poly == 0] <- x.extremes[1]
                  x.poly[!is.finite(x.poly)]<-x.extremes[2]
                  y.poly <- rep(y.extremes, each = 2)
                  polygon(x.poly, y.poly, col = colors[l])
              }
          }
          abline(v = 1 / unique(unlist(borders.wl)), lty = 2)

          ## add points and threshholds for better visibility
          points(frequency, variance, pch = 1, cex = 0.8)
          points(frequency[pairs.inband.step], variance[pairs.inband.step], cex = 1, pch = 20, col = 'blue')
          points(frequency[pairs.valid.step], variance[pairs.valid.step], cex = 1, pch = 20, col = 'green')
          abline(h = var.total * var.thresh.ratio, lty = 2)

          ## add x axis at top and bottom of plot
          if (h == 1 &second.axis)
              plotAdditionalAxis(side = 3, label = '', trans.fun = function(x)1 / x)
          if(h == n.steps)
              axis(1)
          if (open.plot) {
              mtext(text = 'log(variance SSA eig.trp.)', outer = TRUE, side = 2, line = 1)
              mtext(text = 'raw periodogram', outer = TRUE, side = 4, line = 1)
              mtext(text = 'wavelength / period [timesteps]', side = 3, line = 1, outer = TRUE)
              mtext(side = 1, 'frequency [1 / timestep]', line = 1, outer = TRUE)
              text(10^par()$usr[1], 10^par()$usr[4], adj = c(0, 1), labels = paste('step ', h, sep = ''), cex = 1.5)
              
              ## add legend
              if (h == 1) {
                  legend('topright', legend = c('all SSA eigtrp.', 'SSA eigtrp. in band',
                                         'valid SSA eigtrp. in band', 'Fourier periodogram'),
                         col = c('black', 'blue', 'green', 'gray'), pch = c(20, 20, 20, NA), lty = c(NA, NA, NA, 1),
                         merge = TRUE, cex = 0.8, ncol = 2, x.intersp = 0, y.intersp = 0.6, bg = 'white')
              }
          }
      }
      if (break.h.loop)
          break
  }

    ## return results
    if (print.stat)
      printStatus('Spectral decomposition: Process finished!')

    ## add subtracted mean to trend component
    if (center.series) {
      row.trend                   <- which(apply(borders.results, MARGIN = 1, function(x)is.element(Inf, x)))
      series.results[row.trend, ] <- series.results[row.trend, ] + mean.series
    }

    ## calculate high frequency residual
    if (recstr.type  == 'subtraction' && sum(borders.wl[[length(borders.wl)]] == 0) == 1) {
      row.high.freq <- which(apply(borders.results, MARGIN = 1, function(x)is.element(0, x)))
      series.results[row.high.freq, ]<- series  - colSums(matrix(series.results[-row.high.freq, ], ncol = n.datapts))
    }
  } else {                                                                        ## if series constant
    row.trend      <- which(apply(borders.results, MARGIN = 1, function(x)is.element(Inf, x)))
    series.results <- matrix(0,  ncol = length(series),  nrow = n.series.res)
    series.results[row.trend, ] <- mean(series,  na.rm = TRUE)
    g              <- 1
  }

  ## remove padding
  if (sum(pad.series>0)>0)
    series.results <- series.results[, -ind.padded]

  ## add outer na margins
  if (length(na.margins) > 0) {
    series.final   <- array(NA, dim = dim(series.results) + c(0, length(na.margins)))
    series.final[, - na.margins] <- series.results
  } else {
    series.final <- series.results
  }
  
  ## save settings
  settings           <- list(borders.wl = borders.wl, M = M, n.comp = n.comp, harmonics = harmonics,
                             tolerance.harmonics = tolerance.harmonics, var.thresh.ratio = var.thresh.ratio,
                             grouping = grouping, recstr.type = recstr.type, pad.series = pad.series,
                             SSA.method = method.used, repeat.extr = repeat.extr)
  settings$char.string <- paste(paste(names(settings), sapply(settings, function(x)paste(x, collapse = ',  '))
                                      , sep = ': '), collapse = '; ')
  ##value<< list with components
  list(dec.series = series.final ##<< numeric matrix: decomposed timeseries. Each row of this
                                 ## matrix contains one spectrally filtered component. There
                                 ## are as many rows as period borders (borders.wl) were
                                 ## defined. 
      , borders = borders.results##<< numeric matrix: The lower (first column) and upper
                                 ## (second column) period borders of each filtered component.
                                 ## The rows here correspond to the rows in "dec.series".
      , settings = settings      ##<< list: Settings used to perform the calculation.
  ##end<<
  )
}, ex = function(){
  #create series consisting of two oscillations and noise
  series.ex <- sin(2 * pi * 1:1000 / 100) +  0.7 * sin(2 * pi * 1:1000 / 10)  +
    rnorm(n = 1000, sd = 0.4)

  #prepare graphics
  layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8), ncol = 2))
  par(tcl = 0.2, mgp = c(2, 0, 0), mar = c(0, 4, 0, 0), oma = c(2, 0, 2, 0),
      ps = 10, cex = 1)
  plot.new()

  #perform decomposition
  data.decomposed <- filterTSeriesSSA(series = series.ex,
      borders.wl = list(a = c(8, 12), b = c(80, 120)
          , c = c(0, 10, 100, Inf)),
      M = c(30, 200, 100),
      n.comp = c(10, 20, 20),
      harmonics = c(1, 0, 0),
      plot.spectra = TRUE, open.plot = FALSE)

  #plot series and spectral parts
  plot(series.ex)
  plot(data.decomposed$dec.series[1, ], ylab = '')
  plot(data.decomposed$dec.series[2, ], ylab = '')
  plot(colSums(data.decomposed$dec.series[-c(1:2), ]), ylab = '')
  mtext(side = 2, outer = TRUE, at = -(1 / 8) + ((4:1) * (1 / 4)),
      c('orig.series', '1.step', '2.step', '3.step'), las = 3, cex = 1.5, line = -1)
  mtext(side = 3, outer = TRUE, at = -(1 / 4) + ((1:2) * (1 / 2)),
      c('pseudospectra', 'identified components'), las = 1, cex = 1.5, line = 1)

})
