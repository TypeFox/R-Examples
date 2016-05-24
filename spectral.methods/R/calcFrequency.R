
calcFrequency = function(
##title<< Estimate the frequency of a periodic signal
        series                             ##<< numeric vector: input vector (time series)
        , plot.periodogram = FALSE         ##<< logical: whether to plot a periodogram
)
  ##description<<
  ##Function to estimate the "dominant" frequency of a periodic time series signal.
  ##details<<
  ##This function uses Fourier decomposition to determine the 'major' frequency of a time
  ##series. Technically this is the frequency of the Fourier component with the highest
  ##variance. The function is used by filterTSeriesSSA to determine the frequencies
  ##of the individual SSA components.
  ##seealso<<
  ##\code{\link{fft}},\code{\link{filterTSeriesSSA}}
{
    n.series       <- length(series)
    fft.results    <- fft(series) / n.series
    periodogram    <- abs(fft.results[2:(ceiling(n.series / 2) + 1)]) ^ 2
    frequencies    <- (1:(ceiling(n.series / 2))) / n.series
    ind.highest    <- order(periodogram, decreasing = TRUE)[1:2]
    frequency.main <- sum(frequencies[ind.highest] * periodogram[ind.highest]) /
                          sum(periodogram[ind.highest])

    if(plot.periodogram)
    {
        dev.new()
        plot(frequencies, periodogram, log='xy')
        abline(v = frequency.main, col = 'red', lty = 2)
        text(frequency.main, 10 ^ par()$usr[4], round(1 / frequency.main, digits = 2), col = 'red',adj = c(0, 1))
        xaxis.values <- axTicks(side = 1)
        side.new     <- 3
        axis(side = side.new, at = xaxis.values, labels = round(1 / xaxis.values, digits = 2))
        mtext(side = 3, 'period', line = 2)
    }

    ##value<< Frequency of the Fourier component with the highest variance [1/time steps]
    return(frequency.main)
}
