DrawChromatogram <- function(time, intensity, range = list(start, stop), color = "blue", 
  xlab = "retention time", ylab = "intensity", ylim = c(0, max(intensity) * 1.1), las = 1, ...) { 

  if(length(range$start) != length(range$stop)) stop("start and stop vectors in the range list must be equal")
  if(length(color) != 1 & length(color) != length(range$start))
    stop("color vector length must be 1 or equal to the vectors in the range list")

  plot(time, intensity, type = "l", ylim = ylim, 
    xlab = xlab, ylab = ylab, las = las, ...)

  # calculate the retention time, area, and apex intensity of each peak 
  retentionTime <- vector(mode = "numeric", length = length(range$start))
  peakArea <- vector(mode = "numeric", length = length(range$start))
  apexIntensity <- vector(mode = "numeric", length = length(range$start))
  for(i in 1:length(range$start)) {
    peakTime <- time[time >= range$start[i] & time <= range$stop[i]]
    peakIntensity  <- intensity[time >= range$start[i] & time <= range$stop[i]]
    if(length(color) == 1) peakColor <- color
    else peakColor <- color[i]
    polygon(peakTime, peakIntensity, col = peakColor)

    # calculate polygon (peak) area
    n <- length(peakTime)
    x <- vector(mode = "numeric", length = n)
    for(j in 1:n) {
      k <- (j %% n) + 1
      x[j] <- peakTime[j] * peakIntensity[k] - peakTime[k] * peakIntensity[j]
    }

    peakArea[i] <- abs(sum(x) / 2)

    retentionTime[i] <- peakTime[peakIntensity == max(peakIntensity)]
    apexIntensity[i] <- max(peakIntensity)

  }

  return(data.frame(retentionTime, peakArea, apexIntensity))

}

