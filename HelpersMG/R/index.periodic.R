#' index.periodic estimate indices in periodic timeseries based on anchored minimum and maximum
#' @title Estimate indices in periodic timeseries based on anchored minimum and maximum
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @return A data.frame with a column time and a column index
#' @param minmax A data.frame returned by minmax.periodic
#' @param time	The time at which produced the estimate
#' @param replicates Number of replicates to estimate SD
#' @param progressbar Does a progression bar must be shown
#' @family Periodic patterns of indices
#' @description Estimate indices in periodic timeseries based on anchored minimum and maximum.\cr
#' The data.frame minmax can be generated manually. It should have three columns (time, index, SD), 
#' with all the successive minimum and maximum indices.\cr
#' It can be used with sun.info() to get the time of minimum and maximum air temperature or
#' with getTide() to reconstruct the sea level.
#' @examples
#' \dontrun{
#' # Generate a timeserie of time
#' time.obs <- NULL
#' for (i in 0:9) time.obs <- c(time.obs, c(0, 6, 12, 18)+i*24)
#' # For these time, generate a timeseries of temperatures
#' temp.obs <- rep(NA, length(time.obs))
#' temp.obs[3+(0:9)*4] <- rnorm(10, 25, 3)
#' temp.obs[1+(0:9)*4] <- rnorm(10, 10, 3)
#' for (i in 1:(length(time.obs)-1)) 
#'   if (is.na(temp.obs[i])) 
#'   temp.obs[i] <- mean(c(temp.obs[i-1], temp.obs[i+1]))
#'   if (is.na(temp.obs[length(time.obs)])) 
#'   temp.obs[length(time.obs)] <- temp.obs[length(time.obs)-1]/2
#' observed <- data.frame(time=time.obs, temperature=temp.obs)
#' # Search for the minimum and maximum values
#' r <- minmax.periodic(time.minmax.daily=c(Min=2, Max=15), 
#' observed=observed, period=24, colname.index="temperature")
#' 
#' # Estimate all the temperatures for these values
#' t <- index.periodic(minmax=r)
#' 
#' plot_errbar(x=t[,"time"], y=t[,"index"],
#' errbar.y=ifelse(is.na(t[,"sd"]), 0, 2*t[,"sd"]),
#' type="l", las=1, bty="n", errbar.y.polygon = TRUE, 
#' xlab="hours", ylab="Temperatures", ylim=c(0, 35), 
#' errbar.y.polygon.list = list(col="grey"))
#' 
#' plot_add(x=t[,"time"], y=t[,"index"], type="l")
#' 
#' plot_add(observed$time, observed$temperature, pch=19, cex=0.5)
#' }
#' @export

index.periodic <- function(minmax, time=NULL, 
                                 replicates=100, progressbar=FALSE) {
  if (class(minmax) !="data.frame") {
    warning("minmax parameter must be a data.frame with at least the three columns: 'time', 'index' and 'SD'")
    return()
  }
  if (sum(match(colnames(minmax), c("index", "time", "SD")), na.rm=TRUE) !=6) {
    warning("minmax parameter must be a data.frame with at least the three columns: 'time', 'index' and 'SD'")
    return()
  }
  if (is.null(time)) {
    time <- seq(from=minmax$time[1], to=tail(minmax$time, n=1), by=1)
  }
  if (progressbar) pb<-txtProgressBar(min=1, max=length(time), style=3)
  
  dt <- data.frame(time=time, index=rep(NA, length(time)), sd=rep(NA, length(time)))
  for (i in 1:length(time)) {
    
    if (progressbar) setTxtProgressBar(pb, i)
    
    
   g <- NULL
    tec <- time[i]
    dif <- minmax[,"time"]-tec
    change.sign <- which(sign(dif[1:(length(dif)-1)])*sign(dif[2:length(dif)])==-1 | 
    	sign(dif[1:(length(dif)-1)])*sign(dif[2:length(dif)])==0)[1]
    
    
    if (!identical(change.sign, integer(0)) & !is.na(change.sign)) {
  x0 <- minmax[change.sign, "time"]
  x1 <- minmax[change.sign+1, 'time']
  y0 <- minmax[change.sign, 'index']
  y1 <- minmax[change.sign+1, 'index']
  a <- pi/(x1-x0)
  b <- -x0*a
  p <- (y0-y1)/2
  q <- y0-p

  x <- tec
  y <- cos(x*a+b)*p+q
  dt[i, "index"] <- y
  
  
    y0 <- rnorm(replicates, mean=minmax[change.sign, 'index'], 
                sd=ifelse(is.na(minmax[change.sign, 'SD']), 0, minmax[change.sign, 'SD']))
  y1 <- rnorm(replicates, mean=minmax[change.sign+1, 'index'], 
              sd=ifelse(is.na(minmax[change.sign+1, 'SD']), 0, minmax[change.sign+1, 'SD']))

  for (j in 1:replicates) {
  
  a <- pi/(x1-x0)
  b <- -x0*a
  p <- (y0[j]-y1[j])/2
  q <- y0[j]-p

  x <- tec
  y <- cos(x*a+b)*p+q

	g <- c(g, y)
	}
	
	dt[i, "sd"] <- sd(g)

    
    }
  }
  return(dt)
  
}
