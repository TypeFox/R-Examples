#' minmax.periodic search for minimum and maximum indices (temperatures or levels) in periodic timeseries
#' @title Search for minimum and maximum indices in periodic timeseries
#' @author Marc Girondot
#' @return A data.frame with a column time, a column index and a column SD
#' @param time.minmax.daily A named vector with Min and Max being the time in the day with minimum and maximum indices (temperature or level)
#' @param time.minmax A named vector daily with time in the day at which minimum and maximum indices are observed
#' @param observed	A dataframe with at least two columns: time and temperatures. A third column SD can indicate the know error in index
#' @param period The unit of day period (24 for hours, 24*60 for minutes)
#' @param progressbar Tell if a progression bar must be shown
#' @param colname.time The name of the column for time in observed
#' @param colname.index The name of the column for indices in observed
#' @param colname.SD The name of the column for SD in observed
#' @param plot If TRUE, show a plot with the different estimates
#' @family Periodic patterns of indices
#' @description Search for minimum and maximum for periodic timeseries when only intermediate values are known.\cr
#' For each couple of value with an increasing or decreasing segment of 
#' the sinusoid function, it is possible to estimate a minimum and maximum 
#' values using analytical algebra.\cr
#' Then the average and standard deviations of all minima and maxima are evaluated.\cr
#' It should be noted that any extremum can be estimated at least twice, one by
#' increasing segment and one by decreasing segment. Both are used here to produce SD.\cr
#' \code{time.minmax.daily} should be used when the time at which maximum and minimum indices are regular and
#' \code{time.minmax} permits to define this time day by day.
#' @examples
#' \dontrun{
#' library("HelpersMG")
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

minmax.periodic <- function(time.minmax.daily=NULL, time.minmax=NULL, progressbar=FALSE, 
                            observed=stop("data.frame with observed indices"), 
                            period=24, 
                            colname.time="time", colname.index="index", colname.SD="SD", 
                            plot=FALSE) {

if (is.null(time.minmax) & is.null(time.minmax.daily)) {
  stop("Or time.minmax or time.minmax.daily must be provided")
}
if (!is.null(time.minmax) & !is.null(time.minmax.daily)) {
  stop("time.minmax and time.minmax.daily cannot be both provided")
}


observed <- observed[order(observed[,colname.time]), ]

temp.obs <- observed[,colname.index]
time.obs <- observed[,colname.time]
sd.obs <- if(any(colnames(observed)==colname.SD)) observed[,colname.SD] else NULL

# SI j'ai un sd sur les observations, je prends + ou - 1.96 SD
if (!is.null(sd.obs)) {
  temp.obs <- c(temp.obs-1.96*sd.obs, temp.obs+1.96*sd.obs)
  time.obs <- c(time.obs, time.obs)
}


# J'ai juste le time.minmax.daily, je dois reconstruire le time.minmax
if (is.null(time.minmax)) {
  # nombre de périodes
  maxday <- max(time.obs)/period+2
  time.minmax <- NULL
  for (i in 1:maxday) time.minmax <- c(time.minmax, 
                                     Min=unname(period*(i-1)+time.minmax.daily["Min"]), 
                                     Max=unname(period*(i-1)+time.minmax.daily["Max"]))
}

# Je ne garde que ceux dont j'ai besoin
if (any(time.minmax>max(time.obs))) time.minmax <- time.minmax[1:which(time.minmax>max(time.obs))[1]]

if (is.null(names(time.minmax)) | any(names(time.minmax)=="")) {
  stop("time.minmax must be named vector using names Min or Max")
}

# je classe time.minmax
time.minmax <- sort(time.minmax)

res <- data.frame(time=numeric(), temperature=numeric(), extrema=character(), stringsAsFactors=FALSE)

if (progressbar) pb<-txtProgressBar(min=1, max=(length(time.minmax)-1), style=3)

if (plot) {
  plot(x=time.obs, y=temp.obs, pch=4, bty="n", xlab="Time", ylab="Index")
  for (i in seq_along(time.minmax)) segments(x0=time.minmax[i], x1=time.minmax[i], 
                                            y0=ScalePreviousPlot()$ylim[1], 
                                            y1=ScalePreviousPlot()$ylim[2], lty=2)
}

for (i in 1:(length(time.minmax)-1)) {
  if (progressbar) setTxtProgressBar(pb, i)
  
  # Je cherche la valeur y0 et y1 a x0 et x1
  x0 <- time.minmax[i]
  x1 <- time.minmax[i+1]
  sen0 <- names(time.minmax)[i]
  sen1 <- names(time.minmax)[i+1]
  
  # Si je n'ai pas une observation pour exactement ces temps
  # Je ne comprends pas pourquoi je fais ce test
  if (all(time.obs!=x0) | all(time.obs!=x1)) {
    #  print(paste(1, i,length(res)))
    # Je prends les temps avec des observations entre
    lim <- (time.obs>=x0 & time.obs<=x1)
    wlim <- which(lim)
    if (length(wlim)>=2) {
      # Toutes ls combinaisons possibles
      couple <- t(combn(wlim, 2))
      for (j in 1:nrow(couple)) {
        x <- time.obs[couple[j, 1]]
        xp <- time.obs[couple[j, 2]]
        if (x!=xp) {
        y <- temp.obs[couple[j, 1]]
        yp <- temp.obs[couple[j, 2]]
        C <- cos(x*pi/(x1-x0)-x0*pi/(x1-x0))
        Cp <- cos(xp*pi/(x1-x0)-x0*pi/(x1-x0))
        L <- 1/2*(1-C)
        K <- 1/2*(Cp+1) 
        y0 <- (yp/K -y/(2*K*L) +(Cp*y)/(2*K*L)) / 
          (1 -C/(4*K*L) -1/(4*K*L) +(Cp*C)/(4*K*L) +(Cp)/(4*K*L))
        y1 <- y/L -C*y0/(2*L) -y0/(2*L)
        
        if (plot) {
          a <- pi/(x1-x0)
          b <- -x0*a
          p <- (y0-y1)/2
          q <- y0-p
          
          x <- seq(from=x0, to=x1, length=15)
          y <- cos(x*a+b)*p+q
          
          plot_add(x, y, type="l")
        }
        
        if (all(time.obs!=x0)) res <- rbind(res, data.frame(time=x0, temperature=y0, extrema=sen0))
        if (all(time.obs!=x1)) res <- rbind(res, data.frame(time=x1, temperature=y1, extrema=sen1))
        }
      }
    }
  }
  #    print(paste(2, i,length(res)))
}

##### fin du premier tour
##### je fais le deuxième

mn <- aggregate(temperature ~ time, subset(res, res$extrema=="Min"), mean)
mx <- aggregate(temperature ~ time, subset(res, res$extrema=="Max"), mean)
mn2 <- aggregate(temperature ~ time, subset(res, res$extrema=="Min"), sd)
mx2 <- aggregate(temperature ~ time, subset(res, res$extrema=="Max"), sd)

mx <- cbind(mx, mx2[,2])
mx <- cbind(mx, extrema=rep("Max", nrow(mx)), stringsAsFactors=FALSE)
mn <- cbind(mn, mn2[,2])
mn <- cbind(mn, extrema=rep("Min", nrow(mn)), stringsAsFactors=FALSE)

names(mx) <- c("time", "index", "SD", "extrema")
names(mn) <- c("time", "index", "SD", "extrema")

mxmn <- rbind(mx, mn)

o <- order(mxmn$time)
mxmn2 <- mxmn[o,]

rownames(mxmn2) <- 1:dim(mxmn2)[1]

return(mxmn2)

}


