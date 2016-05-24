plotScatterAcd <- function(fitModel, x = "muHats", y = "residuals", xlag = 0, ylag = 0,
                           colour = NULL, xlim = NULL, ylim = NULL, alpha = 1/10,
                           smoothMethod = "auto"){  
  
  x <- match.arg(x, c("muHats", "residuals", "durations", "adjDur", "dayTime", "time", "index"))
  y <- match.arg(y, c("muHats", "residuals", "durations", "adjDur", "dayTime", "time", "index"))
  if(length(colour) != 0) colour <- match.arg(colour, c("muHats", "residuals", "durations", "adjDur", "dayTime", "time", "index"))
  contTime = TRUE
  
  xData <- switch(x,
                  muHats = fitModel$muHats,
                  residuals = fitModel$residuals,
                  durations = fitModel$durations$durations,
                  adjDur = fitModel$durations$adjDur,
                  dayTime = fitModel$durations$time$min / 60 + fitModel$durations$time$hour,
                  time = {if(contTime) fitModel$durations$time
                          else fitModel$durations$time$yday * (60*8 + 25) + fitModel$durations$time$min + fitModel$durations$time$hour * 60},
                  index = 1:fitModel$N)  
  yData <- switch(y,
                  muHats = fitModel$muHats,
                  residuals = fitModel$residuals,
                  durations = fitModel$durations$durations,
                  adjDur = fitModel$durations$adjDur,
                  dayTime = fitModel$durations$time$min / 60 + fitModel$durations$time$hour,
                  time = {if(contTime) fitModel$durations$time
                          else fitModel$durations$time$yday * (60*8 + 25) + fitModel$durations$time$min + fitModel$durations$time$hour * 60},
                  index = 1:fitModel$N)
  if(length(colour) != 0){
    colourData <- switch(colour,
                         muHats = fitModel$muHats,
                         residuals = fitModel$residuals,
                         durations = fitModel$durations$durations,
                         adjDur = fitModel$durations$adjDur,
                         dayTime = fitModel$durations$time$min / 60 + fitModel$durations$time$hour,
                         time = {if(contTime) fitModel$durations$time
                                 else fitModel$durations$time$yday * (60*8 + 25) + fitModel$durations$time$min + fitModel$durations$time$hour * 60},
                         index = 1:fitModel$N,
                         NULL = NULL)
    colourData <- colourData[(1+ylag):length(colourData)]
  }
  
  yData <- yData[(1+xlag):(length(yData)-ylag)]
  xData <- xData[(1+ylag):(length(xData)-xlag)]
  
  if(ylag != 0) y <- paste("lagged ", y, " (i-", ylag, ")", sep = "")
  if(xlag != 0) x <- paste("lagged ", x, " (i-", xlag, ")", sep = "")

  if(length(colour) == 0){ 
    g <- ggplot(data.frame(x = xData, y = yData), aes(x = x, y = y))
  } else{
    g <- ggplot(data.frame(x = xData, y = yData, colour = colourData), aes(x = x, y = y, colour = colour)) + scale_colour_continuous(name = colour)
  }    
  g <- g + geom_point(alpha = alpha) + geom_smooth(colour="red", size=1.5, fill = "blue", alpha = .2, method = smoothMethod)
  if(x == "muHats" && y == "residuals") g <- g + scale_y_continuous(breaks = seq(1, max(yData), 1)) #+ geom_hline(yintercept = 1, colour = "red") 
  if(length(xlim) != 0) g <- g + xlim(xlim)
  if(length(ylim) != 0 ) g <- g + ylim(ylim) 
  g + ylab(y) + xlab(x) 
}