diurnalAdj <- function(dur, method = "cubicSpline", nodes = c(seq(600, 1105, 60), 1105), aggregation = "all", span = "cv", spar = 0, Q = 4, returnSplineFnc = FALSE){
  
  durations <- spline.x <- spline.y <- day <- x <- y <- time <- NULL
  
  if(!("POSIXlt" %in% class(dur$time))) dur$time <- as.POSIXlt(dur$time)
  
  #provides the possibility of entering truncated arguments:
  method <- match.arg(method, c("cubicSpline", "supsmu", "smoothSpline", "FFF"))
  aggregation <- match.arg(aggregation, c("weekdays", "all", "none"))
  if(nodes[length(nodes)] == nodes[length(nodes) - 1]) nodes <- nodes[-length(nodes)] #if the last two nodes are the same, the last is removed
  if(returnSplineFnc == TRUE) rtSpline <- list()
  
  if(method == "cubicSpline" | method == "smoothSpline"){
    
    if(method == "cubicSpline") splnFnc <- function(x, y) splines::interpSpline(x, y)
    else splnFnc <- function(x, y) stats::smooth.spline(x, y, all.knots = TRUE, spar = spar)
    
    timeInMinutes <- dur$time$hour * 60 + dur$time$min
    timeInterval <- numeric(nrow(dur))
    
    if (any(timeInMinutes < nodes[1] |
            timeInMinutes > nodes[length(nodes)]))
      stop(
        "\nAt least one of the durations occured outside of the nodes. \nThe smallest and largest nodes should be at the opening and closing time."
      )
    if (nodes[length(nodes) - 1] > max(timeInMinutes))
      warning(
        "no durations occured at the latest interval. Check if the 'node' argument is correctly specified"
      )
    if (nodes[1] < min(timeInMinutes))
      warning(
        "no durations occured at the first interval. Check if the 'node' argument is correctly specified"
      )
    
    for(i in 1:(length(nodes)-1)){  
      timeInterval <- timeInterval + ifelse((timeInMinutes>=nodes[i] & timeInMinutes<nodes[i+1]),(nodes[i] + nodes[i+1])/2,0)  #all observations are given its mid interval value 
    }  
    
    timeInterval <- timeInterval + ifelse(timeInMinutes == nodes[length(nodes)], (nodes[length(nodes)-1] + nodes[length(nodes)])/2,0)
    
    if(aggregation == "all"){
      
      meandur <- plyr::ddply(data.frame(durations = dur$durations, timeInterval = timeInterval), plyr::.(timeInterval), plyr::summarize,              
                     mean = round(mean(durations, na.rm=TRUE), 2))
      
      if(nrow(meandur) < 4) stop("Needs data in at least 3 nodes")
      spline <- splnFnc(meandur$timeInterval, meandur$mean)      
      if(returnSplineFnc == TRUE) rtSpline <- spline
      adjDur <- dur$durations/stats::predict(spline, timeInMinutes)$y
      
      df = data.frame(spline=stats::predict(spline, seq(nodes[1], nodes[length(nodes)], 1)))
      
      g <- ggplot(df, aes(x=spline.x/60,y=spline.y))
      g <- g + geom_hline(yintercept = 0, color="red")
      graphics::plot(g + geom_line()+ylab("Durations (seconds)")+xlab("time of the day")+ggtitle("Diurnal pattern estimated by a cubic spline function"))
      
    } else if(aggregation == "weekdays"){
      
      dur$timeInterval <- timeInterval
      
      adjDur <- numeric(nrow(dur))   
      df <- data.frame()
      for(j in 1:5){  
        meandur <- plyr::ddply(dur[dur$time$wday == j, -1], plyr::.(timeInterval), plyr::summarize,              
                         mean = round(mean(durations, na.rm=TRUE), 2))
        meandur <- cbind(meandur,day=j)        
        if(nrow(meandur) < 4) stop("Needs data in at least 3 nodes (too few observations on day ", day, " of the week)")
        spline <- splnFnc(meandur$timeInterval, meandur$mean)
        if(returnSplineFnc == TRUE) rtSpline <- c(rtSpline , list(spline))
        
        df = rbind(df, data.frame(spline=stats::predict(spline, seq(nodes[1], nodes[length(nodes)], 1)), day = j))
        adjDur <- adjDur + ifelse(dur$time$wday == j, dur$durations/stats::predict(spline, timeInMinutes)$y, 0)        
      }
      
      df$day <- factor(df$day, labels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday"))
      
      g <- ggplot(df , aes(x=spline.x/60,y=spline.y))
      g <- g + geom_line()
      g <- g + ylab("Durations (seconds)")
      g <- g + xlab("time of the day")
      g <- g + ggtitle("Diurnal pattern estimated by a cubic spline function")
      g <- g + facet_wrap(~day,ncol=5)
      g <- g + geom_hline(yintercept = 0, color="red")
      graphics::plot(g)
      
    } else if(aggregation == "none"){
      
      dur <- cbind(dur, timeInterval, date = strftime(dur$time, format = "%Y-%m-%d"))
      
      df <- data.frame()   
      adjDur <- numeric(nrow(dur))
      allDates <- unique(dur$date)
      for(date in allDates){
        meandur <- plyr::ddply(dur[dur$date == date, -1], plyr::.(timeInterval), plyr::summarize,              
                       mean = mean(durations, na.rm=TRUE))
        meandur <- cbind(meandur,date=date)
        if(nrow(meandur) < 4) stop("Needs data in at least 3 nodes (too few observations on ", date, ")")
        spline <- splnFnc(meandur$timeInterval, meandur$mean)
        if(returnSplineFnc == TRUE) rtSpline <- c(rtSpline , list(spline))
        
        df = rbind(df, data.frame(spline = stats::predict(spline, seq(nodes[1], nodes[length(nodes)], 1)), date = date))
        adjDur <- ifelse(dur$date == date, dur$durations/stats::predict(spline, timeInMinutes)$y, adjDur)         
      }
      
      
      #fix for ggplot2 in case the data are not starting on a monday (adds rows with empty dates in start):
      space = ""
      i <- 1
      while(i < dur$time[1]$wday){
        df <- rbind(data.frame(spline.x=rep(0, 1000), spline.y=rep(0, 1000), date = (space <- paste(space, " "))), df)
        i <- i +1 
      }
      
      g <- ggplot(df, aes(x = spline.x / 60, y = spline.y))
      g <- g + geom_line()
      g <- g + ylab("Durations (seconds)")
      g <- g + xlab("time of the day")
      g <- g + ggtitle("Diurnal pattern estimated by a cubic spline function")
      g <- g + facet_wrap(~date,ncol=5)
      g <- g + geom_hline(yintercept = 0, color="red")
      graphics::plot(g)
    } else stop("The aggregation argument is not supported")
  } else if(method == "supsmu"){
    if(aggregation == "all"){
      smooth <- stats::supsmu((dur$time$hour * 3600 + dur$time$min * 60 + dur$time$sec), dur$durations, span = span)
      df <- data.frame(smooth)
      smooth <- data.frame(smooth, row.names = 1)
      
      adjDur <- dur$durations/smooth[paste(dur$time$hour * 3600 + dur$time$min * 60 + dur$time$sec), ]
      
      g <- ggplot(df, aes(x=x/3600,y=y))
      g <- g + geom_line()
      g <- g + ylab("Durations (seconds)")
      g <- g + xlab("time of the day")
      g <- g + ggtitle("Diurnal pattern estimated by \"super smoother\"")
      g <- g + geom_hline(yintercept = 0, color="red")
      graphics::plot(g)
    } else if(aggregation == "weekdays"){
      
      df <- data.frame()
      adjDur <- numeric(nrow(dur))
      for(j in 1:5){ 
        tempTime <- dur$time[dur$time$wday == j]
        smooth <- stats::supsmu((tempTime$hour * 3600 + tempTime$min * 60 + tempTime$sec), dur$durations[dur$time$wday == j], span = span)
        df <- rbind(df, data.frame(smooth, day = j))
        smooth <- data.frame(smooth, row.names = 1)
        adjDur[dur$time$wday == j] <- dur$durations[dur$time$wday == j]/smooth[paste(tempTime$hour * 3600 + tempTime$min * 60 + tempTime$sec), ]
      }         
      
      df$day <- factor(df$day, labels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday"))
      
      g <- ggplot(df, aes(x=x/3600,y=y))
      g <- g + geom_line()
      g <- g + ylab("Durations (seconds)")
      g <- g + xlab("time of the day")
      g <- g + ggtitle("Diurnal pattern estimated by \"super smoother\"")
      g <- g + facet_wrap(~day,ncol=6)
      g <- g + facet_grid(~day)
      g <- g + geom_hline(yintercept = 0, color="red")
      graphics::plot(g)
      
    } else if(aggregation == "none"){
      
      df <- data.frame()            
      adjDur <- numeric(nrow(dur))
      dateTemp <- strftime(dur$time, format = "%Y-%m-%d")
      for(date in unique(dateTemp)){
        tempTime <- dur$time[dateTemp == date]
        smooth <- stats::supsmu((tempTime$hour * 3600 + tempTime$min * 60 + tempTime$sec), dur$durations[dateTemp == date], span = span)        
        df <- rbind(df, data.frame(smooth, date = date))
        smooth <- data.frame(smooth, row.names = 1)
        adjDur[dateTemp == date] <- dur$durations[dateTemp == date]/smooth[paste(tempTime$hour * 3600 + tempTime$min * 60 + tempTime$sec), ]
      }
      
      #fix for ggplot2 in case the data are not starting on a monday (adds rows with empty dates in start):
      space = ""
      i <- 1
      while(i < dur$time[1]$wday){
        df <- rbind(data.frame(x=rep(43200, 2), y=rep(0, 2), date = (space <- paste(space, ""))), df)
        i <- i +1 
      }
      
      g <- ggplot(df, aes(x = x / 3600, y = y))
      g <- g + geom_line()
      g <- g + ylab("Durations (seconds)")
      g <- g + xlab("time of the day")
      g <- g + ggtitle("Diurnal pattern estimated by \"super smoother\"")
      g <- g + facet_wrap(~date,ncol=5)
      g <- g + geom_hline(yintercept = 0, color="red")
      graphics::plot(g)
      
    } else stop("The aggregation argument is not supported")
  } else if(method == "FFF"){
    timeInSec <- 3600*dur$time$hour + 60*dur$time$min + dur$time$sec
    range <- max(timeInSec) - min(timeInSec)
    tBar <- (timeInSec - min(timeInSec)) / range
    
    if(aggregation == "all"){
      
      deltaC <- matrix(NA, nrow(dur), Q); deltaS <- matrix(NA, nrow(dur), Q)    
      for(j in 1:Q){
        deltaC[, j] <- cos(tBar * 2*pi*j)
        deltaS[, j] <- sin(tBar * 2*pi*j)
      }
      
      OLSest <- stats::lm(dur$durations ~ tBar + deltaC + deltaS)
      
      adjDur <- stats::predict(OLSest)
      adjDurDF <- data.frame(time = timeInSec, adjDur)
      
      #       fff <- function(x, coff = OLSest$coefficients, Q = Q){
      #         out <- coff[1] + x * coff[2]
      #         for(j in 1:Q){
      #           out <- out + coff[j + 2] * cos(x*2*pi*j) + coff[Q + j + 2] * sin(x*2*pi*j)
      #         }
      #         fff <- out
      #       }    
      #plot(seq(0,1,.01),fff(seq(0,1,.01), coff = OLSest$coefficients, Q = Q), t="l")
      
      g <- ggplot(adjDurDF, aes(x = time/3600 ,y = adjDur))
      g <- g + geom_line()
      g <- g + ylab("Durations (seconds)")
      g <- g + xlab("time of the day")
      g <- g + ggtitle("Diurnal pattern estimated by \"Flexible Fourier Form\"")
      g <- g + geom_hline(yintercept = 0, color="red")
      graphics::plot(g)      
    } else if(aggregation == "weekdays"){
      df <- data.frame()
      adjDur <- numeric(nrow(dur))
      for(k in 1:5){ 
        tempTBar <- tBar[dur$time$wday == k]
        deltaC <- matrix(NA, length(tempTBar), Q); deltaS <- matrix(NA, length(tempTBar), Q)    
        for(j in 1:Q){
          deltaC[, j] <- cos(tempTBar * 2*pi*j)
          deltaS[, j] <- sin(tempTBar * 2*pi*j)
        }
        
        OLSest <- stats::lm(dur$durations[dur$time$wday == k] ~ tempTBar + deltaC + deltaS)
        df <- rbind(df, data.frame(y = stats::predict(OLSest), time = timeInSec[dur$time$wday == k], day = k))
        adjDur[dur$time$wday == k] <- dur$durations[dur$time$wday == k]/stats::predict(OLSest)
      }
      
      
      df$day <- factor(df$day, labels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday"))
      
      g <- ggplot(df, aes(x=time/3600,y=y))
      g <- g + geom_line()
      g <- g + ylab("Durations (seconds)")
      g <- g + xlab("time of the day")
      g <- g + ggtitle("Diurnal pattern estimated by \"Flexible Fourier Form\"")
      g <- g + facet_wrap(~day,ncol=5)
      g <- g + geom_hline(yintercept = 0, color="red")
      graphics::plot(g)      
      
    } else if(aggregation == "none"){
      
      df <- data.frame()            
      adjDur <- numeric(nrow(dur))
      dateTemp <- strftime(dur$time, format = "%Y-%m-%d")
      
      for(date in unique(dateTemp)){
        tempTBar <- tBar[dateTemp == date]
        
        deltaC <- matrix(NA, length(tempTBar), Q); deltaS <- matrix(NA, length(tempTBar), Q)    
        for(j in 1:Q){
          deltaC[, j] <- cos(tempTBar * 2 * pi * j)
          deltaS[, j] <- sin(tempTBar * 2 * pi * j)
        }
        
        OLSest <- stats::lm(dur$durations[dateTemp == date] ~ tempTBar + deltaC + deltaS)
        df <- rbind(df, data.frame(y = stats::predict(OLSest), time = timeInSec[dateTemp == date], date = date))
        adjDur[dateTemp == date] <- dur$durations[dateTemp == date]/stats::predict(OLSest)         
      }
      
      #fix for ggplot2 in case the data are not starting on a monday (adds rows with empty dates in start):
      space = ""
      i <- 1
      while(i < dur$time[1]$wday){
        df <- rbind(data.frame(x=rep(0, 1000), y=rep(0, 1000), date = (space <- paste(space, " "))), df)
        i <- i +1 
      }
      
      g <- ggplot(df, aes(x=time/3600,y=y))
      g <- g + geom_line()
      g <- g + ylab("Durations (seconds)")
      g <- g + xlab("time of the day")
      g <- g + ggtitle("Diurnal pattern estimated by \"Flexible Fourier Form\"")
      g <- g + facet_wrap(~date,ncol=5)
      g <- g + geom_hline(yintercept = 0, color="red")
      graphics::plot(g) 
      
    } else stop("The aggregation argument is not supported")
    
  }  else stop("Method not supported")
  
  if(min(adjDur) <= 0) stop("The method and method arguments returned non-positive adjusted durations. Try a diffrent method or nodes etc.")     
  
  if(method == "cubicSpline" && returnSplineFnc ==  TRUE){
    diurnalAdj <- rtSpline
  } else{
    adjDur <- cbind(dur[, !(names(dur) %in% "adjDur")], adjDur)  #overwrites any previous "adjDur" column
    class(adjDur) <- c("durObj", "data.frame")
    attributes(adjDur$adjDur)$type <- attributes(dur)$type 
    attributes(adjDur$adjDur)$method <- method
    attributes(adjDur$adjDur)$aggregation <- aggregation
    attributes(adjDur$adjDur)$methodArguments <- switch(method,
                                                        cubicSpline = nodes, 
                                                        supsmu = span,
                                                        smoothSpline = c(nodes, spar),
                                                        FFF = Q)
    
    diurnalAdj <- adjDur
  }    
}