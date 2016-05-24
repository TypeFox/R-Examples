# wrapper around the new function find_droughts() including pooling
# just emulates the interface of the old and buggy function streamdef()
streamdef <- function(lfobj,
                      pooling = c("none", "MA", "IT", "IC"),
                      threslevel = 70,
                      thresbreaks = c("fixed", "monthly", "daily","seasonal"),
                      breakdays = c("01/06","01/10"),
                      MAdays = 7,
                      tmin = 5,
                      IClevel = 0.1,
                      mindur = 0,
                      minvol = 0,
                      table = c("all", "volmax", "durmax"),
                      na.rm = TRUE
) {

  message("This is a wrapper around the new function 'find_droughts()' including pooling. Please use 'find_droughts()' directly, as in future versions 'streamdef' will be removed.")

  pooling <- match.arg(pooling)
  thresbreaks <- match.arg(thresbreaks)
  table <- match.arg(table)

  breakdays <- paste(breakdays, "2000", sep = "/")
  breakdays <- as.Date(breakdays, format = "%d/%m/%Y")

  if(thresbreaks == "fixed") thresbreaks <- "constant"
  if(thresbreaks == "seasonal") {
    thresbreaks <- breakdays
  }

  x <- as.xts(lfobj)

  threshold <- vary_threshold(x = x, varying = thresbreaks,
                              fun = quantile, probs = (1 - threslevel/100),
                              na.rm = na.rm)

  x <- find_droughts(x, threshold = threshold)

  if(pooling == "MA") x <- pool_ma(x, n = MAdays)
  if(pooling == "IT") x <- pool_ic(x, tmin = tmin, ratio = Inf)
  if(pooling == "IC") x <- pool_ic(x, tmin = tmin, ratio = IClevel)

  y <- summary(x, drop_minor = c(volume = minvol, duration = mindur))
  qmin <- tapply(x$discharge, x$event.no, min, na.rm = T)

  yy <- with(y, data.frame(d = duration, v = volume, mi = volume/duration,
                           Qmin = qmin[-1],
                           strsplit_date(start, prefix = "start")))

  if(table != "all") {
    yy$hyear <- water_year(y$start, hyear_start(lfobj))
    yy <- by(yy, yy$hyear, function(x) x[which.max(x[, substr(table, 1, 1)]), ])
    yy <- do.call(rbind, yy)
    rownames(yy) <- NULL
  }

  return(yy)
}

#########################
#Streamflow plot        #
#########################

streamdefplot <- function(lfobj, year, threslevel = 70,
                          thresbreaks = c("fixed", "monthly", "daily", "seasonal"),
                          breakdays = c("01/06", "01/10")){
  if(!year %in% c(min(lfobj$hyear):max(lfobj$hyear), "any")){
    stop("'year must be within the range of your data or \"any\" for taking whole series")}

  thresbreaks <- match.arg(thresbreaks)
  threshold <- buildthres(lfobj=lfobj,
                          threslevel = threslevel,
                          thresbreaks = thresbreaks,
                          breakdays = breakdays)
  dummi <- year
  sublf <- subset(lfobj, hyear == dummi)
  full <- merge(sublf, threshold, by = c("day", "month"), sort = FALSE)

  plot(x=1:length(full$flow.x),
       y=full$flow.x,
       type = "n",
       xaxt = "n",
       xaxs = "i",
       xlab = "",
       ylab = lflabel("Flow"))

  #Higher resolution to reduce gaps when crossing the threshold
  if(thresbreaks != "daily"){
    approx.x <- approx(1:length(full$flow.x), full$flow.x, xout = seq(1, length(full$flow.x), by = 0.05))
    #approx.y <- approx(1:length(full$flow.y), full$flow.y, xout = seq(1, length(full$flow.x), by = 0.05))
    approx.y2 <- approx(1:length(full$flow.y), full$flow.y, xout = seq(1, length(full$flow.x), by = 0.05), method = "constant")
    level <- apply(cbind(approx.x$y, approx.y2$y), 1, min)


    polygon(c(approx.x$x, rev(approx.y2$x)), c(level, rev(approx.y2$y)), col = "lightblue")#, border = NA)
    points(1:length(full$flow.y), full$flow.y, type = "s", col = "red")
  }else{
    approx.x <- approx(1:length(full$flow.x), full$flow.x, xout = seq(1, length(full$flow.x), by = 0.05))
    approx.y <- approx(1:length(full$flow.y), full$flow.y, xout = seq(1, length(full$flow.x), by = 0.05))
    level <- apply(cbind(approx.x$y, approx.y$y), 1, min)
    polygon(c(approx.x$x, rev(approx.y$x)), c(level, rev(approx.y$y)), col = "lightblue")#, border = NA)
    points(1:length(full$flow.y), full$flow.y, type = "l", col = "red")
  }
  points(x=1:length(full$flow.x), y=full$flow.x, type = "l")

  months <- which(lfobj$day[lfobj$hyear == year] == 1)
  monthsex <-  which(subset(lfobj, hyear == dummi, month)[months[1], ] ==lfobj$month & lfobj$day == 1 & lfobj$hyear == (year +1))
  lab <-rbind(subset(x = lfobj,
                     subset = hyear == dummi,
                     select = c(month, year))[months, ],
              c(lfobj$month[monthsex], lfobj$year[monthsex]))
  months <- c(months, 366)
  label <- paste(lab$month, lab$year, sep = "/")
  while(length(label) < length(months)){
    label <- c(label, "")}
  axis(1, at = months, labels = label)
}


#BUILDING threshold
buildthres <- function(lfobj,
                       threslevel = 70,
                       thresbreaks = c("fixed", "monthly", "daily", "seasonal"),
                       breakdays = c("01/06", "01/10"),
                       na.rm = TRUE){
  thresbreaks <- match.arg(thresbreaks)
  threshold <- data.frame(day = rep(1:31, 12), month = sort(rep(1:12, 31)))

  if(thresbreaks == "fixed")
  {threshold$flow <- quantile(lfobj$flow, 1-threslevel/100, na.rm = na.rm)}
  if(thresbreaks == "monthly"){
    mon <- aggregate(flow~month, lfobj, quantile, probs = 1-threslevel/100)
    threshold <- merge(threshold, mon, sort = FALSE)
  }
  if(thresbreaks == "daily"){
    mon <- aggregate(flow~day + month, lfobj, quantile, probs = 1-threslevel/100)
    threshold <- merge(threshold, mon, sort = FALSE)
  }
  if(thresbreaks == "seasonal"){
    ii <- subset(lfobj, year != hyear, month)
    if(nrow(ii)==0){hyearstart <-1} else if(max(ii) <5.5){hyearstart <- max(ii)+1}else{hyearstart <- min(ii)}
    if(length(breakdays) == 1){
      breakdays <- c(paste("1/", hyearstart, sep = ""), breakdays)
    }
    bdays <- data.frame(matrix(ncol = 2))
    names(bdays) <- c("day", "month")
    str <- strsplit(breakdays, "/")
    for(ii in seq_along(str)){
      bdays[ii, ] <- as.numeric(c(str[[ii]]))
    }
    threshold$breakp <- FALSE
    for(ii in seq_along(bdays$day)){
      threshold$breakp[threshold$day == bdays[ii, "day"] & threshold$month == bdays[ii, "month"]]<-TRUE
    }
    threshold$season = cumsum(threshold$breakp)
    threshold$season[threshold$season == 0] <- max(threshold$season)
    threshold$breakp <- NULL
    something <- merge(lfobj, threshold, by = c("day", "month"), sort = FALSE)
    sea <- aggregate(flow ~ season, something, quantile, probs = 1-threslevel/100)
    threshold <- merge(threshold, sea, by = "season", sort = FALSE)
    threshold$season <- NULL
  }
  threshold}
