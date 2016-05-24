context("replacement functions for streamdef()")

# importing the discharges at gauge "Wildungsmauer (Danube)"
# truncating the data to complete hydrological years
infile <- readlfdata("QMittelTag207373.txt", type="HZB", hyearstart = 4,
                     baseflow = FALSE)
wild <- subset(infile, hyear %in% 1996:2011)


test_that("94% quantile for Wildungsmauer is correct", {
  expect_equal(unname(Qxx(wild, 94)), 988.74, tolerance = 1e-2)
})


arglist <- list(list(),
                list(pooling = "MA"),
                list(pooling = "MA", MAdays = 20),
                #list(pooling = "IC"),
                #list(pooling = "IC", IClevel = 0.5),
                #list(pooling = "IT"),
                #list(pooling = "IT", tmin = 10),
                list(threslevel = 80),
                list(threslevel = 95),
                #list(thresbreaks = "daily"), new version uses %j
                list(thresbreaks = "monthly")
                #list(thresbreaks = "seasonal"),
                #list(thresbreaks = "seasonal", breakdays = c("01/09")),
                #list(thresbreaks = "seasonal", breakdays = c("01/09", "08/08", "22/02")) # new version uses %j
)

streamdef.old <- function(lfobj,
                       pooling = c("none", "MA", "IT", "IC"),
                       threslevel = 70,
                       thresbreaks = c("fixed", "monthly", "daily", "seasonal"),
                       breakdays = c("01/06", "01/10"),
                       MAdays = 7,
                       tmin = 5,
                       IClevel = 0.1,
                       mindur = 0, #min days for IT/IC
                       minvol = 0, #min volume for IT/IC
                       table = c("all", "volmax", "durmax"),
                       na.rm = TRUE){
  lfcheck(lfobj)
  pooling <- match.arg(pooling)
  thresbreaks <- match.arg(thresbreaks)
  table <- match.arg(table)

  threshold <- buildthres(lfobj = lfobj, threslevel = threslevel,
                          thresbreaks=thresbreaks, breakdays = breakdays,
                          na.rm = na.rm)

  #Calculation FLOW
  if(pooling == "MA"){
    ma2 <- function(x, n){
      a <- ma(x=x, n = n, sides=2)

      # to avoid NA's -> first and last values are not averaged
      a[is.na(a)]<-x[is.na(a)]

      return(a)
    }
  }

  flow <- switch(pooling,
                 "none" = lfobj$flow,
                 "IT" = lfobj$flow,
                 "IC" = lfobj$flow,
                 "MA" = as.vector(ma2(x = lfobj$flow, n = MAdays)))

  #Comparing FLOW and THRESHOLD
  temp <- merge(x=lfobj, y=threshold, by = c("day", "month"), sort = FALSE)
  temp <- temp[order(temp$year, temp$month, temp$day), ]

  defrun <-  rle(flow <= temp$flow.y)
  pos <- c(cumsum(c(1, defrun$lengths)))
  streamdef <- data.frame(matrix(ncol = 9))
  names(streamdef) <-  c("pos", "d", "v", "mi", "Qmin", "startyear", "startmonth", "startday", "hyear")
  #a <-  pos[defrun$values] # length(defrun$lengths) != length(pos)
  a <-  pos[which(defrun$values)]
  a <- a[!is.na(a)]
  # a <- a[-length(a)] # switched off (GL)
  d <- defrun$lengths[defrun$values]
  d <- d[!is.na(d)]
  def <- flow - temp$flow.y
  for(ii in seq_along(a)){
    streamdef[ii, 1] <- a[ii]
    streamdef[ii, 6] <- lfobj[a[ii], "year"]
    streamdef[ii, 7] <- lfobj[a[ii], "month"]
    streamdef[ii, 8] <- lfobj[a[ii], "day"]
    streamdef[ii, 9] <- lfobj[a[ii], "hyear"]
    streamdef[ii, 2] <- d[ii]
    # bug: a[ii] can be out of range in the following line
    streamdef[ii, 3] <- sum(def[a[ii]:(a[ii]+d[ii]-1)])
    #  streamdef[ii, 4] <- streamdef[ii, 3]/d[ii]
    #  streamdef[ii, 5] <- min(flow[a[ii]:(a[ii]+d[ii]-1)])
  }

  if(pooling == "IT"){
    # buggy! inter-event time is neglected
    for(ii in nrow(streamdef):2){
      if(streamdef[ii, "pos"]-streamdef[ii-1, "d"]-tmin < streamdef[ii-1, "pos"]){
        streamdef[ii-1, "d"] <- streamdef[ii, "d"]+streamdef[ii-1, "d"]
        streamdef[ii-1, "v"] <- streamdef[ii, "v"]+streamdef[ii-1, "v"]
        streamdef <- streamdef[-ii, ]
      }
    }
  }

  if(pooling == "IC"){
    # buggy! max 2 events are merged
    for(ii in nrow(streamdef):2){
      if(streamdef[ii, "pos"]-streamdef[ii-1, "d"]-tmin < streamdef[ii-1, "pos"]){
        s <- sum((def[(streamdef[ii-1, "pos"]+streamdef[ii-1, "d"]):(streamdef[ii, "pos"]-1)]))
        if(-s/streamdef[ii-1, "v"]<IClevel){
          streamdef[ii-1, "d"] <- streamdef[ii, "pos"]-streamdef[ii-1, "pos"]+streamdef[ii, "d"]
          streamdef <- streamdef[-ii, ]
        }
      }
    }
  }

  for(ii in seq_along(streamdef$pos)){
    if(pooling == "IC"){
      streamdef[ii, 3] <- sum(def[streamdef[ii, "pos"]:(streamdef[ii, "pos"]+streamdef[ii, "d"]-1)])}
    streamdef[ii, 4] <- streamdef[ii, "v"]/streamdef[ii, "d"]
    streamdef[ii, 5] <- min(flow[streamdef[ii, "pos"]:(streamdef[ii, "pos"]+streamdef[ii, "d"]-1)])
  }

  #excluding minor deficits in IC/IT methods
  if(pooling %in% c("IC", "IT")){
    streamdef <- streamdef[streamdef$d >mindur & streamdef$v < -minvol, ]
  }

  streamdef$v <- -streamdef$v*60**2*24
  streamdef$mi <- -streamdef$mi*60**2*24
  if(table == "all") {
    return(streamdef[, -c(1, 9)])
  }else{
    if(table == "volmax"){
      agg <- aggregate(v ~ hyear, streamdef, max)
      lin <- NULL
      for(ii in seq_along(agg$hyear)){
        lin[ii] <- min(which(streamdef$hyear == agg$hyear[ii] & streamdef$v == agg$v[ii]))
      }
      streamdef[lin, -c(1)]
    }else{
      agg <- aggregate(d ~ hyear, streamdef, max)
      lin <- NULL
      for(ii in seq_along(agg$hyear)){
        lin[ii] <- min(which(streamdef$hyear == agg$hyear[ii] & streamdef$d == agg$d[ii]))
      }
      streamdef[lin, -c(1)]
    }
  }
}


test_streamdef <- function(args, data) {
  old <- do.call(streamdef.old, c(list(lfobj = data), args))
  new <- suppressMessages(do.call(streamdef, c(list(lfobj = data), args)))

  n <- nrow(old) == nrow(new)
  d <- all(abs(old$d - new$duration) < 1)
  v <- all(abs(old$v - old$volume) < 1)
  time <- all(with(old, as.Date(paste(startyear, startmonth, startday, sep = "-"))) == new$start)

  return(n & d & v & time)
}


# dataset without NAs
data(ray)
sapply(arglist, test_streamdef, data = ray)
sapply(arglist, test_streamdef, data = wild)

# dataset with NAs
#data(ngaruroro)
#sapply(arglist, test_streamdef, data = ngaruroro)
