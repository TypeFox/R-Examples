plotDescTrans <- function(trans, windowunit = "hours", window = 1){
  
  volume <- price <- n <- transactions <- time <- NULL
  
  windowunit <- match.arg(windowunit, c("secs", "mins", "hours", "days"))
  
  if(!("POSIXlt" %in% class(trans$time))) trans$time <- as.POSIXlt(trans$time)
  timeInterval <- switch(windowunit,
                       secs = as.numeric(trunc(trans$time, "secs") - trans$time$sec %% window),
                       mins = as.numeric(trunc(trans$time, "mins") - (trans$time$min %% window) * 60),
                       hours = as.numeric(trunc(trans$time, "hours") - (trans$time$hour %% window) * 3600),
                       days = as.numeric(trunc(trans$time, "days") - ((trans$time$yday - trans[1,1]$yday) %% window) * 86400))
  
  
  df <- cbind(dplyr::select(trans, volume, price), timeInterval = timeInterval)
  by_timeInterval <- dplyr::group_by(df, timeInterval)
  sumvol <- dplyr::summarise(by_timeInterval,
                             sumvol = sum(volume, na.rm = TRUE),
                             transactions = n())
  sumvol <- as.data.frame(sumvol)
  sumvol$time  <- as.POSIXlt(sumvol$timeInterval, origin = "1970-01-01")  
  
  g1 <- ggplot(trans, aes(x=time,y=price, group = time$year+time$yday)) + geom_line()+ylab("price") + ggtitle("Price")
  
  g2 <- ggplot(sumvol, aes(x=time,y=sumvol)) + geom_bar(stat = "identity")+ylab("volume") + ggtitle(paste("Volume traded per", window, windowunit))
  
  g3 <- ggplot(sumvol, aes(x=time,y=transactions)) + geom_bar(stat = "identity")+ylab("transactions") + ggtitle(paste("Number of transactions per", window, windowunit))
  
  print(g1)
  graphics::par(ask = TRUE)
  print(g2)
  print(g3)
  graphics::par(ask = FALSE)
#   grid::grid.newpage()
#   grid::pushViewport(grid::viewport(layout = grid::grid.layout(3, 1)))  
#   
#   print(g1, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
#   print(g2, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
#   print(g3, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1))
}