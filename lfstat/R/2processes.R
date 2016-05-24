#########################
#Seasonal barchart      #
#########################

#Works with NAs!
sbplot <-  function(lfobj,hyearorder = TRUE){
  lfcheck(lfobj)

  if(hyearorder){
    ii <- subset(lfobj, year != hyear, month)
    if(nrow(ii)==0){hyearstart <-1} else if(max(ii) <5.5){hyearstart <- max(ii)+1}else{hyearstart <- min(ii)}

    if(hyearstart == 1){
      monthreorder <-  1:12} else{
        monthreorder <- c((14-hyearstart):12,1:(13-hyearstart))}
  } else{monthreorder <- 1:12}


  mm <- aggregate(flow ~ month, data = lfobj, FUN = mean)
  mquan <- aggregate(flow ~ month, data = lfobj, FUN = quantile,probs = .1)
  mquan$monthreorder <- monthreorder
  mmplot <- barchart(flow~reorder(month,monthreorder),
                     data = mm,
                     box.ratio=1,
                     horizontal = FALSE,
                     origin = 0,
                     xlab = "Month",
                     ylab =lflabel("Streamflow"))+
    layer(panel.barchart(x=reorder(month,monthreorder),
                         y=flow,
                         box.ratio = .3,
                         origin = 0,
                         horizontal = FALSE,
                         col = "black"),
          data = mquan)
  mmplot
}
#########################
#Double-Mass curve      #
#########################

dmcurve <- function(x,y,year = "any",namex = substitute(x),namey = substitute(y),na.rm = TRUE){
  lfcheck(x)
  lfcheck(y)

  if(length(year) > 1){
    if(!(is.wholenumber(min(year))&is.wholenumber(max(year)))){
      stop("Year must be a number")} else{
        dummi <- year
        xlfobj <- subset(x, subset = hyear %in% min(dummi):max(dummi))
        ylfobj <- subset(y, subset = hyear %in% min(dummi):max(dummi))
      }
  }
  if(length(year) == 1){
    if(!(year == "any" | is.wholenumber(year))){
      stop('Year must be "any", a whole number or a vector of two whole numbers to define a range')
    }
    if(is.wholenumber(year)){
      dummi <- year
      xlfobj <- subset(x, subset  = hyear == dummi)
      ylfobj <- subset(y, subset  = hyear == dummi)
    }else{
      xlfobj <- x
      ylfobj <- y}
  }


  c <- seq(0,1,by = 0.01)
  xq <- quantile(xlfobj$flow, probs = c,na.rm = na.rm)
  yq <- quantile(ylfobj$flow, probs = c,na.rm = na.rm)

  plot(xq,yq,type = "l", xlab = lflabel(paste("Streamgauge", namex, "flow")),
       ylab = lflabel(paste("Streamgauge", namey, "flow")))
}

#########################
#Hydrograph             #
#########################

hydrograph <- function(lfobj, startdate = NULL, enddate = NULL, amin = FALSE,...){
  lfcheck(lfobj)

  if(is.null(startdate)){
    startdate <- paste(lfobj$day[1],lfobj$month[1],lfobj$year[1], sep = "/")
  }

  if(startdate %in% 0:3000){
    first <- min(which(lfobj$hyear == startdate))
    if(length(first)==0){stop("Startyear out of series range")}
    startdate <- paste("1",lfobj$month[first],lfobj$year[first], sep = "/")}

  if(is.null(enddate)){
    end <- nrow(lfobj)
    enddate <- paste(lfobj$day[end],lfobj$month[end],lfobj$year[end], sep = "/")
  }

  if(enddate %in% 0:3000){
    last <- max(which(lfobj$hyear == enddate))
    if(length(last)==0){stop("Endyear out of series range")}
    enddate <- paste(lfobj$day[last],lfobj$month[last],enddate, sep = "/")}

  start <- strsplit(startdate, split = "/")
  end <- strsplit(enddate, split = "/")

  a <- c(which(lfobj$day == as.numeric(start[[1]][1]) &
                 lfobj$month == as.numeric(start[[1]][2]) &
                 lfobj$year == as.numeric(start[[1]][3])),
         which(lfobj$day == as.numeric(end[[1]][1]) &
                 lfobj$month == as.numeric(end[[1]][2]) &
                 lfobj$year == as.numeric(end[[1]][3])))
  if(length(a) < 2) stop("Invalid dates, maybe out of series range")
  #+1 is for getting the last tick/label
  plusone <- max(a)<nrow(lfobj)
  sub <- lfobj[min(a):(max(a)+plusone),]
  xachse <- as.numeric(row.names(sub))
  monstart <- which(sub$day == 1)
  if(length(monstart) > 60){
    monstart <-which(sub$day == 1 & sub$month ==1)}
  year <- sub$year[monstart]
  months <- sub$month[monstart]

  label <- paste(months,year, sep = "/")
  #plot(x=1:(nrow(sub)-plusone), #-1 for label, see prev. comment
  #If last date in Plot, skipp last tickmark!
  if(plusone){
    xplot <- xachse[-length(xachse)]
    yflow <- sub$flow[-length(xachse)]}
  else{
    xplot <- xachse
    yflow <- sub$flow
  }

  plot(x=xplot,
       y= yflow,
       type = "l",
       xaxt = "n",
       xlim = c(min(xplot),max(xplot) + 1.2),
       xlab = "",
       xaxs = "i",
       ylab = lflabel("Flow"),
       ...)
  axis(1, at = xachse[monstart], labels = label)

  if(amin){
    table <- aggregate(flow ~ hyear, data = sub, min)
    amini <- NULL
    for(ii in seq_along(table$hyear)){
      amini[ii] <-  min(which(table[ii,1] == sub$hyear & table[ii,2] == sub$flow))
    }

    points(amini,table$flow,type = "p",col = 4,cex = .6)
  }
}
