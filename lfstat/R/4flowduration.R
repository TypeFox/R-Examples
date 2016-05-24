is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
  if(is.character(x)){return(FALSE)}
  abs(x - round(x)) < tol}

#########################
# Flow Duration Curve   #
#########################


fdc <- function(lfobj,
                year = "any",
                breakdays = NULL,
                colors = TRUE,
                xnorm = FALSE,
                ylog = TRUE,
                legend = TRUE,
                separate = FALSE,
                ...){
  lfcheck(lfobj)

  
  if(length(year) > 1){
    if(!(is.wholenumber(min(year))&is.wholenumber(max(year)))){
      stop("Year must be a number")} else{
        dummi <- year
        lfobj <- subset(lfobj, subset = hyear %in% min(dummi):max(dummi))
      }
  }
  if(length(year) == 1){
    if(!(year == "any" | is.wholenumber(year))){
      stop('Year must be "any", a whole number or a vector of two whole numbers to define a range')
    }
    if(is.wholenumber(year)){
      dummi <- year
      lfobj <- subset(lfobj, subset  = hyear == dummi)
    }  
  }
  
#   if(!(year %in% c(min(lfobj$hyear):max(lfobj$hyear),"any"))){
#    stop("'year must be within the range of your data or \"any\" for taking whole series")}

   bdays <- NULL
  
   #If there is a single breakday, take hyearstart as the second!  
   if(!is.null(breakdays)){
     ii <- subset(lfobj, year != hyear, month)
      if(nrow(ii)==0){hyearstart <-1} else if(max(ii) <5.5){hyearstart <- max(ii)+1}else{hyearstart <- min(ii)}
     if(length(breakdays) == 1){
        days = c(paste("1/",hyearstart,sep = ""),breakdays)
      } else {
        days = breakdays}
    bdays <- data.frame(matrix(ncol = 2))
    names(bdays) <- c("day", "month")
    str <- strsplit(days,"/")
    for(ii in seq_along(str)){
      bdays[ii,] <- as.numeric(c(str[[ii]]))
      }
   }
   #Building a Year/Season-Table
   Year <- data.frame(day = rep(1:31,12),month = sort(rep(1:12,31)))
   Year$breakp <- FALSE
   for(ii in seq_along(bdays$day)){
    Year$breakp[Year$day == bdays[ii,"day"] & Year$month == bdays[ii,"month"]]<-TRUE
    }
   Year$season = cumsum(Year$breakp)
   Year$season[Year$season == 0] <- max(Year$season)
   lfobj <- merge(lfobj, Year, by = c("day", "month"),sort = FALSE)
 
   plotdata <- aggregate(flow ~ season, data = lfobj, quantile, probs = seq(1,0,by = -0.01))[,-1]

  if(is.null(breakdays))
    {rownames(plotdata) <- "FDC-quantiles"} else{
    #bdays$order <- bdays$month + 12*(hyearstart > bdays$month)
    #bdays <- bdays[order(bdays$order),] #DK, 21.11.2011
    bdays <- bdays[order(bdays$month),] #DK 2.5.2012
    rownames(plotdata) <- 1:nrow(plotdata)
    bdays[nrow(bdays)+1,] <- bdays[1,]
    for (ii in 1:nrow(plotdata)){
     rownames(plotdata)[ii] <- paste("Season from ", bdays$day[ii],"/",bdays$month[ii], " to ", bdays$day[ii+1],"/",bdays$month[ii+1],sep = "")
    }}

  if(xnorm){
    x <- qnorm(seq(0.01,0.99,by = 0.01))
    plotdata <- subset(plotdata,select = 2:100)
    ranx = qnorm(c(0.01,.99))}
          else {
               x <- 0:100
               ranx <- c(0,100)}
    
 if(!separate){
  plot(rep(x,nrow(plotdata)),plotdata,
     type = "n",
     log = if(ylog){"y"}else{""},
     xlab = "Exceedance frequency (%)",
     ylab = lflabel("Flow"),
     xlim = ranx,
     xaxs = "i",
     xaxt = "n",
     ...)

  if(xnorm){
    axis(1,at = qnorm(c(.01,.05,seq(.1,.9,by = .1),.95,.99)), labels = c(1,5,seq(10,90,by = 10),95,99))
  } else{axis(1)}
    
    
  if(colors){
   for(ii in 1:nrow(plotdata)){
     points(x,plotdata[ii,],type = "l", col = ii)}
     if(legend){
     legend(x = "topright", legend = rownames(plotdata),lty = 1, col = 1:nrow(plotdata))}
   } else {
     for(ii in 1:nrow(plotdata)){
     points(x,plotdata[ii,],type = "l", lty = ii)
     }
        if(legend){
     legend(x = "topright", legend = rownames(plotdata), lty = 1:nrow(plotdata))}}

}else{#separate = TRUE
       xpar <- par(ask = TRUE)
       
       for(ii in 1:nrow(plotdata)){
           plot(x,plotdata[ii,],
                type = "n",
                log = if(ylog){"y"}else{""},
                xlab = "Exceedance frequency (%)",
                ylab = lflabel("Flow"),
                xlim = ranx,
                xaxs = "i",
                xaxt = "n",
                ...)
        if(xnorm){
          axis(1,at = qnorm(c(.01,.05,seq(.1,.9,by = .1),.95,.99)), labels = c(1,5,seq(10,90,by = 10),95,99))
        } else{axis(1)}
    
        if(colors){
        points(x,plotdata[ii,],type = "l", col = 1)
           if(legend){
             legend(x = "topright", legend = rownames(plotdata)[ii],lty = 1, col = 1)}
         } else {
          points(x,plotdata[ii,],type = "l", lty = 1)
             if(legend){
     legend(x = "topright", legend = rownames(plotdata)[ii], lty = 1)}}    
         }
           par(xpar)
         }
  
  plotdata
   }
