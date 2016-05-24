which.min.na <- function(x) {
  idx <- which.min(x)
  if(!length(idx)) idx <- NA

  return(idx)
}


baseflow <- function(x, tp.factor = 0.9, block.len = 5) {
  x <- as.numeric(x)

  # filling a matrix with ncol = block.len (column = one block)
  # to prevent recycling, pad x with NAs
  # using "negative modulo" to obtain the needed number of NAs
  y <- matrix(c(x, rep(NA, -length(x) %% block.len)), nrow = block.len)

  # absoulte position of block minima (idx.min)
  # matrix + apply() is 30% faster than tapply(), because data is already sorted
  offset <- seq.int(0, ncol(y) - 1) * block.len
  idx.min <- apply(y, 2, which.min.na) + offset
  block.min <- x[idx.min]

  # towards high flows, allow the central value of three consecutive minimas
  # only to be of a factor (1-tp.factor) higher than the surrounding values
  # e.g. modify the central value
  cv.mod <- tp.factor * tail(head(block.min, -1), -1)

  # check if value is a turning point, shift by +/- 1
  # first value is never a turnig point because there is no observation before
  is.tp <- cv.mod <= tail(block.min, -2) & cv.mod <= head(block.min, -2)
  is.tp <- c(F, is.tp)

  # need at least two non-NA turning points to interpolate
  if (sum(is.finite(block.min[is.tp])) < 2) return(rep_len(NA, length(x)))

  # interpolate base flow only between turning points
  # values outside enclosing turning points become NA, because approx(rule = 1)
  bf <- approx(x = idx.min[is.tp], y = block.min[is.tp], xout = seq_along(x))$y

  # baseflow must be lower than actual discharge
  return (pmin(bf, x))
}

#Calculating BFI
BFI <- function(lfobj, year = "any",breakdays = NULL,yearly = FALSE){
  calcbfi <- function(lfobj){
    sum(lfobj$baseflow[!is.na(lfobj$baseflow)&!is.na(lfobj$flow)])/sum(lfobj$flow[!is.na(lfobj$baseflow)&!is.na(lfobj$flow)])
  }

  lfcheck(lfobj)
  if(!all(year %in% c(min(lfobj$hyear):max(lfobj$hyear),"any"))){
    stop("\"year\" must be within the range of your data or \"any\" for calculating the base flow index of the whole series")}
  dummi <- year
  if(!any(dummi == "any")){
    lfobj <- subset(lfobj, hyear %in% dummi)}

  if(!is.null(breakdays)){
    lfobj2 <- usebreakdays(lfobj,breakdays)

    if(!yearly){
      lfobj3 <- split(lfobj2,lfobj2$seasonname)
    } else{
      lfobj3 <- split(lfobj2,list(lfobj2$seasonname,lfobj2$hyear))
    }
    return(sapply(lfobj3,calcbfi))
  }
  if(!yearly){
    return(calcbfi(lfobj))}else{
      return(sapply(split(lfobj,list(lfobj$hyear)),calcbfi))
    }}

#Plotting
bfplot <- function(lfobj,
                   year = "any",
                   col = "green",
                   bfcol = "blue",
                   ylog = FALSE){
  lfcheck(lfobj)
  if(ylog){ln <- "y"}else{ln <- ""}

  if(year == "any"){
    plot(lfobj$flow,
         type = "l",
         col = col,
         xaxt = "n",
         xlab = "",
         ylab = lflabel("Flow"),
         log = ln,
         xaxs = "i",
         mgp = c(1.5,0.5,0),
         tcl = -.3)
    lines(lfobj$baseflow,
          col = bfcol)
    monstart <- which(lfobj$day == 1)
    if (length(monstart) > 60) {
      monstart <- which(lfobj$day == 1 & lfobj$month == 1)
    }
    year <- lfobj$year[monstart]
    months <- lfobj$month[monstart]
    label <- paste(months, year, sep = "/")
    axis(1, at = monstart, labels = label)
    grid(nx = NA, ny = NULL)
    abline(v = monstart,col = "lightgray", lty ="dotted")
    return()
  }

  plot(lfobj$flow[lfobj$hyear == year],
       type = "l",
       col = col,
       xaxt = "n",
       xlab = "",
       ylab = lflabel("Flow"),
       log = ln,
       xaxs = "i",
       mgp = c(1.5,0.5,0),
       tcl = -.3)
  dummi <- year
  months <- which(lfobj$day[lfobj$hyear == year] == 1)
  monthsex <-  which(subset(lfobj, hyear == dummi, month)[months[1], ] == lfobj$month & lfobj$day == 1 & lfobj$hyear == (year + 1))
  lab <-rbind(subset(x = lfobj, hyear == dummi, select = c(month, year))[months, ],
              c(lfobj$month[monthsex],lfobj$year[monthsex]))
  months <- c(months,366)
  label <- paste(lab$month,lab$year,sep = "/")

  axis(1,
       at = months,
       labels = label,
       mgp = c(1.5,0.5,0),
       tcl = -.3,cex.axis = 0.8)
  grid(nx = NA, ny = NULL)
  abline(v = months,col = "lightgray", lty ="dotted")
  lines(lfobj$baseflow[lfobj$hyear == year],
        col = bfcol)
}
#




#Recession constant (MRC and IRS method)

#segselect takes the flow-time series and gives back a data-frame of the first "seglenth" days after each >"seglength" days shortfall under a given "threshold"-level (Q70)
#Used in:
#MRC und #IRS


#Calculation "rainpeaks" (in fact, this peaks are not always rain peaks, but they need to be excluded for the calculation of the recession constant.

#GOOD METHOD?????
pcheck <- function(p){
  if((0<p&1>p)||is.logical(p)){
    p}else{
      stop("p must be in (0,1) or a logical vector")
    }
}


rainpeak <- function(x,p=0.95){
  if(is.logical(p)){
    if(length(x) == length(p)){
      return(p)}  else{
        (stop("p and lfobj$flow differ in length"))}
  }
  rp <- FALSE
  for(ii in 2:(length(x)-1)){
    rp[ii] <- (x[ii]*p >= x[ii-1] & x[ii]*p >= x[ii+1])
    #to remove NAs
    rp[is.na(rp)] <- FALSE
  }
  rp[length(x)] <- FALSE
  rp
}

recessionplot <- function(lfobj,
                          peaklevel = 0.95,
                          plot = TRUE,
                          peakreturn = FALSE,
                          thresplot = TRUE,
                          threscol = "blue",
                          threshold = 70,
                          thresbreaks = c("fixed","monthly","seasonal"),
                          thresbreakdays = c("01/06","01/10"),
                          recessionperiod = TRUE,
                          recessioncol = "darkblue",
                          seglength = 7,
                          ...){
  threslevel <- threshold
  rainpeaklevel <- peaklevel
  pcheck(rainpeaklevel)

  peaks <- rainpeak(lfobj$flow, p = rainpeaklevel)
  if(plot){
    hydrograph(lfobj, ...)
    sub <- lfobj[peaks,]
    xachse <- as.numeric(row.names(sub))
    points(xachse, sub$flow)
    if(thresplot){
      thresbreaks <- match.arg(thresbreaks)
      threshold <- buildthres(lfobj=lfobj,
                              threslevel = threslevel,
                              thresbreaks = thresbreaks,
                              breakdays = thresbreakdays)
      lfobj$row <- as.numeric(row.names(lfobj))
      full <- merge(lfobj, threshold, by = c("day","month"),sort = FALSE)
      full <- full[order(full$row),]
      points(x = full$row,full$flow.y,type = "l", col = threscol)
      if(recessionperiod){
        rain2day <- peaks
        temp <- full
        startpoint <- which(1 == diff(lfobj$flow < temp$flow.y &
                                        !(c(FALSE,rain2day[-length(rain2day)])&c(FALSE,(lfobj$flow > temp$flow.y)[-length(rain2day)])) & #the day before was no rainday > threshold
                                        !(c(FALSE,FALSE,rain2day[-c(length(rain2day)-1,length(rain2day))])& c(FALSE,FALSE,(lfobj$flow > temp$flow.y)[-c(length(rain2day)-1,length(rain2day))]))))

        segment <- rep(FALSE,length(lfobj$flow))
        segment[startpoint] <- TRUE
        dif <- diff(lfobj$flow)
        #Series goes on, if next value is smaller
        for(ii in 1:(length(segment)-1)){
          if(segment[ii])  segment[ii+1] <- (dif[ii] < 0)}

        for(ii in startpoint){
          if(!all(segment[ii:(ii+seglength-1)])){
            a <- TRUE
            for(jj in ii:(ii+seglength-1)){
              segment[jj] <- FALSE
              if(!segment[jj+1]) break
            }}}

        linecol = rep(recessioncol, length(temp$flow.x))
        linecol[!segment] <- 0
        points(temp$row,temp$flow.x,type = "p", col = linecol,pch = 18)
      }}}
  if(peakreturn)
    return(peaks)
}

seglenplot <- function(lfobj,
                       threslevel =70,
                       thresbreaks = c("fixed","monthly","seasonal"),
                       thresbreakdays = NULL,
                       rainpeaklevel = 0.95,
                       na.rm = TRUE){
  lfcheck(lfobj)
  rain2day <- rainpeak(lfobj$flow,rainpeaklevel)
  thresbreaks <- match.arg(thresbreaks)

  if(thresbreaks != "seasonal"){
    thres <- buildthres(lfobj = lfobj, threslevel = threslevel, thresbreaks = thresbreaks, na.rm = na.rm)
  }else{
    if(is.null(thresbreakdays)) stop("No thresbreakdays specified!")
    thres <- buildthres(lfobj = lfobj, threslevel = threslevel, thresbreaks = thresbreaks, breakdays = thresbreakdays, na.rm = na.rm)
  }
  check <- NULL
  temp <- merge(x=lfobj, y=thres, by = c("day","month"), sort = FALSE)
  temp <- temp[order(temp$year,temp$month,temp$day),]
  #  run <- rle(lfobj$flow < temp$flow.y & !rain2day & !c(FALSE,rain2day[-length(rain2day)]#) & !c(FALSE,FALSE,rain2day[-c(length(rain2day)-1,length(rain2day))]) &c(TRUE,diff(lfobj#$flow)<0))


  #Select Startpoints
  startpoint <- which(1 == diff(lfobj$flow < temp$flow.y &
                                  !(c(FALSE,rain2day[-length(rain2day)])&c(FALSE,(lfobj$flow > temp$flow.y)[-length(rain2day)])) & #the day before was no rainday > threslevel
                                  !(c(FALSE,FALSE,rain2day[-c(length(rain2day)-1,length(rain2day))])&
                                      c(FALSE,FALSE,(lfobj$flow > temp$flow.y)[-c(length(rain2day)-1,length(rain2day))]))))
  segment <- rep(FALSE,length(lfobj$flow))
  segment[startpoint] <- TRUE
  dif <- diff(lfobj$flow)
  #Series goes on, if next value is smaller
  for(ii in 1:(length(segment)-1)){
    if(segment[ii])  segment[ii+1] <- (dif[ii] < 0)}

  run <- rle(segment)

  #x11(width = 14, height = 7, title = "Recession duration (days)")
  tab <- table(run$length[run$value])
  splot<-barchart(tab[!(names(tab) %in% c("1","2","3"))],main = "Recession duration", xlab = paste("Days using Q", threslevel, " as threshold",sep = ""),horizontal = FALSE)
  splot
}



segselect <- function(lfobj,
                      threshold = 70,
                      seglength,
                      thresbreaks = NULL,
                      thresbreakdays = NULL,
                      p,
                      na.rm = TRUE,
                      season = FALSE){ #seglength "AUTO" to be solved!!!
  if(!seglength){stop("seglength missing")} ####!!!!! Better stop cond!!!
  #Excluding 2 days after rainpeak as defined by the function "rainpeak"
  rain2day <- rainpeak(lfobj$flow,p)

  if(season){
    thres <- buildthres(lfobj = lfobj, threslevel = threshold, thresbreaks = "fixed", na.rm = na.rm)
  }else{
    thres <- buildthres(lfobj = lfobj, threslevel = threshold, thresbreaks = thresbreaks, breakdays = thresbreakdays, na.rm = na.rm)
  }
  check <- NULL
  temp <- merge(x=lfobj, y=thres, by = c("day","month"), sort = FALSE)
  temp <- temp[order(temp$year,temp$month,temp$day),]
  #  run <- rle(lfobj$flow < temp$flow.y & !rain2day & !c(FALSE,rain2day[-length(rain2day)]#) & !c(FALSE,FALSE,rain2day[-c(length(rain2day)-1,length(rain2day))]) &c(TRUE,diff(lfobj#$flow)<0))


  #Select Startpoints
  startpoint <- which(1 == diff(lfobj$flow < temp$flow.y &
                                  !(c(FALSE,rain2day[-length(rain2day)])&c(FALSE,(lfobj$flow > temp$flow.y)[-length(rain2day)])) & #the day before was no rainday > threshold
                                  !(c(FALSE,FALSE,rain2day[-c(length(rain2day)-1,length(rain2day))])&
                                      c(FALSE,FALSE,(lfobj$flow > temp$flow.y)[-c(length(rain2day)-1,length(rain2day))]))))
  segment <- rep(FALSE,length(lfobj$flow))
  segment[startpoint] <- TRUE
  dif <- diff(lfobj$flow)
  #Series goes on, if next value is smaller
  for(ii in 1:(length(segment)-1)){
    if(segment[ii])  segment[ii+1] <- (dif[ii] < 0)}

  run <- rle(segment)


  pos <- c(cumsum(c(1,run$lengths)))
  lfs <- data.frame(matrix(ncol = seglength))
  a <-  pos[which(run$lengths >=seglength & run$values == TRUE)]
  for(ii in seq_along(a)){
    lfs[ii,] <- c(lfobj$flow[a[ii]:(a[ii]+seglength-1)])
  }
  if(all(is.na(lfs))){return("E")}
  lfs
}

recession <- function(lfobj,
                      method = c("MRC","IRS"),
                      seglength,
                      threshold,
                      peaklevel = 0.95,
                      seasonbreakdays = NULL,
                      thresbreaks = c("fixed","monthly","seasonal"),
                      thresbreakdays = NULL,
                      plotMRC = TRUE,
                      trimIRS = 0,
                      na.rm = TRUE){
  rainpeaklevel <- peaklevel
  pcheck(p = rainpeaklevel)
  lfcheck(lfobj)
  meth <-match.arg(method)
  thresbreaks = match.arg(thresbreaks)

  if(!is.null(seasonbreakdays)){
    b <- usebreakdays(lfobj,seasonbreakdays)
    b<- b[order(b$year,b$month,b$day),]
    ding <- split(b,b$seasonname)
    if(meth == "MRC" & plotMRC){
      message("No plot available for separated seasons")
    }
    res <- switch(meth,
                  MRC = sapply(ding,MRC, seglength = seglength, threshold = threshold,rainpeaklevel = rainpeaklevel, thresbreaks = thresbreaks, thresbreakdays = thresbreakdays,
                               plot = FALSE,na.rm = na.rm,season = TRUE),
                  IRS = sapply(ding,IRS, seglength = seglength, threshold = threshold, rainpeaklevel = rainpeaklevel, thresbreaks = thresbreaks, thresbreakdays = thresbreakdays,
                               na.rm = na.rm,trimlevel = trimIRS,season = TRUE))
  } else{
    res<-switch(meth,
                MRC = MRC(lfobj,threshold = threshold, seglength = seglength,
                          rainpeaklevel = rainpeaklevel,
                          thresbreaks = thresbreaks, thresbreakdays = thresbreakdays,
                          plot = plotMRC,na.rm = na.rm),
                IRS = IRS(lfobj,threshold = threshold, seglength = seglength,
                          rainpeaklevel = rainpeaklevel,
                          thresbreaks = thresbreaks,  thresbreakdays = thresbreakdays,
                          trimlevel = trimIRS,na.rm = na.rm)
    )
  }
  res}

#MRC-Method:
MRC <- function(lfobj,
                threshold,
                seglength,
                thresbreaks,
                thresbreakdays,
                rainpeaklevel,
                plot = TRUE,
                na.rm = TRUE,
                season = FALSE){

  lfs <- segselect(lfobj,threshold = threshold,seglength = seglength,
                   thresbreaks = thresbreaks,thresbreakdays = thresbreakdays,
                   p = rainpeaklevel,na.rm = na.rm,season = season)
  if(any(lfs == "E")){warning("There are no valid recession periodes, returning 'NA'");return(NA)}
  n <- data.frame(Qt = matrix(as.matrix(lfs[,1:(ncol(lfs)-1)])), Qtminus1 = matrix(as.matrix(lfs[,2:ncol(lfs)])))

  slope <- coef(lm(Qtminus1~Qt - 1,data = n))
  if(plot){
    plot(n,xlab = expression(Q[t-1]),ylab=expression(Q[t]),pch = 20)
    abline(a = 0,b = slope)
    legend("topleft", paste("y = ",round(slope,4),"x",sep = ""))
  }
  C <- -1/log(slope)
  names(C) <- NULL
  C
}

#Maybe using a 10% trimmed mean?! thinking about measure of spread
#Excluding C values < 0!?
IRS <- function(lfobj,
                threshold,
                seglength,
                trimlevel,
                thresbreaks,
                thresbreakdays,
                season = FALSE,
                rainpeaklevel,
                na.rm = TRUE){
  lfs <- segselect(lfobj,threshold = threshold,seglength = seglength,p = rainpeaklevel,
                   thresbreaks = thresbreaks,thresbreakdays = thresbreakdays, na.rm = na.rm,season = season)
  if(any(lfs == "E")){warning("There are no valid recession periodes, returning 'NA'");return(NA)}
  k <- NULL
  x <- as.numeric(1:(ncol(lfs)-1))
  for(ii in 1:nrow(lfs)){
    y <- as.numeric(log(lfs[ii,2:ncol(lfs)]/lfs[ii,1]))
    a <- try(coef(try(lm(y ~ x-1))))
    if(inherits(a, "try-error")){
      return(NA)}
    k[ii] <- a
  }
  C <- -1/k
  mean(C[C>0],trim = trimlevel)
}


#Mean Flow
meanflow <- function(lfobj,year = "any",monthly = FALSE,yearly = FALSE,breakdays = NULL,na.rm = TRUE){
  lfcheck(lfobj)

  if(!all(year %in% c(min(lfobj$hyear):max(lfobj$hyear),"any"))){
    stop("\"year\" must be within the range of your data or \"any\" for calculating the mean flow of the whole series")}
  dummi <- year
  if(!any(dummi == "any")){
    lfobj <- subset(lfobj, hyear %in% dummi)}

  if(!is.null(breakdays)){
    if(!yearly){
      return(aggregate(flow~seasonname,usebreakdays(lfobj,breakdays),mean,na.rm = na.rm))} else{
        return(aggregate(flow~seasonname+hyear,usebreakdays(lfobj,breakdays),mean,na.rm = na.rm))}
  }
  #Reordering the table according to the hyear
  if(monthly){
    ii <- subset(lfobj, year != hyear, month)
    if(nrow(ii)==0){hyearstart <-1} else if(max(ii) <5.5){hyearstart <- max(ii)+1}else{hyearstart <- min(ii)}

    if(hyearstart == 1){
      monthorder <-  1:12} else{
        monthorder <- c(hyearstart:12,1:(hyearstart-1))}
  }


  #  dummi <- year
  #  if(any(year == "any")){
  if(monthly){
    if(!yearly){
      aggregate(flow~month,lfobj,mean,na.rm = na.rm)[monthorder,]} else{
        aggregate(flow~month+hyear,lfobj,mean,na.rm = na.rm)}}else{
          if(!yearly){
            mean(lfobj$flow,na.rm = na.rm)} else{
              aggregate(flow~hyear,lfobj,mean,na.rm = na.rm)}}
}
#########################
#Q95                    #
#########################
Q95 <-  function(lfobj,year = "any",monthly = FALSE, yearly = FALSE, breakdays = NULL,na.rm = TRUE){
  Qxx(lfobj = lfobj, Qxx = 95,year = year, monthly = monthly, yearly = yearly, breakdays = breakdays, na.rm = na.rm)}

Q90 <-  function(lfobj,year = "any",monthly = FALSE, yearly = FALSE, breakdays = NULL,na.rm = TRUE){
  Qxx(lfobj = lfobj, Qxx = 90, year = year, monthly = monthly, yearly = yearly, breakdays = breakdays, na.rm = na.rm)}

Q70 <-  function(lfobj,year = "any",monthly = FALSE, yearly = FALSE, breakdays = NULL,na.rm = TRUE){
  Qxx(lfobj = lfobj, Qxx = 70,year = year, monthly = monthly, yearly = yearly, breakdays = breakdays, na.rm = na.rm)}




Qxx <- function(lfobj, Qxx, year = "any",
                monthly = FALSE, yearly = FALSE, breakdays = NULL,
                na.rm = TRUE){

  lfcheck(lfobj)

  if(!all(year %in% c(min(lfobj$hyear):max(lfobj$hyear),"any"))){
    stop("'year must be within the range of your data or \"any\" for calculating the Q95 of the whole series")
  }

  prob <- 1 - Qxx/100
  dummi <- year

  if(!any(dummi == "any")){
    lfobj <- subset(lfobj,hyear %in% dummi)
  }

  if(!is.null(breakdays)){
    if(!yearly){
      res <- aggregate(flow ~ seasonname, usebreakdays(lfobj, breakdays),
                       quantile, probs = prob, na.rm = na.rm)
    } else {
      res <- aggregate(flow ~ seasonname + hyear, usebreakdays(lfobj,breakdays),
                       quantile, probs = prob, na.rm = na.rm)
    }
  }

  #Reordering the table according to the hyear
  if(monthly){
    ii <- subset(lfobj, year != hyear, month)
    if(nrow(ii)==0){
      hyearstart <-1
    } else if(max(ii) <5.5) {
      hyearstart <- max(ii)+1
    } else {
      hyearstart <- min(ii)
    }


    if(hyearstart == 1){
      monthorder <-  1:12
    } else {
      monthorder <- c(hyearstart:12,1:(hyearstart-1))
    }
  }

  if(monthly){
    if(!yearly){
      res <- aggregate(flow~month,lfobj,quantile,probs = prob,na.rm = na.rm)[monthorder,]
    } else {
      res <- aggregate(flow~month+hyear,lfobj,quantile,probs = prob,na.rm = na.rm)
    }
  } else {
    if(!yearly){
      res <- quantile(lfobj$flow,probs = prob,na.rm = na.rm)
      names(res) <- paste0("Q", Qxx)
    } else {
      res <- aggregate(flow~hyear,lfobj,quantile,probs = prob,na.rm = na.rm)
    }
  }

  return(res)
}

#########################
#MAM                    #
#########################
MAannual <- function(lfobj, n=7, breakdays = NULL, year = "any"){
  lfobj$MAn <- ma(x = lfobj$flow, n = n)

  if(any(year != "any")){
    lfobj <- subset(lfobj, hyear %in% year)
  }

  if(is.null(breakdays)){
    annual <- aggregate(MAn ~ hyear, data = lfobj, FUN = min)
  } else {
    annual <- aggregate(MAn ~ seasonname + hyear,
                        data = usebreakdays(lfobj, breakdays),
                        FUN = min)
  }

  return (annual)
}

MAM <- function(lfobj,n=7,year = "any",breakdays = NULL,yearly = FALSE){
  lfcheck(lfobj)

  if(!is.null(breakdays)){
    if(!yearly){
      return(aggregate(MAn~seasonname,data = MAannual(lfobj=lfobj,n=n,breakdays = breakdays,year = year),FUN = mean))}else{
        return(MAannual(lfobj,n,breakdays,year = year))}
  }

  if(yearly){
    MAannual(lfobj,n,year = year)}else{
      mean(MAannual(lfobj,n,year = year)$MAn)}
}
#########################
#Seasonality Ratio      #
#########################

seasratio <- function(lfobj, breakdays, Q = 95){
  if(length(breakdays) == 1){
    ii <- subset(lfobj, year != hyear, month)
    if(nrow(ii)==0){hyearstart <-1} else if(max(ii) <5.5){hyearstart <- max(ii)+1}else{hyearstart <- min(ii)}
    breakdays = c(paste(1,hyearstart,sep = "/"),breakdays)
  }
  if(length(breakdays) >2){stop("Maximum number of breakdays allowed is 2")}

  threshold <- buildthres(lfobj=lfobj, threslevel = Q, thresbreaks = "seasonal", breakdays = breakdays)
  str <- strsplit(breakdays,"/")
  Q95z <- threshold[threshold$day == as.numeric(str[[1]][1]) & threshold$month ==as.numeric(str[[1]][2]),]$flow
  Q95n <- threshold[threshold$day == as.numeric(str[[2]][1]) & threshold$month ==as.numeric(str[[2]][2]),]$flow
  ratio <- Q95z/Q95n
  names(ratio) <- paste('Q95(',breakdays[1],'-',breakdays[2],')/Q95(',breakdays[2],'-',breakdays[1],')',sep = "")
  ratio
}

seasindex <- function(lfobj, Q=95,na.rm = TRUE){
  if(Q<0 | Q>100){stop("Q must be in [0,100]")}
  quan <- quantile(lfobj$flow,probs = 1-Q/100,na.rm = na.rm)
  ind <- which(lfobj$flow <= quan)
  a<-0
  theta <- NULL
  for(ii in ind){
    a <- a+1
    temp <- lfobj[ii,]
    t <- as.Date(paste(temp$year,temp$month,temp$day,sep = "-"))
    dz <- as.numeric(t - as.Date(paste(temp$year,01,01,sep = "-")))+1
    dn <- as.numeric(as.Date(paste(temp$year,12,31,sep = "-"))-as.Date(paste(temp$year,01,01,sep = "-")))+1
    theta[a] <- dz/dn*2*pi}
  xtheta <- mean(cos(theta))
  ytheta <- mean(sin(theta))
  thetafull <- atan2(ytheta,xtheta)
  if(thetafull < 0) {thetafull <- thetafull + 2*pi}
  D = thetafull * 365/(2*pi)
  r = sqrt(xtheta**2+ytheta**2)
  result <- list(theta = thetafull, D = D, r = r)
  result
}



#######################
#Season

usebreakdays <- function(lfobj,breakdays){
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
  #Season-Names
  bdays <- bdays[order(bdays$month),]
  names <- NULL
  bdays[nrow(bdays)+1,] <- bdays[1,]
  for (ii in seq_along(rownames(bdays))){
    lfobj$seasonname[lfobj$season==ii] <- paste("Season from ", bdays$day[ii],"/",bdays$month[ii], " to ", bdays$day[ii+1],"/",bdays$month[ii+1],sep = "")
  }
  lfobj[order(lfobj$season),]
}
