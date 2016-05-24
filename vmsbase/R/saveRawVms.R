
#' Raw VMS editing function
#' 
#' \code{saveRawVms} implements the routines that converts raw values
#'  to standard data.
#' 
#' @param rawfile    The raw VMS dataset.
#' @param widget    The widget list that contains the editing infos.
#'
#' @return The function returns the standardized VMS data.
#' 
#' @usage saveRawVms(rawfile, widget)
#' 
#' @export saveRawVms
#'@seealso \code{\link{gui_vms_editraw}}

saveRawVms <- function (rawfile, widget)
{
  
  widget <- gsub("\"", "", as.character(widget))
  vessIdSel <- widget[1]
  latModeSel <- widget[2]
  latDeg <- widget[3]
  latMin <- widget[4]
  latSec <- widget[5]
  latDir <- widget[6]
  latDec <- widget[7]
  lonModeSel <- widget[8]
  lonDeg <- widget[9]
  lonMin <- widget[10]
  lonSec <- widget[11]
  lonDir <- widget[12]
  lonDec <- widget[13]
  timeMode <- widget[14]
  timeUtc <- widget[15]
  
  timeDate2a <- widget[16]
  timeDate2b <- widget[17]
  
  timeDate <- widget[18]
  timeHour <- widget[19]
  timeMinute <- widget[20]
  timeSecond <- widget[21]
  speedModeSel <- widget[22]
  speedCol <- widget[23]
  headModeSel <- widget[24]
  headCol <- widget[25]
  
  if(widget[26] == "DD/MM/YYYY")
  {data_frm <-  c(dates = "d/m/y", times = "h:m:s")
  }else{
    data_frm <- c(dates = "m/d/y", times = "h:m:s")
  }
  
  numlines <- length(rawfile$data[,1])
  cat("\n\n   ---   Raw Pings Editing Started!   ---\n",
      "\nProcessing ", numlines, " raw pings...\n", sep = "")
  #to_rep <- paste(to_rep, "\nProcessed ", numlines, " raw pings.\n", sep = "")
  vmsdata <- data.frame("I_NCEE" = numeric(numlines),
                        "LAT" = numeric(numlines),
                        "LON" = numeric(numlines),
                        "DATE" = numeric(numlines),
                        "SPE" = numeric(numlines),
                        "HEA" = numeric(numlines))
  
  vessid <- which(colnames(rawfile$data) == vessIdSel)
  vmsdata[,"I_NCEE"] <- rawfile$data[,vessid]
  
  if (latModeSel == "sex")
  {
    
    latdeg <- which(colnames(rawfile$data) == latDeg)
    latmin <- which(colnames(rawfile$data) == latMin)
    latsec <- which(colnames(rawfile$data) == latSec)
    latdir <- which(colnames(rawfile$data) == latDir)
    tole_ladeg <- which(is.na(rawfile$data[,latdeg]))
    cat("\n   -   ", length(tole_ladeg), " NAs found in latitude degrees...", sep = "")
    #to_rep <- paste(to_rep, "\n   -   ", length(tole_ladeg), " NAs found in latitude degrees...\n", sep = "")
    if(length(tole_ladeg) > 0)
    {
      rawfile$data <- rawfile$data[-tole_ladeg,]
      vmsdata <- vmsdata[-tole_ladeg,]
    }
    tole_lamin <- which(is.na(rawfile$data[,latmin]))
    cat("\n   -   ", length(tole_lamin), " NAs found in latitude minutes...", sep = "")
    #to_rep <- paste(to_rep, "\n   -   ", length(tole_lamin), " NAs found in latitude minutes...\n", sep = "")
    if(length(tole_lamin) > 0)
    {
      rawfile$data <- rawfile$data[-tole_lamin,]
      vmsdata <- vmsdata[-tole_lamin,]
    }
    tole_lasec <- which(is.na(rawfile$data[,latsec]))
    cat("\n   -   ", length(tole_lasec), " NAs found in latitude seconds...", sep = "")
    #to_rep <- paste(to_rep, "\n   -   ", length(tole_lasec), " NAs found in latitude seconds...\n", sep = "")
    if(length(tole_lasec) > 0)
    {
      rawfile$data <- rawfile$data[-tole_lasec,]
      vmsdata <- vmsdata[-tole_lasec,]
    }
    tole_ladir <- which(is.na(rawfile$data[,latdir]))
    cat("\n   -   ", length(tole_ladir), " NAs found in latitude direction...", sep = "")
    #to_rep <- paste(to_rep, "\n   -   ", length(tole_ladir), " NAs found in latitude direction...\n", sep = "")
    if(length(tole_ladir) > 0)
    {
      rawfile$data <- rawfile$data[-tole_ladir,]
      vmsdata <- vmsdata[-tole_ladir,]
    }
    pre_lat <- latsex2dec(rawfile$data[,latdeg],rawfile$data[,latmin],rawfile$data[,latsec],rawfile$data[,latdir])
    
    tole_prelat <- which(pre_lat < -90 | pre_lat > 90)
    if(length(tole_prelat) > 0)
    {
      pre_lat <- pre_lat[-tole_prelat]
      vmsdata <- vmsdata[-tole_prelat,]
      rawfile$data <- rawfile$data[-tole_prelat,]
    }
    cat("\n   -   ", length(tole_prelat), " Latitudes out of range ( -90 / 90 )", sep = "")
    
    vmsdata[,"LAT"] <- pre_lat
    
  }
  
  if (latModeSel == "dec")
  {
    
    latdec <- which(colnames(rawfile$data) == latDec)
    tole_lade <- which(is.na(rawfile$data[,latdec]))
    cat("\n   -   ", length(tole_lade), " NAs found in decimal latitude...", sep = "")
    #to_rep <- paste(to_rep, "\n   -   ", length(tole_lade), " NAs found in decimal latitude...\n", sep = "")
    if(length(tole_lade) > 0)
    {
      rawfile$data <- rawfile$data[-tole_lade,]
      vmsdata <- vmsdata[-tole_lade,]
    }
    
    pre_lat <- rawfile$data[,latdec]
    
    tole_prelat <- which(pre_lat < -90 | pre_lat > 90)
    if(length(tole_prelat) > 0)
    {
      pre_lat <- pre_lat[-tole_prelat]
      vmsdata <- vmsdata[-tole_prelat,]
      rawfile$data <- rawfile$data[-tole_prelat,]
    }
    cat("\n   -   ", length(tole_prelat), " Latitudes out of range ( -90 / 90 )", sep = "")
    
    vmsdata[,"LAT"] <- pre_lat
    
  }
  
  
  if (lonModeSel == "sex")
  {
    
    londeg <- which(colnames(rawfile$data) == lonDeg)
    lonmin <- which(colnames(rawfile$data) == lonMin)
    lonsec <- which(colnames(rawfile$data) == lonSec)
    londir <- which(colnames(rawfile$data) == lonDir)
    tole_lodeg <- which(is.na(rawfile$data[,londeg]))
    cat("\n   -   ", length(tole_lodeg), " NAs found in longitude degrees...", sep = "")
    #to_rep <- paste(to_rep, "\n   -   ", length(tole_lodeg), " NAs found in longitude degrees...\n", sep = "")
    if(length(tole_lodeg) > 0)
    {
      rawfile$data <- rawfile$data[-tole_lodeg,]
      vmsdata <- vmsdata[-tole_lodeg,]
    }
    tole_lomin <- which(is.na(rawfile$data[,lonmin]))
    cat("\n   -   ", length(tole_lomin), " NAs found in longitude minutes...", sep = "")
    #to_rep <- paste(to_rep, "\n   -   ", length(tole_lomin), " NAs found in longitude minutes...\n", sep = "")
    if(length(tole_lomin) > 0)
    {
      rawfile$data <- rawfile$data[-tole_lomin,]
      vmsdata <- vmsdata[-tole_lomin,]
    }
    tole_losec <- which(is.na(rawfile$data[,lonsec]))
    cat("\n   -   ", length(tole_losec), " NAs found in longitude seconds...", sep = "")
    #to_rep <- paste(to_rep, "\n   -   ", length(tole_losec), " NAs found in longitude seconds...\n", sep = "")
    if(length(tole_losec) > 0)
    {
      rawfile$data <- rawfile$data[-tole_losec,]
      vmsdata <- vmsdata[-tole_losec,]
    }
    tole_lodir <- which(is.na(rawfile$data[,londir]))
    cat("\n   -   ", length(tole_lodir), " NAs found in longitude direction...", sep = "")
    #to_rep <- paste(to_rep, "\n   -   ", length(tole_lodir), " NAs found in longitude direction...\n", sep = "")
    if(length(tole_lodir) > 0)
    {
      rawfile$data <- rawfile$data[-tole_lodir,]
      vmsdata <- vmsdata[-tole_lodir,]
    }
    
    pre_lon <- lonsex2dec(rawfile$data[,londeg],rawfile$data[,lonmin],rawfile$data[,lonsec],rawfile$data[,londir])
    
    tole_prelon <- which(pre_lon < -180 | pre_lon > 180)
    if(length(tole_prelon) > 0)
    {
      pre_lon <- pre_lon[-tole_prelon]
      vmsdata <- vmsdata[-tole_prelon,]
      rawfile$data <- rawfile$data[-tole_prelon,]
    }
    cat("\n   -   ", length(tole_prelon), " Longitudes out of range ( -180 / 180 )", sep = "")
    
    vmsdata[,"LON"] <- pre_lon
    
  }
  
  if (lonModeSel == "dec")
  {
    
    londec <- which(colnames(rawfile$data) == lonDec)
    tole_lode <- which(is.na(rawfile$data[,londec]))
    cat("\n   -   ", length(tole_lode), " NAs found in decimal longitude...", sep = "")
    if(length(tole_lode) > 0)
    {
      rawfile$data <- rawfile$data[-tole_lode,]
      vmsdata <- vmsdata[-tole_lode,]
    }
    
    pre_lon <- rawfile$data[,londec]
    
    tole_prelon <- which(pre_lon < -180 | pre_lon > 180)
    if(length(tole_prelon) > 0)
    {
      pre_lon <- pre_lon[-tole_prelon]
      vmsdata <- vmsdata[-tole_prelon,]
      rawfile$data <- rawfile$data[-tole_prelon,]
    }
    cat("\n   -   ", length(tole_prelon), " Longitudes out of range ( -180 / 180 )", sep = "")
    
    vmsdata[,"LON"] <- pre_lon
    
  }
  
  
  if (timeMode == "UTC")
  {
    
    timeutc <- which(colnames(rawfile$data) == timeUtc)
    tole_ti <- which(is.na(rawfile$data[,timeutc]))
    cat("\n   -   ", length(tole_ti), " NAs found in UTC times...", sep = "")
    #to_rep <- paste(to_rep, "\n   -   ", length(tole_ti), " NAs found in UTC times...\n", sep = "")
    if(length(tole_ti) > 0)
    {
      rawfile$data <- rawfile$data[-tole_ti,]
      vmsdata <- vmsdata[-tole_ti,]
    }
    vmsdata[,"DATE"] <- rawfile$data[,timeutc]
    
  }
  
  if (timeMode == "Date + Time")
  {
    
    date <- which(colnames(rawfile$data) == timeDate2a)
    time <- which(colnames(rawfile$data) == timeDate2b)
    tole_tida <- which(is.na(rawfile$data[,time]))
    cat("\n   -   ", length(tole_tida), " NAs found in times...", sep = "")
    #to_rep <- paste(to_rep, "\n   -   ", length(tole_tida), " NAs found in times...\n", sep = "")
    if(length(tole_tida) > 0)
    {
      rawfile$data <- rawfile$data[-tole_tida,]
      vmsdata <- vmsdata[-tole_tida,]
    }
    tole_dati <- which(is.na(rawfile$data[,date]))
    cat("\n   -   ", length(tole_dati), " NAs found in dates...", sep = "")
    #to_rep <- paste(to_rep, "\n   -   ", length(tole_dati), " NAs found in dates...\n", sep = "")
    if(length(tole_dati) > 0)
    {
      rawfile$data <- rawfile$data[-tole_dati,]
      vmsdata <- vmsdata[-tole_dati,]
    }
    vmsdata[,"DATE"] <- as.numeric(chron(as.character(rawfile$data[,date]),
                                        as.character(rawfile$data[,time]),
                                        data_frm))    
  }
  
  if (timeMode == "Date + H M S")
  {
    date <- which(colnames(rawfile$data) == timeDate)
    hour <- which(colnames(rawfile$data) == timeHour)
    minute <- which(colnames(rawfile$data) == timeMinute)
    second <- which(colnames(rawfile$data) == timeSecond)
    tole_dat <- which(is.na(rawfile$data[,date]))
    cat("\n   -   ", length(tole_dat), " NAs found in dates...", sep = "")
    #to_rep <- paste(to_rep, "\n   -   ", length(tole_dat), " NAs found in dates...\n", sep = "")
    if(length(tole_dat) > 0)
    {
      rawfile$data <- rawfile$data[-tole_dat,]
      vmsdata <- vmsdata[-tole_dat,]
    }
    tole_hou <- which(is.na(rawfile$data[,hour]))
    cat("\n   -   ", length(tole_hou), " NAs found in hours...", sep = "")
    #to_rep <- paste(to_rep, "\n   -   ", length(tole_hou), " NAs found in hours...\n", sep = "")
    if(length(tole_hou) > 0)
    {
      rawfile$data <- rawfile$data[-tole_hou,]
      vmsdata <- vmsdata[-tole_hou,]
    }
    tole_min <- which(is.na(rawfile$data[,minute]))
    cat("\n   -   ", length(tole_min), " NAs found in minutes...", sep = "")
    #to_rep <- paste(to_rep, "\n   -   ", length(tole_min), " NAs found in minutes...\n", sep = "")
    if(length(tole_min) > 0)
    {
      rawfile$data <- rawfile$data[-tole_min,]
      vmsdata <- vmsdata[-tole_min,]
    }
    tole_sec <- which(is.na(rawfile$data[,second]))
    cat("\n   -   ", length(tole_sec), " NAs found in seconds...", sep = "")
    #to_rep <- paste(to_rep, "\n   -   ", length(tole_sec), " NAs found in seconds...\n", sep = "")
    if(length(tole_sec) > 0)
    {
      rawfile$data <- rawfile$data[-tole_sec,]
      vmsdata <- vmsdata[-tole_sec,]
    }
    time <- paste(rawfile$data[,hour], rawfile$data[,minute], rawfile$data[,second], sep = ":")
    vmsdata[,"DATE"] <- as.numeric(chron(as.character(rawfile$data[,date]),
                                        as.character(time),
                                        data_frm))
  }
  
  tole_vda <- which(is.na(vmsdata[,"DATE"]))
  cat("\n   -   ", length(tole_vda), " dates found with bad format...", sep = "")
  #to_rep <- paste(to_rep, "\n   -   ", length(tole_vda), " dates found with bad format...\n", sep = "")
  if(length(tole_vda) > 0)
  {
    rawfile$data <- rawfile$data[-tole_vda,]
    vmsdata <- vmsdata[-tole_vda,]
  }
  
  if (speedModeSel == "Knots")
  {
    
    knots <- which(colnames(rawfile$data) == speedCol)
    tole_kn <- which(is.na(rawfile$data[,knots]))
    cat("\n   -   ", length(tole_kn), " NAs found in knots speed...", sep = "")
    #to_rep <- paste(to_rep, "\n   -   ", length(tole_kn), " NAs found in knots speed...\n", sep = "")
    if(length(tole_kn) > 0)
    {
      rawfile$data <- rawfile$data[-tole_kn,]
      vmsdata <- vmsdata[-tole_kn,]
    }
    vmsdata[,"SPE"] <- kno2kmh(rawfile$data[,knots])
    
  }
  
  if (speedModeSel == "Km/h")
  {
    
    kmh <- which(colnames(rawfile$data) == speedCol)
    tole_km <- which(is.na(rawfile$data[,kmh]))
    cat("\n   -   ", length(tole_km), " NAs found in km/h speed...", sep = "")
    #to_rep <- paste(to_rep, "\n   -   ", length(tole_km), " NAs found in km/h speed...\n", sep = "")
    if(length(tole_km) > 0)
    {
      rawfile$data <- rawfile$data[-tole_km,]
      vmsdata <- vmsdata[-tole_km,]
    }
    vmsdata[,"SPE"] <- rawfile$data[,kmh]
    
  }
  
  
  if (headModeSel == "Rad")
  {
    
    rad <- which(colnames(rawfile$data) == headCol)
    tole_ra <- which(is.na(rawfile$data[,rad]))
    cat("\n   -   ", length(tole_ra), " NAs found in radiants heading...", sep = "")
    #to_rep <- paste(to_rep, "\n   -   ", length(tole_ra), " NAs found in radiants heading...\n", sep = "")
    if(length(tole_ra) > 0)
    {
      rawfile$data <- rawfile$data[-tole_ra,]
      vmsdata <- vmsdata[-tole_ra,]
    }
    vmsdata[,"HEA"] <- rad2deg(rawfile$data[,rad])
    
  }
  
  if (headModeSel == "Deg")
  {
    
    deg <- which(colnames(rawfile$data) == headCol)
    tole_de <- which(is.na(rawfile$data[,deg]))
    cat("\n   -   ", length(tole_de), " NAs found in degrees heading...", sep = "")
    #to_rep <- paste(to_rep, "\n   -   ", length(tole_de), " NAs found in degrees heading...\n", sep = "")
    if(length(tole_de) > 0)
    {
      rawfile$data <- rawfile$data[-tole_de,]
      vmsdata <- vmsdata[-tole_de,]
    }
    vmsdata[,"HEA"] <- rawfile$data[,deg]
    
  }
  cat("\nRemoved ", round((100/numlines) * (numlines-nrow(vmsdata)), 2), "% of data, that is ", numlines-nrow(vmsdata)," pings\n",
      "\n\n   ---   Raw Pings Editing Complete!   ---\n\n", sep = "")
  
  return(vmsdata)
}