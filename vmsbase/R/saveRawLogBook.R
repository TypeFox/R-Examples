
#' Raw LogBook editing function
#' 
#' \code{saveRawLogBook} implements the routines that converts raw values
#'  to standard data.
#' 
#' @param rawfile    The raw LogBook dataset.
#' @param widget    The widget list that contains the editing infos.
#'
#' @return The function returns the standardized VMS data.
#' 
#' @usage saveRawLogBook(rawfile, widget)
#' 
#' @export saveRawLogBook
#'@seealso \code{\link{gui_lb_editraw}}

saveRawLogBook <- function (rawfile,
                            widget)
{
  vess <- widget[1]
  th_met_s <- widget[2]
  th_meti <- widget[3]
  th_gea_s <- widget[4]
  th_gear <- widget[5]
  #   sta_sel <- widget[2]
  sdate <- widget[6]
  stime <- widget[7]
  #   end_sel <- widget[5]
  edate <- widget[8]
  etime <- widget[9]
  species <- widget[10]
  qty <- widget[11]
  if(widget[12] == "DD/MM/YYYY")
  {data_frm <-  c(dates = "d/m/y", times = "h:m:s")
  }else{
    data_frm <- c(dates = "m/d/y", times = "h:m:s")
  }
  
  numlines <- length(rawfile$data[,1])
  
  cat("\n\n   ---   Raw Logbooks Editing Started!   ---\n",
      "\nProcessing ", numlines, " raw logbooks...\n", sep = "")
  
  logbook <- data.frame("vessUE" = numeric(numlines),
                        "s_utc" = numeric(numlines),
                        "e_utc" = numeric(numlines),
                        "specie" = character(numlines),
                        "qty" = numeric(numlines),
                        "gear" = numeric(numlines),
                        "metier" = numeric(numlines))
  
  
  vess_id <- which(colnames(rawfile$data) == vess)
  logbook["vessUE"] <- rawfile$data[,vess_id]
  
  #   if(sta_sel == "Yes")
  #   {
  log_sdate <- which(colnames(rawfile$data) == sdate)
  log_stime <- which(colnames(rawfile$data) == stime)
  sdts <- gsub("\"", "", as.character(rawfile$data[,log_sdate]))
  stms <- gsub("\"", "", as.character(rawfile$data[,log_stime]))
  toes <- which(nchar(stms) == 5)
  if(length(toes) > 0){stms[toes] <- paste(stms[toes], ":00", sep = "")}
  sta_utc <- as.numeric(chron(dates. = sdts, times. = stms, format = data_frm))
  
  tole_sutc <- which(is.na(sta_utc))
  cat("\n   -   ", length(tole_sutc), " NAs found in Start Times...", sep = "")
  if(length(tole_sutc) > 0)
  {
    rawfile$data <- rawfile$data[-tole_sutc,]
    sta_utc <- sta_utc[-tole_sutc]
    logbook <- logbook[-tole_sutc,]
  }
  
  logbook["s_utc"] <- sta_utc
  #   }else{logbook["s_utc"] <- 0}
  
  #   if(end_sel == "Yes")
  #   {
  log_edate <- which(colnames(rawfile$data) == edate)
  log_etime <- which(colnames(rawfile$data) == etime)
  edts <- gsub("\"", "", as.character(rawfile$data[,log_edate]))
  etms <- gsub("\"", "", as.character(rawfile$data[,log_etime]))
  toee <- which(nchar(etms) == 5)
  if(length(toee) > 0){etms[toee] <- paste(etms[toee], ":00", sep = "")}
  end_utc <- as.numeric(chron(dates. = edts, times. = etms, format = data_frm))
  
  tole_eutc <- which(is.na(end_utc))
  cat("\n   -   ", length(tole_eutc), " NAs found in End Times...", sep = "")
  if(length(tole_eutc) > 0)
  {
    rawfile$data <- rawfile$data[-tole_eutc,]
    end_utc <- end_utc[-tole_eutc]
    logbook <- logbook[-tole_eutc,]
  }
  
  logbook["e_utc"] <- end_utc
  #   }else{logbook["e_utc"] <- 0}
  
  
  spcs <- which(colnames(rawfile$data) == species)
  
  tole_spc <- which(is.na(rawfile$data[,spcs]))
  cat("\n   -   ", length(tole_spc), " NAs found in Species...", sep = "")
  if(length(tole_spc) > 0)
  {
    rawfile$data <- rawfile$data[-tole_spc,]
    logbook <- logbook[-tole_spc,]
  }
  
  logbook["specie"] <- gsub("\"", "", rawfile$data[,spcs])
  
  sp_qty <- which(colnames(rawfile$data) == qty)
  
  tole_qty <- which(is.na(rawfile$data[,sp_qty]))
  cat("\n   -   ", length(tole_qty), " NAs found in Quantity...", sep = "")
  if(length(tole_qty) > 0)
  {
    rawfile$data <- rawfile$data[-tole_qty,]
    logbook <- logbook[-tole_qty,]
  }
  
  logbook["qty"] <- as.numeric(rawfile$data[,sp_qty])
  
  
  if(th_met_s == "Yes")
  {
    log_meti <- which(colnames(rawfile$data) == th_meti)
    met_inf <- gsub("\"", "", as.character(rawfile$data[,log_meti]))
    tole_met <- which(is.na(met_inf))
    cat("\n   -   ", length(tole_met), " NAs found in Metier...", sep = "")
    if(length(tole_met) > 0)
    {
      rawfile$data <- rawfile$data[-tole_met,]
      met_inf <- met_inf[-tole_met]
      logbook <- logbook[-tole_met,]
    }
    logbook["metier"] <- met_inf
  }else{logbook["metier"] <- 0}
  
  if(th_gea_s == "Yes")
  {
    log_gear <- which(colnames(rawfile$data) == th_gear)
    gea_inf <- gsub("\"", "", as.character(rawfile$data[,log_gear]))
    tole_gea <- which(is.na(gea_inf))
    cat("\n   -   ", length(tole_gea), " NAs found in Gear...", sep = "")
    if(length(tole_gea) > 0)
    {
      rawfile$data <- rawfile$data[-tole_gea,]
      gea_inf <- gea_inf[-tole_gea]
      logbook <- logbook[-tole_gea,]
    }
    logbook["gear"] <- gea_inf
  }else{logbook["gear"] <- 0}
  
  
  cat("\n\nRemoved ",
      round((100/numlines) * (numlines-nrow(logbook)), 2), "% of data, that is ",
      numlines-nrow(logbook)," logbooks\n",
      "\n\n   ---   Raw Logbooks Editing Complete!   ---\n\n", sep = "")
  
  
  return(logbook)
}




saveRawEflalo <- function (rawfile)
{
  numlines <- length(rawfile$data[,1])
  
  logbook <- data.frame("vessUE" = numeric(numlines),
                        "s_utc" = numeric(numlines),
                        "e_utc" = numeric(numlines),
                        "gear" = character(numlines),
                        "metier" = character(numlines))
  
  vess_id <- which(colnames(rawfile$data) == "VE_REF")
  logbook["vessUE"] <- rawfile$data[,vess_id]
  data_frm <-  c(dates = "d/m/y", times = "h:m:s")
  log_sdate <- which(colnames(rawfile$data) == "FT_DDAT")
  log_stime <- which(colnames(rawfile$data) == "FT_DTIME")
  sdts <- gsub("\"", "", as.character(rawfile$data[,log_sdate]))
  stms <- gsub("\"", "", as.character(rawfile$data[,log_stime]))
  #toes <- which(nchar(stms) == 5)
  #if(length(toes) > 0){stms[toes] <- paste(stms[toes], ":00", sep = "")}
  sta_utc <- as.numeric(chron(dates. = sdts, times. = stms, format = data_frm))
  
  tole_sutc <- which(is.na(sta_utc))
  cat("\n   -   ", length(tole_sutc), " NAs found in Start Times...", sep = "")
  if(length(tole_sutc) > 0)
  {
    rawfile$data <- rawfile$data[-tole_sutc,]
    sta_utc <- sta_utc[-tole_sutc]
    logbook <- logbook[-tole_sutc,]
  }
  
  logbook["s_utc"] <- sta_utc
  
  log_edate <- which(colnames(rawfile$data) == "FT_LDAT")
  log_etime <- which(colnames(rawfile$data) == "FT_LTIME")
  edts <- gsub("\"", "", as.character(rawfile$data[,log_edate]))
  etms <- gsub("\"", "", as.character(rawfile$data[,log_etime]))
  end_utc <- as.numeric(chron(dates. = edts, times. = etms, format = data_frm))
  
  tole_eutc <- which(is.na(end_utc))
  cat("\n   -   ", length(tole_eutc), " NAs found in End Times...", sep = "")
  if(length(tole_eutc) > 0)
  {
    rawfile$data <- rawfile$data[-tole_eutc,]
    end_utc <- sta_utc[-tole_eutc]
    logbook <- logbook[-tole_eutc,]
  }
  
  logbook["e_utc"] <- end_utc
  
  metco <- grep("LE_MET", colnames(rawfile$data))  
  logbook[,"metier"] <- as.character(rawfile$data[,metco])
  
  geaco <- grep("LE_GEAR", colnames(rawfile$data))  
  logbook[,"gear"] <- as.character(rawfile$data[,geaco])
  
  speco <- grep("LE_KG", colnames(rawfile$data))  
  spequ <- rawfile$data[,speco]
  
  logbook <- cbind(logbook, spequ)
  colnames(logbook) <- gsub("LE_KG", "FAO", colnames(logbook))
  
  cat("\nRemoved ", round((100/numlines) * (numlines-nrow(logbook)), 2), "% of data, that is ", numlines-nrow(logbook)," logbooks\n",
      "\n\n   ---   Raw Logbooks Editing Complete!   ---\n\n", sep = "")
  
  
  return(logbook)
}
