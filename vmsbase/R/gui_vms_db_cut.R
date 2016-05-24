
#' VMS DB Cut GUI
#'  
#' The \code{gui_vms_db_cut} function implements the graphical user interface for the
#'  VMS Track cutting routine.
#' 
#' This function, with a VMS cleaned database (see \code{\link{gui_vms_db_clean}}),
#'  assign track numbers to the raw pings of each vessel in the VMS DB.
#'  
#' @param vms_db_name The path of a VMS DataBase
#' 
#' @return This function does not return a value. 
#' 
#' @usage gui_vms_db_cut(vms_db_name = "")
#' 
#' @export gui_vms_db_cut
#'
#'@seealso \code{\link{gui_vms_db_clean}}

gui_vms_db_cut <- function(vms_db_name = "")
{
  vms_DB <- vms_DB$new()
  vms_DB$db <- vms_db_name
  
  vms_db_cut_win <- gwindow("VMS Track Cutter Utility", visible = FALSE)
  
  # TRACK CUTTER
  cut_g <- gframe(horizontal = FALSE, container = vms_db_cut_win)
  cut_g2 <- ggroup(horizontal = TRUE, container = cut_g)
  addSpring(cut_g2)
  gimage(system.file("ico/retroshare1.png", package="vmsbase"), container = cut_g2)
  proglab_cut <- glabel("Track Cutter" , container = cut_g2)
  addSpring(cut_g2)
  
  cut_g3 <- ggroup(horizontal = TRUE, container = cut_g)
  
  #################
  addSpring(cut_g3)
  vms_db_f <- gframe(text = "VMS DB file", horizontal = TRUE, container = cut_g3)
  addSpring(vms_db_f)
  sel_vms_f <- glabel("Select VMS DB file", container = vms_db_f)
  addSpring(vms_db_f)
  gimage(system.file("ico/folder-blue.png", package="vmsbase"), container = vms_db_f,
         handler = function(h,...){
           vms_DB$db <- gfile(text = "Select VMS DataBase file",
                              type = "open",
                              filter = list("VMS DB file" = list(patterns = c("*.vms.sqlite"))))
           cat(vms_DB$db)
           svalue(sel_vms_f) <- ifelse(.Platform$OS.type == "windows", strsplit(vms_DB$db, "\\\\")[[1]][length(strsplit(vms_DB$db, "\\\\")[[1]])],strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])])
           enabled(start_b) <- TRUE
         })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = vms_db_f,
         handler = function(h,...){
           vms_DB$db <- ""
           enabled(start_b) <- FALSE
           svalue(sel_vms_f) <- "Select VMS DB file"
         })
  addSpring(cut_g3)
  ################
  
  cut_opt <- ggroup(horizontal = TRUE, container = cut_g)
  
  addSpring(cut_opt)
  opt_f <- gframe(text = "Split Options", horizontal = FALSE, container = cut_opt, expand = TRUE)
  addSpace(opt_f, 20, horizontal = TRUE)
  opt_out_g <- ggroup(horizontal = TRUE, container = opt_f)
  addSpring(opt_out_g)
  opt_out_c <- gcheckbox("Use Outliers", checked = TRUE, container = opt_out_g)
  addSpring(opt_out_g)
  addSpace(opt_f, 20, horizontal = TRUE)
  opt_med_g <- ggroup(horizontal = TRUE, container = opt_f)
  addSpring(opt_med_g)
  opt_med_c <- gcheckbox("Lag > Median", checked = TRUE, container = opt_med_g)
  addSpring(opt_med_g)
  addSpace(opt_f, 20, horizontal = FALSE)
  glabel(text = "Lag > Minutes", container = opt_f)
  opt_hou_s <- gslider(from = 30, to = 360, by = 5, 
                       value = 120, horizontal = TRUE, container = opt_f)
  addSpace(opt_f, 20, horizontal = TRUE)
  addSpring(cut_opt)
  
  ################
  addSpring(cut_g)
  infolab_cut <- glabel("" , container = cut_g)
  addSpring(cut_g)
  start_b <- gbutton("Start cutting", container = cut_g, handler = function(h,...)
  {
    enabled(vms_db_f) <- FALSE
    enabled(start_b) <- FALSE
    
    if(sqldf("select count(*) from warn", dbname = vms_DB$db) > 0)
    {
      svalue(infolab_cut) <- "Updating...\nVMS DataBase"
      sqldf("drop table if exists track", dbname = vms_DB$db)
      sqldf("CREATE TABLE track(I_NCEE INT, LAT REAL, LON REAL, DATE REAL, SPE REAL, HEA REAL, W_HARB INT, T_NUM INT, P_ID INT)", dbname = vms_DB$db)
      cat("\n\n   ---   Track Cutting Started   ---\n")
      incee <- sqldf("select distinct I_NCEE from ping", dbname = vms_DB$db)
      num_incee <- nrow(incee)
      for ( v in 1:num_incee )
      {
        svalue(infolab_cut) <- paste("Processing...\nVessel: ", v," of ", num_incee, spe = "")
        cat("\nVessel: ", v," of ", num_incee, spe = "")
        vessel <- fn$sqldf("select I_NCEE, LAT, LON, DATE, SPE, HEA, warn.* from ping, warn where ping.ROWID = warn.p_id and I_NCEE = `incee[v,1]` and W_DUPL = 0 and W_COHE != 0 and W_LAND = 0 order by DATE ", dbname = vms_DB$db)
        if(nrow(vessel) == 0)
        {
          cat(" - Skipped, no pings", sep = "")
        }else{
          numlines <- nrow(vessel)
          cat(" with ", numlines, " pings ", sep = "")
          track_data <- data.frame("I_NCEE" = numeric(numlines),
                                   "LAT" = numeric(numlines),
                                   "LON" = numeric(numlines),
                                   "DATE" = numeric(numlines),
                                   "SPE" = numeric(numlines),
                                   "HEA" = numeric(numlines),
                                   "W_HARB" = integer(numlines),
                                   "T_NUM" = numeric(numlines),
                                   "ROWID" = numeric(numlines))
          
          track_data["ROWID"] <- vessel["p_id"]
          track_data["I_NCEE"] <- vessel["I_NCEE"]
          track_data["LAT"] <- vessel["LAT"]
          track_data["LON"] <- vessel["LON"]
          track_data["DATE"] <- vessel["DATE"]
          track_data["SPE"] <- vessel["SPE"]
          track_data["HEA"] <- vessel["HEA"]
          track_data["W_HARB"] <- vessel["W_HARB"]
          spe3qua = quantile(as.numeric(vessel[, "SPE"]), probs = 0.75)
          useharb <- unique(cbind(vessel[which(vessel$W_HARB == 1),"LON"], 
                                  vessel[which(vessel$W_HARB == 1),"LAT"]))
          #assign raw track numbers and add in harbour points
          tr_num <- 1
          for (j in 1:numlines)
          {
            tr_lin <- nrow(track_data)
            delta <- tr_lin - numlines
            if(vessel[j,"W_HARB"] == 1 | vessel[j,"W_LAND"] == 1)
            {
              track_data[j+delta,"T_NUM"] <- 0
              next
            }else{
              track_data[j+delta,"T_NUM"] <- tr_num
              if(nrow(useharb) == 0)
              {next}
              if(nrow(useharb) > 1)
              {
                dist <- spDistsN1(useharb[,1:2], 
                                  as.matrix(c(vessel[j,"LON"],
                                              vessel[j,"LAT"])),
                                  longlat = TRUE)
              }
              if(nrow(useharb) == 1)
              {
                dist <- spDists(cbind(useharb[,1],useharb[,2]), 
                                (cbind(vessel[j,"LON"], vessel[j,"LAT"])),
                                longlat = TRUE)
              }
              hdist <- min(dist)
              nearh <- which(dist == hdist)
              if(j > 1)
              {
                if(vessel[j-1,"W_HARB"] == 1)
                {
                  newdate <- vessel[j,"DATE"]-((1/24)*(hdist/spe3qua))
                  track_data <- rbind(track_data[1:(j+delta-1),],
                                      c(vessel[j,"I_NCEE"], 
                                        useharb[nearh[1],2], 
                                        useharb[nearh[1],1], 
                                        ifelse(vessel[j-1,"DATE"] < newdate & newdate != Inf & newdate != -Inf, newdate, vessel[j-1,"DATE"]),
                                        0,
                                        0,
                                        1,
                                        tr_num), 
                                      track_data[(j+delta):nrow(track_data),])
                }
              }
              tr_lin <- nrow(track_data)
              delta <- tr_lin - numlines
              if(j < nrow(vessel))
              {
                if(vessel[j+1,"W_HARB"] == 1)
                {
                  newdate <- vessel[j,"DATE"]+((1/24)*(hdist/spe3qua))
                  track_data <- rbind(track_data[1:(j+delta),], 
                                      c(vessel[j,"I_NCEE"],
                                        useharb[nearh[1],2],
                                        useharb[nearh[1],1],
                                        ifelse(vessel[j+1,"DATE"] > newdate & newdate != Inf & newdate != -Inf, newdate, vessel[j+1,"DATE"]),
                                        0, 0, 1, tr_num), 
                                      track_data[(j+delta+1):nrow(track_data),])
                  tr_num <- tr_num + 1
                }
              }
            }
          }
          track_data <- track_data[which(track_data$T_NUM > 0),]
          if(nrow(track_data) > 0)
          {
            # pingfreq <- diff(track_data$DATE[which(track_data$W_HARB == 0)])
            pingfreq <- diff(track_data$DATE)
            
            if(svalue(opt_out_c) == TRUE){
              
              out = boxplot(diff(track_data$DATE[which(track_data$W_HARB == 0)]),plot=F)$out
              
              if(length(out) > 0)
              {
                if(svalue(opt_med_c) == TRUE)
                {
                  da_to <- which(out > median(pingfreq) & out > (1/24)*(svalue(opt_hou_s)/60))
                }else{
                  da_to <- which(out > (1/24)*(svalue(opt_hou_s)/60)) 
                }
                
                if(length(da_to) > 0)
                {
                  outmin = min(out[da_to])
                  outliers = which((diff(track_data$DATE) >= outmin | diff(track_data$DATE) >= 0.125) & track_data$W_HARB[-c(nrow(track_data))] != 1 )
                  if(length(outliers) > 0)
                  {
                    cat(" - Splitting Tracks ", sep = "")
                    for (k in 1:length(outliers))
                    {    
                      cat(":", sep = "")
                      if(nrow(useharb) == 0)
                      {
                        newdate1 <- track_data[(outliers[k]),"DATE"]+0.0002
                        newdate2 <- newdate1
                        track_data <- rbind(track_data[1:(outliers[k]), ], 
                                            c(track_data[(outliers[k]),"I_NCEE"], track_data[(outliers[k]),"LAT"], track_data[(outliers[k]),"LON"], newdate1, 0, 0, 0, (track_data[(outliers[k]), "T_NUM"])),
                                            c(track_data[(outliers[k]+1),"I_NCEE"], track_data[(outliers[k]+1),"LAT"], track_data[(outliers[k]+1),"LON"], newdate2, 0, 0, 0, (track_data[(outliers[k]+1), "T_NUM"])),
                                            track_data[(outliers[k]+1):nrow(track_data),])
                        track_data$T_NUM[(outliers[k]+2):nrow(track_data)] <- track_data$T_NUM[(outliers[k]+2):nrow(track_data)]+1
                        outliers <- outliers + 2
                        next
                      }
                      if(nrow(useharb) > 1)
                      {
                        dist1 <-spDistsN1(useharb[,1:2], 
                                          as.matrix(c(track_data[outliers[k],"LON"], track_data[outliers[k],"LAT"])), 
                                          longlat = TRUE)
                      }
                      
                      if(nrow(useharb) == 1)
                      {
                        dist1 <- spDists(cbind(useharb[,1],useharb[,2]), 
                                         (cbind(track_data[outliers[k],"LON"], track_data[outliers[k],"LAT"])),
                                         longlat = TRUE)
                      }
                      hdist1 <- min(dist1)
                      nearh1 <- which(dist1 == hdist1)
                      if(nrow(useharb) > 1)
                      {
                        dist2 <-spDistsN1(useharb[,1:2],
                                          as.matrix(c(track_data[(outliers[k]+1),"LON"], track_data[(outliers[k]+1),"LAT"])),
                                          longlat = TRUE)
                      }
                      if(nrow(useharb) == 1)
                      {
                        dist2 <- spDists(cbind(useharb[,1],useharb[,2]), 
                                         (cbind(track_data[(outliers[k]+1),"LON"], track_data[(outliers[k]+1),"LAT"])),
                                         longlat = TRUE)
                      }
                      hdist2 <- min(dist2)
                      nearh2 <- which(dist2 == hdist2)
                      newdate1 <- track_data[(outliers[k]),"DATE"]+((1/24)*(spe3qua/hdist1))
                      newdate1 <- ifelse(track_data[(outliers[k]+1),"DATE"] > newdate1 & newdate1 != Inf & newdate1 != -Inf, 
                                         newdate1,
                                         track_data[(outliers[k]),"DATE"]+0.0002)
                      newdate2 <- track_data[(outliers[k]+1),"DATE"]-((1/24)*(spe3qua/hdist2))
                      newdate2 <- ifelse(newdate1 < newdate2 & newdate2 != Inf & newdate2 != -Inf,
                                         newdate2, 
                                         newdate1)
                      track_data <- rbind(track_data[1:(outliers[k]), ], 
                                          c(track_data[(outliers[k]),"I_NCEE"], useharb[nearh[1],2], useharb[nearh[1],1], newdate1, 0, 0, 1, (track_data[(outliers[k]), "T_NUM"]), NA),
                                          c(track_data[(outliers[k]+1),"I_NCEE"], useharb[nearh2[1],2], useharb[nearh2[1],1], newdate2, 0, 0, 1, (track_data[(outliers[k]+1), "T_NUM"]), NA),
                                          track_data[(outliers[k]+1):nrow(track_data),])
                      track_data$T_NUM[(outliers[k]+2):nrow(track_data)] <- track_data$T_NUM[(outliers[k]+2):nrow(track_data)]+1
                      outliers <- outliers + 2
                    }
                  }
                }else{
                  cat(" -  No Ping Frequency outliers found ", sep = "")
                }
              }else{
                cat(" -  No Ping Frequency outliers found ", sep = "")
              }
            }else{
              if(svalue(opt_med_c) == TRUE)
              {
                outliers <- which(pingfreq > median(pingfreq) | pingfreq > (1/24)*(svalue(opt_hou_s)/60))
              }else{
                outliers <- which(pingfreq > (1/24)*(svalue(opt_hou_s)/60)) 
              }
              out_harb <- which(track_data$W_HARB[outliers] == 1)
              if(length(out_harb) > 0){
                outliers <- outliers[-out_harb]
              }
              
              if(length(outliers) > 0)
              {
                cat(" - Splitting Tracks ", sep = "")
                for (k in 1:length(outliers))
                {    
                  cat(":", sep = "")
                  if(nrow(useharb) == 0)
                  {
                    newdate1 <- track_data[(outliers[k]),"DATE"]+0.0002
                    newdate2 <- newdate1
                    track_data <- rbind(track_data[1:(outliers[k]), ], 
                                        c(track_data[(outliers[k]),"I_NCEE"], track_data[(outliers[k]),"LAT"], track_data[(outliers[k]),"LON"], newdate1, 0, 0, 0, (track_data[(outliers[k]), "T_NUM"])),
                                        c(track_data[(outliers[k]+1),"I_NCEE"], track_data[(outliers[k]+1),"LAT"], track_data[(outliers[k]+1),"LON"], newdate2, 0, 0, 0, (track_data[(outliers[k]+1), "T_NUM"])),
                                        track_data[(outliers[k]+1):nrow(track_data),])
                    track_data$T_NUM[(outliers[k]+2):nrow(track_data)] <- track_data$T_NUM[(outliers[k]+2):nrow(track_data)]+1
                    outliers <- outliers + 2
                    next
                  }
                  if(nrow(useharb) > 1)
                  {
                    dist1 <-spDistsN1(useharb[,1:2], 
                                      as.matrix(c(track_data[outliers[k],"LON"], track_data[outliers[k],"LAT"])), 
                                      longlat = TRUE)
                  }
                  
                  if(nrow(useharb) == 1)
                  {
                    dist1 <- spDists(cbind(useharb[,1],useharb[,2]), 
                                     (cbind(track_data[outliers[k],"LON"], track_data[outliers[k],"LAT"])),
                                     longlat = TRUE)
                  }
                  hdist1 <- min(dist1)
                  nearh1 <- which(dist1 == hdist1)
                  if(nrow(useharb) > 1)
                  {
                    dist2 <-spDistsN1(useharb[,1:2],
                                      as.matrix(c(track_data[(outliers[k]+1),"LON"], track_data[(outliers[k]+1),"LAT"])),
                                      longlat = TRUE)
                  }
                  if(nrow(useharb) == 1)
                  {
                    dist2 <- spDists(cbind(useharb[,1],useharb[,2]), 
                                     (cbind(track_data[(outliers[k]+1),"LON"], track_data[(outliers[k]+1),"LAT"])),
                                     longlat = TRUE)
                  }
                  hdist2 <- min(dist2)
                  nearh2 <- which(dist2 == hdist2)
                  newdate1 <- track_data[(outliers[k]),"DATE"]+((1/24)*(spe3qua/hdist1))
                  newdate1 <- ifelse(track_data[(outliers[k]+1),"DATE"] > newdate1 & newdate1 != Inf & newdate1 != -Inf, 
                                     newdate1,
                                     track_data[(outliers[k]),"DATE"]+0.0002)
                  newdate2 <- track_data[(outliers[k]+1),"DATE"]-((1/24)*(spe3qua/hdist2))
                  newdate2 <- ifelse(newdate1 < newdate2 & newdate2 != Inf & newdate2 != -Inf,
                                     newdate2, 
                                     newdate1)
                  track_data <- rbind(track_data[1:(outliers[k]), ], 
                                      c(track_data[(outliers[k]),"I_NCEE"], useharb[nearh[1],2], useharb[nearh[1],1], newdate1, 0, 0, 1, (track_data[(outliers[k]), "T_NUM"]), NA),
                                      c(track_data[(outliers[k]+1),"I_NCEE"], useharb[nearh2[1],2], useharb[nearh2[1],1], newdate2, 0, 0, 1, (track_data[(outliers[k]+1), "T_NUM"]), NA),
                                      track_data[(outliers[k]+1):nrow(track_data),])
                  
                  track_data$T_NUM[(outliers[k]+2):nrow(track_data)] <- track_data$T_NUM[(outliers[k]+2):nrow(track_data)]+1
                  outliers <- outliers + 2
                }
              }
            }
            sqldf("insert into track select * from `track_data`", dbname = vms_DB$db)
          }
        }
      }
      cat("\n\n   ---   End Track Cutting   ---\n", sep = "")
      gconfirm("VMS DB Track Cutting Completed!",
               title = "Confirm",
               icon = "info",
               parent = vms_db_cut_win,
               handler = function(h,...){dispose(vms_db_cut_win)})
    }else{
      gconfirm("Warning data not available\n\nExecute DB Cleaning first!",
               title = "Error",
               icon = "error",
               parent = vms_db_cut_win,
               handler = function(h,...){dispose(vms_db_cut_win)})
    } 
  })
  enabled(start_b) <- FALSE
  if(vms_DB$db != "")
  {
    svalue(sel_vms_f) <- ifelse(.Platform$OS.type == "windows", strsplit(vms_DB$db, "\\\\")[[1]][length(strsplit(vms_DB$db, "\\\\")[[1]])],strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])])
    enabled(start_b) <- TRUE
  } 
  visible(vms_db_cut_win) <- TRUE
}