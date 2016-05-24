
#' Mark Fishing Points GUI
#'  
#' 
#' The \code{gui_mark_fis_poi} function implements the graphical user interface for the
#'  Fishing Point Marking routine.
#' 
#' This function, with a VMS database and a shape file with harbours points, performs a filtered search over the whole
#'  db assigning fishing status to the vms interpolated data.
#'   
#' @param vms_db_name The path of a VMS DataBase
#' @param harb_file_name The path of a shape file with harbours point data
#' 
#' @return This function does not return a value. 
#' 
#' @usage gui_mark_fis_poi(vms_db_name = "", harb_file_name = "")
#' 
#' @export gui_mark_fis_poi  
#'

gui_mark_fis_poi <- function(vms_db_name = "", harb_file_name = "")
{
  
  vms_DB <- vms_DB$new()
  vms_DB$db <- vms_db_name
  harb <- harbCoo$new()
  harb$path <- harb_file_name
  
  met_list <- data.frame()
  
  mark_fis_poi_win <- gwindow(title = "Mark Fishing Points Tool",
                              visible = FALSE,
                              width = 950,
                              height = 350)
  
  big_g <- ggroup(horizontal = FALSE, container = mark_fis_poi_win, use.scrollwindow = TRUE)
  
  g_input <- ggroup(horizontal = TRUE, container = big_g)
  
  #### Load VMS DB
  
  addSpring(g_input)
  
  vms_db_f <- gframe(text = "VMS DB file", horizontal = TRUE, container = g_input)
  addSpring(vms_db_f)
  sel_vms_f <- glabel("Select VMS DB file", container = vms_db_f)
  addSpring(vms_db_f)
  
  gimage(system.file("ico/folder-blue.png", package="vmsbase"), container = vms_db_f,
         handler = function(h,...){
           vms_DB$db <- gfile(text = "Select VMS DataBase file",
                              type = "open",
                              filter = list("VMS DB file" = list(patterns = c("*.vms.sqlite"))))
           
           if(vms_DB$db != "")
           {
             svalue(sel_vms_f) <- ifelse(.Platform$OS.type == "windows", strsplit(vms_DB$db, "\\\\")[[1]][length(strsplit(vms_DB$db, "\\\\")[[1]])],strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])])
             
             #              nn_tab <- as.numeric(sqldf("SELECT count(*) FROM sqlite_master WHERE type='table' AND name='pre_nn'", dbname = vms_DB$db))
             #              if(nn_tab == 1)
             #              {
             #                nn_tab <- as.numeric(sqldf("SELECT count(*) FROM sqlite_master WHERE type='table' AND name='nn_clas'", dbname = vms_DB$db))
             #                if(nn_tab == 1)
             #                {
             enabled(gri_g3f4) <- TRUE
             #                }
             
             met_list <<- sqldf("select distinct met_des from vms_lb", dbname = vms_DB$db)
             if(nrow(met_list) == 0){
               nn_tab <- as.numeric(sqldf("SELECT count(*) FROM sqlite_master WHERE type='table' AND name='nn_:clas'", dbname = vms_DB$db))
               if(nn_tab == 1){
                 met_list <<- sqldf("select distinct met_des from nn_clas", dbname = vms_DB$db)
               }
             }
             if(nrow(met_list) > 0)
             {
               delete(big_g, up_g)
               up_g <<- gframe(text = "Filter Parameters", horizontal = TRUE, container = big_g)
               addSpring(up_g)
               lay_win <<- glayout(spacing = 5, container = up_g)
               lay_win[1,1, anchor = 0] <- "Metier"
               lay_win[1,2, anchor = 0] <- "Min Vel Kn"
               lay_win[1,3, anchor = 0] <- "Max Vel Kn"
               lay_win[1,4, anchor = 0] <- "Min Depth Mt"
               lay_win[1,5, anchor = 0] <- "Max Depth Mt"
               lay_win[1,6, anchor = 0] <- "Harb Dist Km"
               addSpring(up_g)
               for(i in 1:nrow(met_list))
               {
                 lay_win[i+1,1] <- as.character(met_list[i,1])
                 lay_win[i+1,2] <- gedit(text = "0", width = 10, container = lay_win)
                 lay_win[i+1,3] <- gedit(text = "Inf", width = 10, container = lay_win)
                 lay_win[i+1,4] <- gedit(text = "0", width = 10, container = lay_win)
                 lay_win[i+1,5] <- gedit(text = "-Inf", width = 10, container = lay_win)
                 lay_win[i+1,6] <- gedit(text = 3*1.85200, width = 10, container = lay_win)
               }
             }
             if(harb$path != "")
             {
               enabled(gri_alg) <- TRUE
               enabled(b_mark_fis_poi) <- TRUE
             }
           }
         })
  
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = vms_db_f,
         handler = function(h,...){
           vms_DB$db <- ""
           
           delete(big_g, up_g)
           
           #######################
           
           up_g <<- gframe(text = "Filter Parameters", horizontal = TRUE, container = big_g)
           #            add(big_g, up_g)
           
           #######################
           
           addSpring(up_g)
           lay_win <<- glayout(spacing = 5, container = up_g)
           lay_win[1,1, anchor = 0] <- "Metier"
           lay_win[1,2, anchor = 0] <- "Min Vel Kn"
           lay_win[1,3, anchor = 0] <- "Max Vel Kn"
           lay_win[1,4, anchor = 0] <- "Min Depth Mt"
           lay_win[1,5, anchor = 0] <- "Max Depth Mt"
           lay_win[1,6, anchor = 0] <- "Harb Dist Km"
           addSpring(up_g)
           
           #######################
           
           enabled(b_mark_fis_poi) <- FALSE
           enabled(gri_g3f4) <- FALSE
           enabled(gri_alg) <- FALSE
           svalue(sel_vms_f) <- "Select VMS DB file"
         })
  addSpring(g_input)
  
  #######################
  
  #### Load Harbours file
  
  cus_har_g <- gframe(text = "Harbours Shape File", horizontal = TRUE, container = g_input)
  addSpring(cus_har_g)
  cus_har_lab <- glabel("Select Harbours Shape File", container = cus_har_g)
  addSpring(cus_har_g)
  
  gimage(system.file("ico/folder-man.png", package="vmsbase"), container = cus_har_g,
         handler = function(h,...){
           harb$path <- gfile(text = "Select ShapePoints map",
                              type = "open",
                              filter = list("shp data" = list(patterns = c("*.shp"))))
           #            svalue(cus_har_lab) <- paste("Harbour: ", strsplit(harb, "/")[[1]][length(strsplit(harb, "/")[[1]])], sep = "")
           svalue(cus_har_lab) <- paste("Harbour: ", ifelse(.Platform$OS.type == "windows", strsplit(harb$path, "\\\\")[[1]][length(strsplit(harb$path, "\\\\")[[1]])],strsplit(harb$path, "/")[[1]][length(strsplit(harb$path, "/")[[1]])]), sep = "")
           
           if(harb$path != "" & vms_DB$db != "")
           {
             enabled(b_mark_fis_poi) <- TRUE
             enabled(gri_alg) <- TRUE
           }
         })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = cus_har_g,
         handler = function(h,...){
           harb$path <- ""
           svalue(cus_har_lab) <- "Select Harbours Shape File"
           enabled(b_mark_fis_poi) <- FALSE
           enabled(gri_alg) <- FALSE
         })
  addSpring(g_input)
  
  #######################
  
  #### Select Data Source
  addSpring(g_input)
  gri_g3f4 <- ggroup(horizontal = TRUE, container = g_input)
  addSpring(gri_g3f4)
  dat_sel_f <- gframe(text = "Metier Data Source", horizontal=TRUE, container = gri_g3f4) 
  met_da_su <- gdroplist(c("VMS-LB Match", "NN Prediction"), selected = 1, container = dat_sel_f,
                         handler = function(h,...){
                           
                           if(svalue(met_da_su) == "VMS-LB Match")
                           {
                             met_list <<- sqldf("select distinct met_des from vms_lb", dbname = vms_DB$db)
                           }else{
                             met_list <<- sqldf("select distinct met_des from nn_clas", dbname = vms_DB$db)
                           }
                           
                           if(nrow(met_list) > 0)
                           {
                             delete(big_g, up_g)
                             
                             #######################
                             
                             up_g <<- gframe(text = "Filter Parameters", horizontal = TRUE, container = big_g)
                             #                add(big_g, up_g)
                             
                             #######################
                             
                             addSpring(up_g)
                             lay_win <<- glayout(spacing = 5, container = up_g)
                             lay_win[1,1, anchor = 0] <- "Metier"
                             lay_win[1,2, anchor = 0] <- "Min Vel Kn"
                             lay_win[1,3, anchor = 0] <- "Max Vel Kn"
                             lay_win[1,4, anchor = 0] <- "Min Depth Mt"
                             lay_win[1,5, anchor = 0] <- "Max Depth Mt"
                             lay_win[1,6, anchor = 0] <- "Harb Dist Km"
                             addSpring(up_g)
                             
                             #######################
                             
                             for(i in 1:nrow(met_list))
                             {
                               
                               lay_win[i+1,1] <- as.character(met_list[i,1])
                               lay_win[i+1,2] <- gedit(text = "0", width = 10, container = lay_win)
                               lay_win[i+1,3] <- gedit(text = "Inf", width = 10, container = lay_win)
                               lay_win[i+1,4] <- gedit(text = "0", width = 10, container = lay_win)
                               lay_win[i+1,5] <- gedit(text = "-Inf", width = 10, container = lay_win)
                               lay_win[i+1,6] <- gedit(text = 3*1.85200, width = 10, container = lay_win)
                               
                             }
                           }
                         }) 
  addSpring(gri_g3f4)
  enabled(gri_g3f4) <- FALSE
  
  #### Select Algorithm
  addSpring(g_input)
  gri_alg <- ggroup(horizontal = TRUE, container = g_input)
  addSpring(gri_alg)
  alg_sel_f <- gframe(text = "Method", horizontal=TRUE, container = gri_alg) 
  sel_alg <- gdroplist(c("Slow & Ligth", "Fast & Heavy"), selected = 1, container = alg_sel_f) 
  addSpring(gri_alg)
  enabled(gri_alg) <- FALSE
  
  
  addSpring(g_input)
  b_mark_fis_poi <- gbutton(text = "\n   Start   \n", container = g_input, handler = function(h,...)
  {
    enabled(big_g) <- FALSE
    
    harbs <- readShapePoints(harb$path)
    
    sta_met <- cbind(met_list, matrix(0, ncol = 5, nrow = nrow(met_list)))
    colnames(sta_met) <- c("Metier", "min_vel_kn",  "max_vel_kn", "min_depth_mt", "max_depth_mt", "harb_dist_km")
    for(m in 1:nrow(sta_met))
    {
      sta_met[m,2:6] <- as.numeric(eval(parse(text = paste("c(", paste("svalue(lay_win[", m+1,",", 2:6, "])", sep = "", collapse = ", "), ")", sep = ""))))
    }
    
    cat("\n\n   ---   Fishing Point Analysis Started   ---\n", sep = "")
    
    sqldf("drop table if exists p_fish", dbname = vms_DB$db)
    sqldf("CREATE TABLE p_fish(i_id INT, F_SPE INT, F_DEP INT, F_DIS INT, FISH INT)", dbname = vms_DB$db)
    sqldf("drop table if exists p_fish_nn", dbname = vms_DB$db)
    sqldf("CREATE TABLE p_fish_nn(i_id INT, F_SPE INT, F_DEP INT, F_DIS INT, FISH INT)", dbname = vms_DB$db)
    
    if(svalue(sel_alg) == "Slow & Ligth")
    {
      incee_vms <- sqldf("select distinct I_NCEE from intrp", dbname = vms_DB$db)
      
      num_vess <- nrow(incee_vms)
      
      cat("\n   -     Metier Data Source: ", svalue(met_da_su),"     -\n", sep = "")
      
      for(i in 1:num_vess)
      {
        
        cat("\nVessel: ", incee_vms[i,1], " - N.", i, " of ", num_vess, sep = "")
        
        if(svalue(met_da_su) == "VMS-LB Match")
        {
          match <- fn$sqldf("select vessel, track, met_des from vms_lb where vessel = `incee_vms[i,1]`", dbname = vms_DB$db)
        }else{
          match <- fn$sqldf("select I_NCEE, T_NUM, met_des from nn_clas where I_NCEE = `incee_vms[i,1]`", dbname = vms_DB$db)
          colnames(match) <- c("vessel", "track", "met_des")
        } 
        
        num_track <- nrow(match)
        
        if(num_track == 0)
        {
          
          cat(" - Skipped, VMS-LogBook Match not found!")
          next
          
        }
        
        for(k in 1:num_track)
        {
          #cat("\nTrack: ", k, " of ", num_track, sep = "")
          cat(".", sep = "")
          eff_tra <- match[k,"track"]
          sin_tra <- fn$sqldf("select * from intrp, p_depth where I_NCEE = `incee_vms[i,1]` and T_NUM = `eff_tra` and intrp.rowid = i_id", dbname = vms_DB$db)
          
          if(nrow(sin_tra) == 0){
            cat("-", sep = "")
            next
          }
          
          tra_met <- match[k,"met_des"]
          sta_met_num <- which(sta_met[,"Metier"] == tra_met)
          
          sin_tra <- cbind(sin_tra, 0, 0, 0, 0)
          colnames(sin_tra)[(ncol(sin_tra)-3):ncol(sin_tra)] <- c("F_SPE", "F_DEP", "F_DIS", "FISH")
          
          #cat(" - Checking Speed... ", sep = "")
          sin_tra[which((sin_tra[,"SPE"]*0.539956803456) > sta_met[sta_met_num, "min_vel_kn"] & (sin_tra[,"SPE"]*0.539956803456) < sta_met[sta_met_num, "max_vel_kn"]), "F_SPE"] <- 1
          #cat("Depth... ", sep = "")
          sin_tra[which(sin_tra[,"DEPTH"] < sta_met[sta_met_num, "min_depth_mt"] & sin_tra[,"DEPTH"] > sta_met[sta_met_num, "max_depth_mt"]), "F_DEP"] <- 1
          #cat("Distance... ", sep = "")
          for(j in 1:nrow(sin_tra))
          {
            nea_har <- which(spDistsN1(harbs, as.numeric(sin_tra[j,c("LON","LAT")]), longlat = TRUE) > (sta_met[sta_met_num, "harb_dist_km"]))
            if(length(nea_har) != 0)
            {
              sin_tra[j, "F_DIS"] <- 1
            }
          }
          
          fis_poi <- which(sin_tra[,"F_SPE"] == 1 & sin_tra[,"F_DEP"] == 1 &  sin_tra[, "F_DIS"] == 1)
          if(length(fis_poi) != 0)
          {
            cat("+", sep = "")
            sin_tra[fis_poi, "FISH"] <- 1
          }
          
          result <- sin_tra[,c("i_id", "F_SPE", "F_DEP", "F_DIS", "FISH")]
          
          rm(sin_tra)
          gc()
          #          cat("\n   -     DB update...    -\n", sep = "")
          if(svalue(met_da_su) == "VMS-LB Match")
          {
            sqldf("insert into p_fish select * from `result`", dbname = vms_DB$db)
          }else{
            sqldf("insert into p_fish_nn select * from `result`", dbname = vms_DB$db)
          }
          rm(result)
          gc()
          
        }
      }
      
    }else{
      
      
      if(svalue(met_da_su) == "VMS-LB Match")
      {
        the_met <- sqldf("select distinct met_des from vms_lb", dbname = vms_DB$db)
      }else{
        the_met <- sqldf("select distinct met_des from nn_clas", dbname = vms_DB$db)
      }
      
      
      to_cle <- which(is.na(the_met[,1]))
      if(length(to_cle) > 0){
        the_met <- the_met[-to_cle,]
      }else{
        the_met <- the_met[,1]
      }
      num_met <- length(the_met)
      
      cat("\n   -     Metier Data Source: ", svalue(met_da_su),"     -\n", sep = "")
      
      for(i in 1:num_met)
      {
        cat("\n   -     Analyzing Metier ", the_met[i], "     -\n", sep = "")
        cat("\n   -     Loading DB data ", sep = "")
        
        if(svalue(met_da_su) == "VMS-LB Match")
        {
          #           cat("from ", svalue(met_da_su), " ", sep = "")
          sin_tra <- fn$sqldf("select * from p_depth join (select intrp.ROWID as i_id, * from intrp join (select vessel as I_NCEE, track as T_NUM, met_des from vms_lb where met_des = '`the_met[i]`')  using (I_NCEE, T_NUM)) using (i_id)", dbname = vms_DB$db)
        }else{
          #           cat("from ", svalue(met_da_su), " ", sep = "")
          sin_tra <- fn$sqldf("select * from p_depth join (select intrp.ROWID as i_id, * from intrp join (select I_NCEE, T_NUM, met_des from nn_clas where met_des = '`the_met[i]`')  using (I_NCEE, T_NUM)) using (i_id)", dbname = vms_DB$db)
        }
        cat(nrow(sin_tra), " points   -\n", sep = "")
        
        if(nrow(sin_tra) == 0)
        {
          
          cat(" - Skipped, no data!")
          next
          
        }
        
        sta_met_num <- which(sta_met[,"Metier"] == the_met[i])
        
        sin_tra <- sin_tra[,-c(2,7,9,10,12,13,14)]
        
        sin_tra <- cbind(sin_tra, 0, 0, 0, 0)
        colnames(sin_tra)[(ncol(sin_tra)-3):ncol(sin_tra)] <- c("F_SPE", "F_DEP", "F_DIS", "FISH")
        
        cat("\n   -     Checking Speed...    -\n", sep = "")
        
        to_spe <- which((sin_tra[,"SPE"]*0.539956803456) > sta_met[sta_met_num, "min_vel_kn"] & (sin_tra[,"SPE"]*0.539956803456) < sta_met[sta_met_num, "max_vel_kn"])
        if(length(to_spe) > 0){sin_tra[to_spe, "F_SPE"] <- 1}
        gc()
        cat("\n   -     Depth...    -\n", sep = "")
        to_dep <- which(sin_tra[,"DEPTH"] < sta_met[sta_met_num, "min_depth_mt"] & sin_tra[,"DEPTH"] > sta_met[sta_met_num, "max_depth_mt"])
        if(length(to_dep) > 0){sin_tra[to_dep, "F_DEP"] <- 1}
        gc()
        cat("\n   -     Distance", sep = "")
        tochk <- which(sin_tra[,"F_SPE"] == 1 & sin_tra[,"F_DEP"] == 1)
        gc()
        if(length(tochk) > 0)
        {
          if(length(tochk)<= 10000){
            cat("...", sep = "")
            dismat <- spDists(as.matrix(sin_tra[tochk,c("LON","LAT")]), as.matrix(harbs@coords), longlat = TRUE)
            tochk <- tochk[-which(dismat < (3*1.85200), arr.ind = TRUE)[,1]]
            rm(dismat)
            gc()
            sin_tra[tochk, "F_DIS"] <- 1
          }else{
            nTock <- ceiling(length(tochk)/10000)
            to_remo <- vector()
            for(nto in 1:nTock)
            {
              cat(".", sep = "")
              r1 <- 10000*(nto-1)+1
              r2 <- min(length(tochk),r1+10000-1)
              dismat <- spDists(as.matrix(sin_tra[tochk[r1:r2],c("LON","LAT")]), as.matrix(harbs@coords), longlat = TRUE)
              to_remo_par <- which(dismat < (3*1.85200), arr.ind = TRUE)[,1]
              to_remo <- c(to_remo, to_remo_par + r1 - 1)
              rm(dismat)
              rm(to_remo_par)
              gc()
            }
            sin_tra[tochk[-to_remo], "F_DIS"] <- 1
            rm(tochk)
            rm(to_remo)
            gc()
          }
        }
        fis_poi <- which(sin_tra[,"F_SPE"] == 1 & sin_tra[,"F_DEP"] == 1 &  sin_tra[, "F_DIS"] == 1)
        if(length(fis_poi) != 0)
        {
          cat("\n   -     ", length(fis_poi)," fishing point founded    -\n", sep = "")
          sin_tra[fis_poi, "FISH"] <- 1
        }
        rm(fis_poi)
        gc()
        
        result <- sin_tra[,c("i_id", "F_SPE", "F_DEP", "F_DIS", "FISH")]
        rm(sin_tra)
        gc()
        cat("\n   -     DB update...    -\n", sep = "")
        if(svalue(met_da_su) == "VMS-LB Match")
        {
          sqldf("insert into p_fish select * from `result`", dbname = vms_DB$db)
        }else{
          sqldf("insert into p_fish_nn select * from `result`", dbname = vms_DB$db)
        }
        rm(result)
        gc()
        
      }
    }
    
    cat("\n\n   ---   END Fishing Point Analysis   ---\n\n", sep = "")
    
    gconfirm("Fishing Point Analysis complete!",
             title = "Confirm",
             icon = "info",
             parent = mark_fis_poi_win,
             handler = function(h,...){dispose(mark_fis_poi_win)})
    
  })
  addSpring(g_input)
  
  
  
  #######################
  
  up_g <- gframe(text = "Filter Parameters", horizontal = TRUE)
  add(big_g, up_g)
  
  #######################
  
  addSpring(up_g)
  lay_win <<- glayout(spacing = 5, container = up_g)
  lay_win[1,1, anchor = 0] <- "Metier"
  lay_win[1,2, anchor = 0] <- "Min Vel Kn"
  lay_win[1,3, anchor = 0] <- "Max Vel Kn"
  lay_win[1,4, anchor = 0] <- "Min Depth Mt"
  lay_win[1,5, anchor = 0] <- "Max Depth Mt"
  lay_win[1,6, anchor = 0] <- "Harb Dist Km"
  addSpring(up_g)
  
  #######################
  
  enabled(b_mark_fis_poi) <- FALSE
  
  if(vms_DB$db != "")
  {
    #     svalue(sel_vms_f) <- strsplit(vms_DB, "/")[[1]][length(strsplit(vms_DB, "/")[[1]])]
    svalue(sel_vms_f) <- ifelse(.Platform$OS.type == "windows", strsplit(vms_DB$db, "\\\\")[[1]][length(strsplit(vms_DB$db, "\\\\")[[1]])],strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])])
    #     nn_tab <- as.numeric(sqldf("SELECT count(*) FROM sqlite_master WHERE type='table' AND name='nn_clas'", dbname = vms_DB$db))
    #     if(nn_tab == 1)
    #     {
    enabled(gri_g3f4) <- TRUE
    #     }
  }
  if(harb$path != "")
  {
    svalue(cus_har_lab) <- paste("Harbour: ", ifelse(.Platform$OS.type == "windows", strsplit(harb$path, "\\\\")[[1]][length(strsplit(harb$path, "\\\\")[[1]])],strsplit(harb$path, "/")[[1]][length(strsplit(harb$path, "/")[[1]])]), sep = "")
  }
  if(harb$path != "" & vms_DB$db != "")
  {
    enabled(b_mark_fis_poi) <- TRUE
    enabled(gri_alg) <- TRUE
  }
  visible(mark_fis_poi_win) <- TRUE
  
}
