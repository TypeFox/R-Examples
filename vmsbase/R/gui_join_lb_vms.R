
#' VMS-LB Match GUI
#'  
#' 
#' The \code{gui_join_lb_vms} function implements the graphical user interface for the
#'  VMS-LB Matching routine.
#' 
#' This function, with both VMS and LogBook databases,
#'  performs a VMS-LB track match based on vessel, track and time conformity.
#'   
#' @param vms_db_name The path of a VMS DataBase
#' @param lb_db_name The path of a LogBook DataBase
#' 
#' @return This function does not return a value. 
#' 
#' 
#' @usage gui_join_lb_vms(lb_db_name = "", vms_db_name = "")
#' 
#' @export gui_join_lb_vms
#'
#'
#'@references free text reference Pointers to the literature related to this object.



gui_join_lb_vms <- function(lb_db_name = "", vms_db_name = "")
{
  
  lb_DB <- log_DB$new()
  lb_DB$db <- lb_db_name
  vms_DB <- vms_DB$new()
  vms_DB$db <- vms_db_name
  
  join_lb_vms_win <- gwindow(title = "LogBook - VMS Matching Tool",
                             visible = FALSE,
                             width = 600,
                             height = 350)
  
  big_g <- ggroup(horizontal = FALSE, container = join_lb_vms_win)
  
  g_input <- ggroup(horizontal = TRUE, container = big_g)
  addSpring(big_g)
  g_go <- ggroup(horizontal = FALSE, container = big_g)
  
  addSpring(g_input)
  lb_db_f <- gframe(text = "LogBook DB file", horizontal = TRUE, container = g_input)
  addSpring(lb_db_f)
  sel_lb_f <- glabel("Select LB DB file", container = lb_db_f)
  addSpring(lb_db_f)
  gimage(system.file("ico/folder-orange.png", package="vmsbase"), container = lb_db_f,
         handler = function(h,...){
           lb_DB$db <<- gfile(text = "Select LB DataBase file",
                              type = "open",
                              filter = list("LB DB file" = list(patterns = c("*.lb.sqlite"))))
           svalue(sel_lb_f) <- ifelse(.Platform$OS.type == "windows", strsplit(lb_DB$db, "\\\\")[[1]][length(strsplit(lb_DB$db, "\\\\")[[1]])],strsplit(lb_DB$db, "/")[[1]][length(strsplit(lb_DB$db, "/")[[1]])])
           
           if(vms_DB$db != "")
           {
             enabled(b_lb_vms_match) <- TRUE
           }
         })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = lb_db_f,
         handler = function(h,...){
           lb_DB$db <<- ""
           enabled(b_lb_vms_match) <- FALSE
           svalue(sel_lb_f) <- "Select LB DB file"
         })
  addSpring(g_input)
  gimage(system.file("ico/insert-link-64.png", package="vmsbase"), container = g_input)
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
           svalue(sel_vms_f) <- ifelse(.Platform$OS.type == "windows", strsplit(vms_DB$db, "\\\\")[[1]][length(strsplit(vms_DB$db, "\\\\")[[1]])],strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])])
           if(lb_DB$db != "")
           {
             enabled(b_lb_vms_match) <- TRUE
           }
         })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = vms_db_f,
         handler = function(h,...){
           vms_DB$db <- ""
           enabled(b_lb_vms_match) <- FALSE
           svalue(sel_vms_f) <- "Select VMS DB file"
         })
  addSpring(g_input)
  addSpring(g_go)
  
  info_lab <- glabel("", container = g_go)
  #   addSpring(g_go)
  b_lb_vms_match <- gbutton(text = "   Start\nMatching", container = g_go, handler = function(h,...)
  {
    enabled(g_input) <- FALSE
    enabled(b_lb_vms_match) <- FALSE
    
    svalue(info_lab) <- paste("VMS-Logbook Matching Started", sep = "")
    cat("\n   ---   VMS-Logbook Matching Started   ---\n", sep = "")
    
    incee_vms <- sqldf("select distinct I_NCEE from track", dbname = vms_DB$db)
    incee_lb <- sqldf("select distinct vessUE from elobo", dbname = lb_DB$db)
    
    sqldf("drop table if exists vms_lb", dbname = vms_DB$db)
    
    query <- "CREATE TABLE vms_lb(vessel INT, track INT, logbook INT, log_id INT, met_des CHAR)"
    
    sqldf(query, dbname = vms_DB$db)
    
    vess <- incee_lb[incee_lb[,1] %in% incee_vms[,1], 1]
    
    num_vess <- length(vess)
    
    if(num_vess > 0 )
    {
      for(i in 1:num_vess)
      {
        svalue(info_lab) <- paste("Vessel: ", vess[i], " - N.", i, " of ", num_vess, sep = "")
        cat("\n   -   Vessel: ", vess[i], " - N.", i, " of ", num_vess, sep = "")
        vms_data <- fn$sqldf("select * from track where I_NCEE = `vess[i]` order by DATE", dbname = vms_DB$db)
        lb_data <- fn$sqldf("select elobo.ROWID, * from elobo, lb_cla where vessUE = `vess[i]` and log_num = elobo.rowid order by s_utc, e_utc", dbname = lb_DB$db)
        
        #                                 if(is.null(nrow(vms_data["T_NUM"])) | nrow(lb_data) == 1)
        if(nrow(vms_data) == 0 | nrow(lb_data) == 0)
        {
          cat(" - Skipped, not enough data!\n")
          next
        }else{
        
          num_track <- max(vms_data["T_NUM"])
          
          res_over <- data.frame(vessel = vess[i],
                                 track = unique(vms_data[,"T_NUM"]),
                                 logbook = numeric(num_track),
                                 log_id = numeric(num_track),
                                 met_fou = character(num_track))
          
          res_over[,"met_fou"] <- NA
          
          cat(" - ", num_track, " tracks ", sep = "")
          for( k in 1:num_track)
          {
            trakap <- which(vms_data["T_NUM"] == k)
            if(length(trakap) == 0)
            {
              cat("-", sep = "")
              next
            }
            cat(".", sep = "")
            min_tr <- min(vms_data[trakap,"DATE"])
            max_tr <- max(vms_data[trakap,"DATE"])
            
            if(!is.na(min_tr) | !is.na(max_tr))
            {
              if(min_tr != max_tr)
              {
                int_tr <- Intervals(c(min_tr, max_tr))
                
                int_lb <- Intervals(cbind(lb_data[,"s_utc"], lb_data[,"e_utc"]))
                
                interval <- interval_intersection(int_tr, int_lb)
                
                overlap <- (interval_overlap(int_tr, int_lb))
                
                best <- which.max(apply(as.matrix(interval),1,diff))
                
                log_num <- unlist(overlap[[1]][best])
                
                if(length(log_num) != 0)
                {
                  cat("+", sep = "")
                  res_over[k, "logbook"] <- log_num
                  res_over[k, "log_id"] <- lb_data[log_num, "rowid"]
                  res_over[k, "met_fou"] <- lb_data[log_num, "met_des"]
                }
              }
            }
          }
          no_lobo <- which(res_over[,"logbook"] == 0 & res_over[,"log_id"] == 0)
          if(length(no_lobo) != 0)
          {
            res_over <- res_over[-no_lobo,]
            if(nrow(res_over) == 0){
              cat(" *", sep = "")
              next
            }
          }
          
          sqldf("INSERT INTO vms_lb SELECT * FROM `res_over`", dbname = vms_DB$db)
        }
      }
    }else{
      cat("\n\n   ---   STOP   -   No VMS-LogBook Vessel Name Match Found!   ---\n\n", sep = "")      
    }
    
    no_ves <- length(sqldf("select distinct vessel from vms_lb", dbname = vms_DB$db)[,1])
    no_mes <- length(sqldf("select distinct met_des from vms_lb", dbname = vms_DB$db)[,1])
    no_elog <- length(sqldf("select distinct log_id from vms_lb", dbname = vms_DB$db)[,1])
    no_mat <- nrow(sqldf("select distinct vessel, track from vms_lb", dbname = vms_DB$db))
    no_cov <- nrow(sqldf("select distinct I_NCEE, T_NUM from track", dbname = vms_DB$db))
    no_lobo <- sqldf("select count(*) from elobo", dbname = lb_DB$db)[1,]

    
    cat("\n\n   -     VMS DB tracks: ",  no_cov, "  -  LB DB logs: ", no_lobo, "  -  estimated VMS Coverage: ", round(100/no_cov*no_lobo, 2), "%     -", sep = "")
    cat("\n   -     VMS-LB match: ",  no_mat, "  -  LB*VMS track: ", round(no_elog/no_mat, 2) , "  -  real VMS Coverage: ", round(100/no_cov*no_mat, 2), "%     -", sep = "")
    cat("\n   -     Matched Metier: ", no_mes, "     -\n", sep = "")
    
    cat("\n\n   ---   END VMS-Logbook Matching   ---\n\n", sep = "")
    svalue(info_lab) <- paste("  ---  END VMS-Logbook Matching  ---",
                              "\n\n   -     VMS DB tracks: ",  no_cov, " - LB DB logs: ", no_lobo, " - estimated VMS Coverage: ", round(100/no_cov*no_lobo, 2), "%   -",
                              "\n   -     VMS-LB match: ",  no_mat, " - LB*VMS track: ", round(no_elog/no_mat, 2) , " - real VMS Coverage: ", round(100/no_cov*no_mat, 2), "%   -",
                              "\n   -     Matched Metier: ", no_mes, "     -\n\n\n", sep = "")
    gconfirm("LogBook - Vms Matching Complete!",
             title = "Confirm",
             icon = "info",
             parent = join_lb_vms_win,
             handler = dispose(join_lb_vms_win))
    
  })
  addSpring(g_go)
  enabled(b_lb_vms_match) <- FALSE
  
  if(lb_DB$db != "")
  {
    svalue(sel_lb_f) <- ifelse(.Platform$OS.type == "windows", strsplit(lb_DB$db, "\\\\")[[1]][length(strsplit(lb_DB$db, "\\\\")[[1]])],strsplit(lb_DB$db, "/")[[1]][length(strsplit(lb_DB$db, "/")[[1]])])
  }
  if(vms_DB$db != "")
  {
    svalue(sel_vms_f) <- ifelse(.Platform$OS.type == "windows", strsplit(vms_DB$db, "\\\\")[[1]][length(strsplit(vms_DB$db, "\\\\")[[1]])],strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])])
  }
  if(lb_DB$db != "" & vms_DB$db != "")
  {
    enabled(b_lb_vms_match) <- TRUE 
  }
  
  visible(join_lb_vms_win) <- TRUE
}