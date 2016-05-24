
#' LogBook DB Status GUI
#'  
#' 
#' The \code{gui_lb_db_stat} function implements the graphical user interface for the
#'  LogBook DB Status viewer.
#' 
#' This function, with a LogBook database,
#'  shows the current LogBook DB status.
#'   
#' @param lb_db_name The path of a LogBook DataBase
#' 
#' @return This function does not return a value. 
#' 
#' @usage gui_lb_db_stat(lb_db_name = "")
#' 
#' @export gui_lb_db_stat
#'

gui_lb_db_stat <- function(lb_db_name = "")
{
  lb_DB <- log_DB$new()
  lb_DB$db <- lb_db_name
  
  lb_stat_win <- gwindow("LogBook DataBase Status", visible = FALSE, height = 600)
  
  big_g <- ggroup(horizontal = FALSE, container = lb_stat_win, spacing = 0)
  chk_g3 <- ggroup(horizontal = TRUE, container = big_g)
  top_g <- gframe(horizontal = TRUE, container = big_g, expand = F)

  addSpring(chk_g3)
  lb_db_f <- gframe(text = "LogBook DB file", horizontal = TRUE, container = chk_g3)
  addSpring(lb_db_f)
  sel_lb_f <- glabel("Select LogBook DB file", container = lb_db_f)
  addSpring(lb_db_f)
  gimage(system.file("ico/folder-orange.png", package="vmsbase"), container = lb_db_f,
         handler = function(h,...){
           lb_DB$db <- gfile(text = "Select LogBook DataBase file",
                             type = "open",
                             filter = list("LogBook DB file" = list(patterns = c("*.lb.sqlite"))))
#            svalue(sel_lb_f) <- strsplit(lb_DB$db, "/")[[1]][length(strsplit(lb_DB$db, "/")[[1]])]
           svalue(sel_lb_f) <- ifelse(.Platform$OS.type == "windows", strsplit(lb_DB$db, "\\\\")[[1]][length(strsplit(lb_DB$db, "\\\\")[[1]])],strsplit(lb_DB$db, "/")[[1]][length(strsplit(lb_DB$db, "/")[[1]])])
           n_log <- sqldf("select count(*) from elobo", dbname = lb_DB$db)
           if (n_log > 0)
           {
             svalue(log_n_l) <- paste(" N. of Logs: ", n_log, sep = "")
             svalue(log_spe_l) <- paste(" N. of Species: ", ncol(sqldf("select * from elobo limit 1", dbname = lb_DB$db))-3, sep = "")
             #svalue(log_not_l) <- paste("Invalid Logs: " , nrow(sqldf("select distinct vessUE, s_utc, e_utc from logbook order by s_utc, e_utc", dbname = lb_DB$db)) - n_log, sep = "")
             log_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
             enabled(log_b) <- TRUE
           }else{
             svalue(log_n_l) <- ""
             #svalue(log_not_l) <- ""
             svalue(log_spe_l) <- ""
             log_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
             enabled(log_b) <- FALSE
           }
           delete(log_g, log_sta)
           log_sta <<- log_sta_n
           add(log_g, log_sta)
         })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = lb_db_f,
         handler = function(h,...){
           lb_DB$db <- ""
           log_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
           svalue(log_n_l) <- ""
           #svalue(log_not_l) <- ""
           svalue(log_spe_l) <- ""
           enabled(log_b) <- FALSE
           delete(log_g, log_sta)
           log_sta <<- log_sta_n
           add(log_g, log_sta)
           svalue(sel_lb_f) <- "Select LogBook DB file"
         })
  ################
  dat_ref <- gimage(system.file("ico/document-quick_restart.png", package="vmsbase"))
  add(lb_db_f, dat_ref)
  addHandlerClicked(dat_ref, handler = function(h,...)
  {
    if(lb_DB$db != "")
    {
      n_log <- sqldf("select count(*) from elobo", dbname = lb_DB$db)
      if (n_log > 0)
      {
        svalue(log_n_l) <- paste(" N. of Logs: ", n_log, sep = "")
        svalue(log_spe_l) <- paste(" N. of Species: ", sqldf("select count(*) from elobo", dbname = lb_DB$db), sep = "")
        #svalue(log_not_l) <- paste("Invalid Logs: " , nrow(sqldf("select distinct vessUE, s_utc, e_utc from logbook order by s_utc, e_utc", dbname = lb_DB$db)) - n_log, sep = "")
        log_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
        enabled(log_b) <- TRUE
      }else{
        svalue(log_n_l) <- ""
        svalue(log_spe_l) <- ""
        #svalue(log_not_l) <- ""
        log_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
        enabled(log_b) <- FALSE
      }
      delete(log_g, log_sta)
      log_sta <<- log_sta_n
      add(log_g, log_sta)
    }
  })
  addSpring(lb_db_f)
  addSpring(chk_g3)
  addSpring(top_g)  
  log_g <- gframe(text = "LogBook", horizontal = FALSE, container = top_g)
  log_n_l <- glabel("", container = log_g)
  #log_not_l <- glabel("", container = log_g)
  log_spe_l <- glabel("", container = log_g)
  addSpring(log_g)
  log_b <- gbutton(text = "Show data", container = log_g, handler = function(h,..){
    
    tb_dat <- sqldf("select * from elobo", dbname = lb_DB$db)
    colnames(tb_dat) <- sub("FAO_", "", colnames(tb_dat))
    logs <- gtable(tb_dat, expand = T)
    add(data_nb, logs, label = "LogBook")
  })
  log_sta <- gimage(system.file("ico/user-invisible.png", package="vmsbase"))
  add(log_g, log_sta)
  addSpring(top_g)
  gbutton("FAO Codes", container = top_g, handler = function(h,..)
    {
    browseURL(url = "www.afma.gov.au/wp-content/uploads/2010/06/fao_species_codes.xls", browser = getOption("browser"))
  })
  addSpring(top_g)
  
  bot_g <- ggroup(container = big_g, expand = T)
  data_nb <- gnotebook(tab.pos = 3, closebuttons = TRUE, container = bot_g, expand = T)
  enabled(log_b) <- FALSE
  visible(lb_stat_win) <- TRUE
  if(lb_DB$db != "") 
  {
#     svalue(sel_lb_f) <- strsplit(lb_DB$db, "/")[[1]][length(strsplit(lb_DB$db, "/")[[1]])]
    svalue(sel_lb_f) <- ifelse(.Platform$OS.type == "windows", strsplit(lb_DB$db, "\\\\")[[1]][length(strsplit(lb_DB$db, "\\\\")[[1]])],strsplit(lb_DB$db, "/")[[1]][length(strsplit(lb_DB$db, "/")[[1]])])
    n_log <- sqldf("select count(*) from elobo", dbname = lb_DB$db)
    if (n_log > 0)
    {
      svalue(log_n_l) <- paste("N. of Logs: ", n_log, sep = "")
      #svalue(log_not_l) <- paste("Invalid Logs: ", nrow(sqldf("select distinct vessUE, s_utc, e_utc from logbook order by s_utc, e_utc", dbname = lb_DB$db)) - n_log, sep = "")
      svalue(log_spe_l) <- paste("N. of Species: ", sqldf("select count(*) from elobo", dbname = lb_DB$db), sep = "")
      log_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
      enabled(log_b) <- TRUE
    }else{
      svalue(log_n_l) <- ""
      #svalue(log_not_l) <- ""
      svalue(log_spe_l) <- ""
      log_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
      enabled(log_b) <- FALSE
    }
    delete(log_g, log_sta)
    log_sta <<- log_sta_n
    add(log_g, log_sta)
  }
}
