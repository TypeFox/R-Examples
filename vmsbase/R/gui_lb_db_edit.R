#' Logbook Editing GUI
#' 
#' The \code{gui_lb_db_edit} function implements the graphical user interface for the
#'  editing of a LogBook Database
#' 
#' In this gui, with a LogBook Database, the user can automatically edit the logbook
#' DB.
#'
#' @param lb_db_name The path of a Logbook DataBase
#'
#' @return This function does not return a value. 
#' After the execution, the Logbook DB will be updated with the edited format.
#' 
#' @usage gui_lb_db_edit(lb_db_name = "")
#' 
#' @export gui_lb_db_edit
#'

gui_lb_db_edit <- function(lb_db_name = "")
{
  lb_DB <- log_DB$new()
  lb_DB$db <- lb_db_name
  
  
  lb_edi_win <- gwindow("Logbook Editing Tool", visible = FALSE, width = 500, height = 300)
  big_g <- ggroup(horizontal = FALSE, container = lb_edi_win)
  
  new_g <- ggroup(horizontal = TRUE, container = big_g)
  addSpring(new_g)
  lb_db_f <- gframe(text = "LogBook DB file", horizontal = TRUE, container = new_g)
  addSpring(lb_db_f)
  sel_lb_f <- glabel("Select LB DB file", container = lb_db_f)
  addSpring(lb_db_f)
  gimage(system.file("ico/document-open-remote.png", package="vmsbase"), container = lb_db_f,
         handler = function(h,...){
           lb_DB$db <<- gfile(text = "Select LB DataBase file",
                              type = "open",
                              filter = list("LB DB file" = list(patterns = c("*.lb.sqlite"))))
#            svalue(sel_lb_f) <- strsplit(lb_DB$db, "/")[[1]][length(strsplit(lb_DB$db, "/")[[1]])]
           svalue(sel_lb_f) <- ifelse(.Platform$OS.type == "windows", strsplit(lb_DB$db, "\\\\")[[1]][length(strsplit(lb_DB$db, "\\\\")[[1]])],strsplit(lb_DB$db, "/")[[1]][length(strsplit(lb_DB$db, "/")[[1]])])
           enabled(start_b) <- TRUE
         })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = lb_db_f,
         handler = function(h,...){
           lb_DB$db <<- ""
           enabled(start_b) <- FALSE
           svalue(sel_lb_f) <- "Select LB DB file"
         })
  addSpring(new_g)
  
  g_sup <- ggroup(horizontal = FALSE, container = big_g)
  addSpring(g_sup)
  
  g_go <- ggroup(horizontal = TRUE, container = g_sup)
  addSpring(g_go)
  start_b <- gbutton(text = "Start\nLogBook Editing", container = g_go, handler = function(h,...)
  {
    distin <- sqldf("select distinct vessUE, s_utc, e_utc from logbook order by s_utc", dbname = lb_DB$db)
    
    cat("\n\n   ---   Logbook Editing Started   ---\n", sep = "")
    
    numrow <- nrow(distin)
    
    species <- sqldf("select distinct specie from logbook order by specie", dbname = lb_DB$db)
    
    spclst <- matrix(data = 0, nrow = 1, ncol = nrow(species))
    
    colnames(spclst) <- paste("FAO_", species[,1], sep = "")
    
    proto_log <- data.frame("vessUE" = numeric(1),
                            "sutc" = numeric(1),
                            "eutc" = numeric(1))
    
    proto_log <- cbind(proto_log, spclst)
    
    hea_str <- paste("vessUE INT,", paste(colnames(proto_log[2:3]), "REAL",  sep = " ", collapse = ", "))
    cre_str <- paste(colnames(proto_log[4:ncol(proto_log)]), "INT", sep = " ", collapse = ", ")
    fin_str <- paste(hea_str, cre_str, sep = ", ")
    
    sqldf("drop table if exists elobo", dbname = lb_DB$db)
    
    query <- paste("CREATE TABLE elobo(", fin_str,")", sep = "")
    
    sqldf(query, dbname = lb_DB$db)
    
    for(n in 1:numrow)
    {
      cat("\nLogbook ", n, " of ", numrow, sep = "")
      new_log <- data.frame("vessUE" = numeric(1), "sutc" = numeric(1), "eutc" = numeric(1))
      
      new_log <- cbind(new_log, spclst)
      
      vess_id <- distin[n,"vessUE"]
      utc_rs  <- distin[n,"s_utc"]
      utc_re  <- distin[n,"e_utc"]
      temp <- fn$sqldf("select specie, qty from logbook where vessUE = `vess_id` and s_utc >= `utc_rs` and e_utc <= `utc_re`", dbname = lb_DB$db)
      
      nclm <- (which(species[,1] %in% temp[,1]))+3
      
      if(max(as.numeric(table(temp[,1])))>1) temp <- as.data.frame(cbind(names(tapply(temp[,2],temp[,1],sum)),as.numeric(tapply(temp[,2],temp[,1],sum))))
      
      new_log["vessUE"] <- vess_id
      new_log["eutc"] <- utc_re
      new_log["sutc"] <- utc_rs
      
      new_log[nclm] = new_log[nclm]+abs(as.numeric(temp[,2]))
      
      sqldf("insert into elobo select * from `new_log`", dbname = lb_DB$db)
    }
    
    cat("\n\n   ---   End Logbook Editing   ---\n\n", sep = "")
    
    
    gconfirm(paste("LogBook DB editing complete!\n\n","Edited ", n, " records", sep = ""),
             title = "Confirm",
             icon = "info")
  })
  enabled(start_b) <- FALSE
  addSpring(g_go)
  
  
  if(lb_DB$db != "")
  {
#     svalue(sel_lb_f) <- strsplit(lb_DB$db, "/")[[1]][length(strsplit(lb_DB$db, "/")[[1]])]
    svalue(sel_lb_f) <- ifelse(.Platform$OS.type == "windows", strsplit(lb_DB$db, "\\\\")[[1]][length(strsplit(lb_DB$db, "\\\\")[[1]])],strsplit(lb_DB$db, "/")[[1]][length(strsplit(lb_DB$db, "/")[[1]])])
    enabled(start_b) <- TRUE
  }
  
  visible(lb_edi_win) <- TRUE
  
}