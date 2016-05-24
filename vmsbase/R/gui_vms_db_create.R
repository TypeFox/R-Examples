

#' VMS DataBase Creation GUI
#' 
#' The \code{gui_vms_db_create} function implements the graphical user interface for the
#' creation of a VMS DataBase.
#' 
#' This function, with an edited vms dataset (see \code{\link{gui_vms_editraw}}),
#'  creates a new VMS DataBase.
#'
#' @return This function does not return a value. 
#' After the execution a VMS DataBase will be deployed.
#' 
#' @usage gui_vms_db_create()
#' 
#' @export gui_vms_db_create
#'
#'@seealso \code{\link{gui_vms_editraw}}

### VMS DataBase manager gui

gui_vms_db_create <- function()
{
  vms_DB <- vms_DB$new()
  vms_db_win <- gwindow("Create New VMS DataBase", visible = FALSE)
  one <- ggroup(horizontal = FALSE, container = vms_db_win)
  glabel("To proceed with the creation of a VMS Database,
  first enter a name for the new DB and press enter...", container = one)
  c_name <- ggroup(horizontal = TRUE, container = one)
  addSpring(c_name)
  dbname <- gedit(initial.msg = "New VMS DB name...",
                  container = c_name, handler = function(h,...)
                  {
                    enabled(data) <- TRUE
                  }) 
  addSpring(c_name)
  gimage(system.file("ico/go-down-3.png", package="vmsbase"), container = one)
  glabel("Then select an edited VMS file 
  clicking on the button below...", container = one)
  data <- gbutton(text = "Select VMS data", container = one, handler = function(h,...)
  {
    svalue(data) <- gfile(text = "Select VMS data",
                          type = "open",
                          filter = list("VMS data" = list(patterns = c("*.vms")),
                                        "All files" = list(patterns = c("*"))))
    enabled(go) <- TRUE
  })
  
  gimage(system.file("ico/go-down-3.png", package="vmsbase"), container = one)
  glabel("At last provide a destination for the VMS DB 
      clicking on the \'Create DB\' button.", container = one)
  go <- gbutton("Create DB", container = one, handler = function(h,...)
  {
    enabled(vms_db_win) <- FALSE
    
    vms_DB$dir <- gfile(text = "Select VMS DB destination",
                       type = "selectdir",
                       filter = list("VMS DB data" = list(patterns = c("*.vms.sqlite"))))
    
    vms_DB$db <- paste(vms_DB$dir, "/", svalue(dbname), ".vms.sqlite", sep = "")  
    
    dbConnect(SQLite(), dbname = vms_DB$db)
    
    read.csv.sql(svalue(data), sql = "CREATE TABLE ping AS SELECT * FROM file", dbname = vms_DB$db, sep = ";", eol = "\n")
    sqldf("CREATE TABLE warn(p_id INT, W_DUPL INT, W_HARB INT, W_LAND INT, W_COHE INT)", dbname = vms_DB$db)
    sqldf("CREATE TABLE track(I_NCEE INT, LAT REAL, LON REAL, DATE REAL, SPE REAL, HEA REAL, W_HARB INT, T_NUM INT, P_ID INT)", dbname = vms_DB$db)
    sqldf("CREATE TABLE intrp(I_NCEE INT, LAT REAL, LON REAL, DATE REAL, SPE REAL, HEA REAL, W_HARB INT, T_NUM INT, P_ID INT, P_INT INT, T_ID INT)", dbname = vms_DB$db)
    sqldf("CREATE TABLE p_depth(i_id INT, vess_id INT, DEPTH REAL)", dbname = vms_DB$db)
    sqldf("CREATE TABLE p_area(vess_id INT, t_num INT, AREA INT)", dbname = vms_DB$db)
    sqldf("CREATE TABLE vms_lb(vessel INT, track INT, logbook INT, log_id INT, met_des CHAR)", dbname = vms_DB$db)
    sqldf("CREATE TABLE p_fish(i_id INT, F_SPE INT, F_DEP INT, F_DIS INT, FISH INT)", dbname = vms_DB$db)
    
          gconfirm("VMS DataBase Deploy Completed!",
             title = "Confirm",
             icon = "info",
             parent = vms_db_win,
             handler = function(h,...){dispose(vms_db_win)})
  })
  
  enabled(data) <- FALSE
  enabled(go) <- FALSE
  visible(vms_db_win) <- TRUE
}