
#' VMS DataBase Creation GUI
#'
#' The \code{gui_vmsdb_mixsou} function implements the graphical user interface for the
#' creation of a VMS DataBase from multiple sources.
#'
#' This function, with an edited vms dataset (see \code{\link{gui_vms_editraw}}),
#' creates a new VMS DataBase.
#'
#' @return This function does not return a value.
#' After the execution a VMS DataBase will be deployed.
#'
#' @usage gui_vmsdb_mixsou()
#'
#' @export gui_vmsdb_mixsou
#'
#'@seealso \code{\link{gui_vms_editraw}}
#'

### VMS DataBase manager gui
gui_vmsdb_mixsou <- function()
{
  cur_sou <- miso_list$new()
  
  dim_plo <- matrix(c(1,1,1,2), 1,4,byrow = TRUE)
  
  vms_db_win <- gwindow("Merging sources", visible = TRUE, width = 900, height= 500, toolkit = guiToolkit())
  g_big <- ggroup(horizontal = FALSE, container = vms_db_win)
  
  g_do <- ggroup(horizontal = TRUE, container = g_big)
  g_do_lef <- ggroup(horizontal = FALSE, container = g_do)
  
  g_lef_up <- ggroup(horizontal = TRUE, container = g_do_lef)
  addSpring(g_lef_up)
  gimage(system.file("ico/edit-add-2.ico", package="vmsbase"), container = g_lef_up, handler = function(h,...){
    svalue(stat_bar) <- "Select file to load"
    new_db <- gsub('\\\\', '/', gfile(text = "Select Fleet Data",type = "open",
                                      filter = list("VMS data" = list(patterns = c("*.vms")),
                                                    "All files" = list(patterns = c("*")))))
    svalue(stat_bar) <- paste("Loading dataset: ", new_db, sep = "")
    if(new_db != ""){
      enabled(vms_db_win) <- FALSE
      new_dataset <- mix_sou$new(db_path = new_db, main_cont = g_do_lef, sel_wid = sel_view, cur_sou)
      
      cur_sou$add_sou(new_dataset)
      
      svalue(stat_bar) <- paste("Plotting ", svalue(sel_view),"...", sep = "")
      enabled(vms_db_win) <- FALSE
      plot_comp(svalue(sel_view), cur_sou$th_lst)
      enabled(vms_db_win) <- TRUE
      enabled(g_rig_f) <- TRUE
      enabled(g_bot_f) <- TRUE
      svalue(stat_bar) <- ""
      
    }
  })
  glabel("Add Dataset", container = g_lef_up)
  addSpring(g_lef_up)
  g_do_rig <- ggroup(horizontal = FALSE, container = g_do)
  
  ggraphics(width = 800, height = 500, dpi = 100, container = g_do_rig)
  
  g_do_bot <- ggroup(horizontal = TRUE, container = g_do_rig)
  
  addSpring(g_do_bot)
  g_rig_f <- gframe("View", container = g_do_bot, horizontal = TRUE)
  sel_view <- gradio(c("Absolute", "Relative", "IDs", "Stats", "Box"), selected = 1, horizontal = TRUE, container = g_rig_f, handler = function(h,...){
    svalue(stat_bar) <- paste("Plotting ", svalue(sel_view),"...", sep = "")
    enabled(vms_db_win) <- FALSE
    plot_comp(svalue(sel_view), cur_sou$th_lst)
    enabled(vms_db_win) <- TRUE
    svalue(stat_bar) <- ""
  })
  addSpring(g_do_bot)
  g_bot_f <- gframe("Database", container = g_do_bot, horizontal = TRUE)
  gimage(system.file("ico/db_add.ico", package="vmsbase"), container = g_bot_f, handler = function(h,...){
    
    enabled(vms_db_win) <- FALSE
    
    new_db_name <- ginput("\nEnter the name of the new database", text = "", title = "New Database", icon = "info", parent = vms_db_win)
    
    db_dir <- gfile(text = "Select VMS DB destination",
                    type = "selectdir",
                    filter = list("VMS DB data" = list(patterns = c("*.vms.sqlite"))))
    
    
    mix_db <- paste(db_dir, "/", new_db_name, ".vms.sqlite", sep = "")
    dbConnect(SQLite(), dbname = mix_db)

    obj_lst <- cur_sou$th_lst
    for(i in 1:length(obj_lst)){
      if(i == 1){
        fn$read.csv.sql(obj_lst[[i]]$dbPath, sql = "CREATE TABLE ping AS SELECT *, '`obj_lst[[i]]$sourceName`' as SOU FROM file", dbname = mix_db, sep = ";", eol = "\n")
      }else{
        tmp_matc <- obj_lst[[i]]$ids_tab
        if(length(tmp_matc) == 0){
          fn$read.csv.sql(obj_lst[[i]]$dbPath, sql = "insert into ping SELECT *, '`obj_lst[[i]]$sourceName`' as SOU FROM file", dbname = mix_db, sep = ";", eol = "\n")
        }else{
          fn$read.csv.sql(obj_lst[[i]]$dbPath, sql = "insert into ping SELECT *, '`obj_lst[[i]]$sourceName`' as SOU FROM (select V2 as I_NCEE, LAT, LON, DATE, SPE, HEA from (SELECT * FROM file join (select V1 as I_NCEE, * from tmp_matc) using (I_NCEE)))", dbname = mix_db, sep = ";", eol = "\n")
          fn$read.csv.sql(obj_lst[[i]]$dbPath, sql = "insert into ping select *, '`obj_lst[[i]]$sourceName`' as SOU from file where I_NCEE not in (select V1 from tmp_matc)", dbname = mix_db, sep = ";", eol = "\n")
        }
      }
    }
    
    sqldf("CREATE TABLE warn(p_id INT, W_DUPL INT, W_HARB INT, W_LAND INT, W_COHE INT)", dbname = mix_db)
    sqldf("CREATE TABLE track(I_NCEE INT, LAT REAL, LON REAL, DATE REAL, SPE REAL, HEA REAL, W_HARB INT, T_NUM INT, P_ID INT)", dbname = mix_db)
    sqldf("CREATE TABLE intrp(I_NCEE INT, LAT REAL, LON REAL, DATE REAL, SPE REAL, HEA REAL, W_HARB INT, T_NUM INT, P_ID INT, P_INT INT, T_ID INT)", dbname = mix_db)
    sqldf("CREATE TABLE p_depth(i_id INT, vess_id INT, DEPTH REAL)", dbname = mix_db)
    sqldf("CREATE TABLE p_area(vess_id INT, t_num INT, AREA INT)", dbname = mix_db)
    sqldf("CREATE TABLE vms_lb(vessel INT, track INT, logbook INT, log_id INT, met_des CHAR)", dbname = mix_db)
    sqldf("CREATE TABLE p_fish(i_id INT, F_SPE INT, F_DEP INT, F_DIS INT, FISH INT)", dbname = mix_db)
    
    gconfirm("VMS DataBase Deploy Completed!",
             title = "Confirm",
             icon = "info",
             parent = vms_db_win,
             handler = function(h,...){dispose(vms_db_win)})
    
  })
  glabel("Create New", container = g_bot_f)
  stat_bar <- gstatusbar(text = "Click on the green plus to add a dataset...", container = g_big)
  enabled(g_rig_f) <- FALSE
  enabled(g_bot_f) <- FALSE
  #   visible(vms_db_win) <- TRUE
}
