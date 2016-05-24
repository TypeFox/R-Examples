
#' VMS DataBase Status GUI
#'  
#' The \code{gui_vms_db_stat} function implements the graphical user interface for the
#'  VMS DataBase Status viewer.
#' 
#' This function, with a VMS DataBase,
#'  shows the current VMS DataBase status.
#'   
#' @param vms_db_name The path of a VMS DataBase
#' 
#' @return This function does not return a value. 
#' 
#' @usage gui_vms_db_stat(vms_db_name = "")
#' 
#' @export gui_vms_db_stat
#'

gui_vms_db_stat <- function(vms_db_name = "")
{
  vms_DB <- vms_DB$new()
  vms_DB$db <- vms_db_name
  
  max_date <- numeric()
  min_date <- numeric()
  
  vms_select_win <- gwindow("VMS DataBase Status", visible = FALSE, height = 600)
  
  ##################
  
  big_g <- ggroup(horizontal = FALSE, container = vms_select_win, spacing = 0)
  
  sup_g <- ggroup(horizontal = TRUE, container = big_g)
  
  top_g <- gframe(horizontal = TRUE, container = big_g, expand = F)
  
  chk_g3 <- ggroup(horizontal = TRUE, container = sup_g)
  
  #################
  addSpring(chk_g3)
  vms_db_f <- gframe(text = "VMS DB file", horizontal = TRUE, container = chk_g3)
  addSpring(vms_db_f)
  sel_vms_f <- glabel("Select VMS DB file", container = vms_db_f)
  addSpring(vms_db_f)
  gimage(system.file("ico/folder-blue.png", package="vmsbase"), container = vms_db_f,
         handler = function(h,...){
           vms_DB$db <- gfile(text = "Select VMS DataBase file",
                              type = "open",
                              filter = list("VMS DB file" = list(patterns = c("*.vms.sqlite"))))
           #            svalue(sel_vms_f) <- strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])]
           svalue(sel_vms_f) <- ifelse(.Platform$OS.type == "windows", strsplit(vms_DB$db, "\\\\")[[1]][length(strsplit(vms_DB$db, "\\\\")[[1]])],strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])])
           #################
           n_png <- sqldf("select count(*) from ping", dbname = vms_DB$db)
           if (n_png > 0)
           {
             png_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
             svalue(num_vess) <- paste(" N. of Vessels: \n", sqldf("select count(distinct I_NCEE) from ping", dbname = vms_DB$db), sep = "")
             svalue(num_png) <- paste(" N. of Pings: \n", sqldf("select count(*) from ping", dbname = vms_DB$db), sep = "")
             
             min_date <<- as.numeric(sqldf("select min(DATE) from ping", dbname = vms_DB$db)[,1])
             max_date <<- as.numeric(sqldf("select max(DATE) from ping", dbname = vms_DB$db)[,1])
             
             svalue(min_date_l) <- paste("First: ", days(min_date), "/", months(min_date), "/", years(min_date), sep = "")
             svalue(max_date_l) <- paste("Last: ", days(max_date), "/", months(max_date), "/", years(max_date), sep = "")
             
             enabled(png_b) <- TRUE
           }else{
             png_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
             svalue(num_vess) <- "   ---"
             svalue(num_png) <- "   ---"
             enabled(png_b) <- FALSE
           }
           delete(png_g, png_sta)
           png_sta <<- png_sta_n
           add(png_g, png_sta)
           
           n_wrn <- sqldf("select count(*) from warn", dbname = vms_DB$db)
           if (n_wrn > 0)
           {
             wrn_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
             svalue(num_wrn) <- paste(" N. of Warnings: \n", sqldf("select count(*) from warn", dbname = vms_DB$db), sep = "")
             enabled(wrn_b) <- TRUE
           }else{
             wrn_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
             svalue(num_wrn) <- "   ---"
             enabled(wrn_b) <- FALSE
           }
           delete(wrn_g, wrn_sta)
           wrn_sta <<- wrn_sta_n
           add(wrn_g, wrn_sta)
           
           n_trk <- sqldf("select count(*) from track", dbname = vms_DB$db)
           if (n_trk > 0)
           {
             trk_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
             svalue(mean_trk) <- paste(" Mean tracks for vessel: \n", round(sqldf("select count(*) from (select distinct I_NCEE, T_NUM from track)", dbname = vms_DB$db)/sqldf("select count(distinct I_NCEE) from track", dbname = vms_DB$db), 1), sep = "")
             svalue(tot_trk) <- paste(" N. of tracks: \n", sqldf("select count(*) from (select distinct I_NCEE, T_NUM from track)", dbname = vms_DB$db), sep = "")
             enabled(trk_b) <- TRUE
             enabled(sel_vms) <- TRUE  
           }else{
             trk_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
             svalue(mean_trk) <- "   ---"
             svalue(tot_trk) <- "   ---"
             enabled(trk_b) <- FALSE
           }
           delete(trk_g, trk_sta)
           trk_sta <<- trk_sta_n
           add(trk_g, trk_sta)
           
           n_ntr <- sqldf("select count(*) from intrp", dbname = vms_DB$db)
           if (n_ntr > 0)
           {
             ntr_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
             svalue(num_ntr) <- paste(" N. of Pings: \n", sqldf("select count(*) from intrp", dbname = vms_DB$db), sep = "")
             enabled(ntr_b) <- TRUE
           }else{
             ntr_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
             svalue(num_ntr) <- "   ---"
             enabled(ntr_b) <- FALSE
           }
           delete(ntr_g, ntr_sta)
           ntr_sta <<- ntr_sta_n
           add(ntr_g, ntr_sta)
           
           n_dep <- sqldf("select count(*) from p_depth", dbname = vms_DB$db)
           if (n_dep > 0)
           {
             dep_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
             svalue(num_dep) <- paste(" N. of Pings: \n", sqldf("select count(*) from p_depth", dbname = vms_DB$db), sep = "")
             enabled(dep_b) <- TRUE
           }else{
             dep_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
             svalue(num_dep) <- "   ---"
             enabled(dep_b) <- FALSE
           }
           delete(dep_g, dep_sta)
           dep_sta <<- dep_sta_n
           add(dep_g, dep_sta)
           
           n_are <- sqldf("select count(*) from p_area", dbname = vms_DB$db)
           if (n_are > 0)
           {
             are_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
             svalue(num_are) <- paste(" N. of Pings: \n", sqldf("select count(*) from p_area", dbname = vms_DB$db), sep = "")
             enabled(are_b) <- TRUE
           }else{
             are_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
             svalue(num_are) <- "   ---"
             enabled(are_b) <- FALSE
           }
           delete(are_g, are_sta)
           are_sta <<- are_sta_n
           add(are_g, are_sta)
           
           n_mat <- sqldf("select count(*) from vms_lb", dbname = vms_DB$db)
           if (n_mat > 0)
           {
             mat_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
             svalue(num_mat) <- paste(" N. of match: \n", sqldf("select count(*) from vms_lb", dbname = vms_DB$db), sep = "")
             enabled(mat_b) <- TRUE
           }else{
             mat_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
             svalue(num_mat) <- "   ---"
             enabled(mat_b) <- FALSE
           }
           delete(mat_g, mat_sta)
           mat_sta <<- mat_sta_n
           add(mat_g, mat_sta)
           
           n_fis <- sqldf("select count(*) from p_fish", dbname = vms_DB$db)
           if (n_fis > 0)
           {
             fis_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
             svalue(num_fis) <- paste(" N. of points: \n", sqldf("select count(*) from p_fish where FISH = 1", dbname = vms_DB$db), sep = "")
             enabled(fis_b) <- TRUE
           }else{
             fis_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
             svalue(num_fis) <- "   ---"
             enabled(fis_b) <- FALSE
           }
           delete(fis_g, fis_sta)
           fis_sta <<- fis_sta_n
           add(fis_g, fis_sta)
           ##################################
           
           enabled(sel_vms) <- TRUE
           
         })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = vms_db_f,
         handler = function(h,...){
           vms_DB$db <- ""
           #################
           
           png_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
           svalue(min_date_l) <- "   ---"
           svalue(max_date_l) <- "   ---"
           svalue(num_vess) <- "   ---"
           svalue(num_png) <- "   ---"
           enabled(png_b) <- FALSE
           delete(png_g, png_sta)
           png_sta <<- png_sta_n
           add(png_g, png_sta)
           
           wrn_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
           svalue(num_wrn) <- "   ---"
           enabled(wrn_b) <- FALSE
           delete(wrn_g, wrn_sta)
           wrn_sta <<- wrn_sta_n
           add(wrn_g, wrn_sta)
           
           trk_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
           svalue(mean_trk) <- "   ---"
           svalue(tot_trk) <- "   ---"
           enabled(trk_b) <- FALSE
           delete(trk_g, trk_sta)
           trk_sta <<- trk_sta_n
           add(trk_g, trk_sta)
           
           ntr_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
           svalue(num_ntr) <- "   ---"
           enabled(ntr_b) <- FALSE
           delete(ntr_g, ntr_sta)
           ntr_sta <<- ntr_sta_n
           add(ntr_g, ntr_sta)
           
           dep_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
           svalue(num_dep) <- "   ---"
           enabled(dep_b) <- FALSE
           delete(dep_g, dep_sta)
           dep_sta <<- dep_sta_n
           add(dep_g, dep_sta)
           
           are_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
           svalue(num_are) <- "   ---"
           enabled(are_b) <- FALSE
           delete(are_g, are_sta)
           are_sta <<- are_sta_n
           add(are_g, are_sta)
           
           mat_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
           svalue(num_mat) <- "   ---"
           enabled(mat_b) <- FALSE
           delete(mat_g, mat_sta)
           mat_sta <<- mat_sta_n
           add(mat_g, mat_sta)
           
           fis_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
           svalue(num_fis) <- "   ---"
           enabled(fis_b) <- FALSE
           delete(fis_g, fis_sta)
           fis_sta <<- fis_sta_n
           add(fis_g, fis_sta)
           ###################
           enabled(sel_vms) <- FALSE
           svalue(sel_vms_f) <- "Select VMS DB file"
         })
  ################
  dat_ref <- gimage(system.file("ico/document-quick_restart.png", package="vmsbase"))
  add(vms_db_f, dat_ref)
  addHandlerClicked(dat_ref, handler = function(h,...)
  {
    if(vms_DB$db != "")
    {
      n_png <- sqldf("select count(*) from ping", dbname = vms_DB$db)
      if (n_png > 0)
      {
        png_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
        svalue(num_vess) <- paste(" N. of Vessels: \n", sqldf("select count(distinct I_NCEE) from ping", dbname = vms_DB$db), sep = "")
        svalue(num_png) <- paste(" N. of Pings: \n", sqldf("select count(*) from ping", dbname = vms_DB$db), sep = "")
        
        min_date <<- as.numeric(sqldf("select min(DATE) from ping", dbname = vms_DB$db)[,1])
        max_date <<- as.numeric(sqldf("select max(DATE) from ping", dbname = vms_DB$db)[,1])
        
        svalue(min_date_l) <- paste("First: ", days(min_date), "/", months(min_date), "/", years(min_date), sep = "")
        svalue(max_date_l) <- paste("Last: ", days(max_date), "/", months(max_date), "/", years(max_date), sep = "")
        
        enabled(png_b) <- TRUE
      }else{
        png_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
        svalue(num_vess) <- "   ---"
        svalue(num_png) <- "   ---"
        enabled(png_b) <- FALSE
      }
      delete(png_g, png_sta)
      png_sta <<- png_sta_n
      add(png_g, png_sta)
      
      n_wrn <- sqldf("select count(*) from warn", dbname = vms_DB$db)
      if (n_wrn > 0)
      {
        wrn_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
        svalue(num_wrn) <- paste(" N. of Warnings: \n", sqldf("select count(*) from warn", dbname = vms_DB$db), sep = "")
        enabled(wrn_b) <- TRUE
      }else{
        wrn_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
        svalue(num_wrn) <- "   ---"
        enabled(wrn_b) <- FALSE
      }
      delete(wrn_g, wrn_sta)
      wrn_sta <<- wrn_sta_n
      add(wrn_g, wrn_sta)
      
      n_trk <- sqldf("select count(*) from track", dbname = vms_DB$db)
      if (n_trk > 0)
      {
        trk_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
        svalue(mean_trk) <- paste(" Mean tracks for vessel: \n", round(sqldf("select count(*) from (select distinct I_NCEE, T_NUM from track)", dbname = vms_DB$db)/sqldf("select count(distinct I_NCEE) from track", dbname = vms_DB$db), 1), sep = "")
        svalue(tot_trk) <- paste(" N. of tracks: \n", sqldf("select count(*) from (select distinct I_NCEE, T_NUM from track)", dbname = vms_DB$db), sep = "")
        enabled(trk_b) <- TRUE
        enabled(sel_vms) <- TRUE  
      }else{
        trk_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
        svalue(mean_trk) <- "   ---"
        svalue(tot_trk) <- "   ---"
        enabled(trk_b) <- FALSE
      }
      delete(trk_g, trk_sta)
      trk_sta <<- trk_sta_n
      add(trk_g, trk_sta)
      
      n_ntr <- sqldf("select count(*) from intrp", dbname = vms_DB$db)
      if (n_ntr > 0)
      {
        ntr_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
        svalue(num_ntr) <- paste(" N. of Pings: \n", sqldf("select count(*) from intrp", dbname = vms_DB$db), sep = "")
        enabled(ntr_b) <- TRUE
      }else{
        ntr_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
        svalue(num_ntr) <- "   ---"
        enabled(ntr_b) <- FALSE
      }
      delete(ntr_g, ntr_sta)
      ntr_sta <<- ntr_sta_n
      add(ntr_g, ntr_sta)
      
      n_dep <- sqldf("select count(*) from p_depth", dbname = vms_DB$db)
      if (n_dep > 0)
      {
        dep_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
        svalue(num_dep) <- paste(" N. of Pings: \n", sqldf("select count(*) from p_depth", dbname = vms_DB$db), sep = "")
        enabled(dep_b) <- TRUE
      }else{
        dep_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
        svalue(num_dep) <- "   ---"
        enabled(dep_b) <- FALSE
      }
      delete(dep_g, dep_sta)
      dep_sta <<- dep_sta_n
      add(dep_g, dep_sta)
      
      n_are <- sqldf("select count(*) from p_area", dbname = vms_DB$db)
      if (n_are > 0)
      {
        are_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
        svalue(num_are) <- paste(" N. of Pings: \n", sqldf("select count(*) from p_area", dbname = vms_DB$db), sep = "")
        enabled(are_b) <- TRUE
      }else{
        are_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
        svalue(num_are) <- "   ---"
        enabled(are_b) <- FALSE
      }
      delete(are_g, are_sta)
      are_sta <<- are_sta_n
      add(are_g, are_sta)
      
      n_mat <- sqldf("select count(*) from vms_lb", dbname = vms_DB$db)
      if (n_mat > 0)
      {
        mat_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
        svalue(num_mat) <- paste(" N. of match: \n", sqldf("select count(*) from vms_lb", dbname = vms_DB$db), sep = "")
        enabled(mat_b) <- TRUE
      }else{
        mat_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
        svalue(num_mat) <- "   ---"
        enabled(mat_b) <- FALSE
      }
      delete(mat_g, mat_sta)
      mat_sta <<- mat_sta_n
      add(mat_g, mat_sta)
      
      n_fis <- sqldf("select count(*) from p_fish", dbname = vms_DB$db)
      if (n_fis > 0)
      {
        fis_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
        svalue(num_fis) <- paste(" N. of points: \n", sqldf("select count(*) from p_fish where FISH = 1", dbname = vms_DB$db), sep = "")
        enabled(fis_b) <- TRUE
      }else{
        fis_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
        svalue(num_fis) <- "   ---"
        enabled(fis_b) <- FALSE
      }
      delete(fis_g, fis_sta)
      fis_sta <<- fis_sta_n
      add(fis_g, fis_sta)
    }
  })
  addSpring(vms_db_f)
  addSpring(chk_g3)  
  addSpring(top_g)  
  png_g <- gframe(text = "Pings", horizontal = FALSE, container = top_g)
  min_date_l <- glabel("", container = png_g)
  max_date_l <- glabel("", container = png_g)  
  num_vess <- glabel("   ---", container = png_g)
  num_png <- glabel("   ---", container = png_g)
  addSpring(png_g)
  png_b <- gbutton(text = "Show data", container = png_g, handler = function(h,..){
    pings <- gtable(sqldf("select * from ping order by DATE", dbname = vms_DB$db), expand = T)
    add(data_nb, pings, label = "Raw Pings")
  })
  png_sta <- gimage(system.file("ico/user-invisible.png", package="vmsbase"))
  add(png_g, png_sta)
  addSpring(top_g)
  enabled(png_b) <- FALSE
  
  wrn_g <- gframe(text = "Warnings", horizontal = F, container = top_g)
  num_wrn <- glabel("   ---", container = wrn_g)
  addSpring(wrn_g)
  wrn_b <- gbutton(text = "Show data", container = wrn_g, handler = function(h,..){
    wrns <- gtable(sqldf("select * from warn", dbname = vms_DB$db), expand = T)
    add(data_nb, wrns, label = "Warnings")
  })
  wrn_sta <- gimage(system.file("ico/user-invisible.png", package="vmsbase"))
  add(wrn_g, wrn_sta)
  addSpring(top_g)
  enabled(wrn_b) <- FALSE
  
  trk_g <- gframe(text = "Track", horizontal = F, container = top_g)
  
  mean_trk <- glabel("   ---", container = trk_g)
  tot_trk <- glabel("   ---", container = trk_g)
  addSpring(trk_g)
  trk_b <- gbutton(text = "Show data", container = trk_g, handler = function(h,...){
    tracks <- gtable(sqldf("select * from track order by I_NCEE", dbname = vms_DB$db), expand = T)
    add(data_nb, tracks, label = "Tracks")
  })
  trk_sta <- gimage(system.file("ico/user-invisible.png", package="vmsbase"))
  add(trk_g, trk_sta)
  addSpring(top_g)
  enabled(trk_b) <- FALSE
  
  ntr_g <- gframe(text = "Interpolated", horizontal = F, container = top_g)
  
  num_ntr <- glabel("   ---", container = ntr_g)
  addSpring(ntr_g)
  ntr_b <- gbutton(text = "Show data", container = ntr_g, handler = function(h,...){
    intrp <- gtable(sqldf("select * from intrp order by I_NCEE", dbname = vms_DB$db), expand = T)
    add(data_nb, intrp, label = "Interpolated")
  })
  ntr_sta <- gimage(system.file("ico/user-invisible.png", package="vmsbase"))
  add(ntr_g, ntr_sta)
  addSpring(top_g)
  enabled(ntr_b) <- FALSE
  
  dep_g <- gframe(text = "Depth", horizontal = F, container = top_g)
  num_dep <- glabel("   ---", container = dep_g)
  addSpring(dep_g)
  dep_b <- gbutton(text = "Show data", container = dep_g, handler = function(h,..){
    deps <- gtable(sqldf("select * from p_depth", dbname = vms_DB$db), expand = T)
    add(data_nb, deps, label = "Depths")
  })
  dep_sta <- gimage(system.file("ico/user-invisible.png", package="vmsbase"))
  add(dep_g, dep_sta)
  addSpring(top_g)
  enabled(dep_b) <- FALSE
  
  are_g <- gframe(text = "Area", horizontal = F, container = top_g)
  num_are <- glabel("   ---", container = are_g)
  addSpring(are_g)
  are_b <- gbutton(text = "Show data", container = are_g, handler = function(h,..){
    ares <- gtable(sqldf("select * from p_area", dbname = vms_DB$db), expand = T)
    add(data_nb, ares, label = "Area")
  })
  are_sta <- gimage(system.file("ico/user-invisible.png", package="vmsbase"))
  add(are_g, are_sta)
  addSpring(top_g)
  enabled(are_b) <- FALSE
  
  mat_g <- gframe(text = "VMS-LogBook Match", horizontal = F, container = top_g)
  num_mat <- glabel("   ---", container = mat_g)
  addSpring(mat_g)
  mat_b <- gbutton(text = "Show data", container = mat_g, handler = function(h,..){
    mats <- gtable(sqldf("select * from vms_lb", dbname = vms_DB$db), expand = T)
    add(data_nb, mats, label = "VMS-LB Match")
  })
  mat_sta <- gimage(system.file("ico/user-invisible.png", package="vmsbase"))
  add(mat_g, mat_sta)
  addSpring(top_g)
  enabled(mat_b) <- FALSE
  
  fis_g <- gframe(text = "Fishing Points", horizontal = F, container = top_g)
  num_fis <- glabel("   ---", container = fis_g)
  addSpring(fis_g)
  fis_b <- gbutton(text = "Show data", container = fis_g, handler = function(h,..){
    fiss <- gtable(sqldf("select * from p_fish", dbname = vms_DB$db), expand = T)
    add(data_nb, fiss, label = "Fishing Points")
  })
  fis_sta <- gimage(system.file("ico/user-invisible.png", package="vmsbase"))
  add(fis_g, fis_sta)
  addSpring(top_g)
  enabled(fis_b) <- FALSE
  
  ####################
  addSpring(sup_g)
  sel_f <- gframe(text = "Data Extraction", horizontal = TRUE, container = sup_g, expand = F)
  addSpring(sel_f)
  sel_vms <- gbutton(text = "Select on\nVMS DB", container = sel_f, handler = function(h,..)
  {
    gui_vms_db_sel(vms_DB$db)
  })
  addSpring(sel_f)
  
  bot_g <- ggroup(container = big_g, expand = T)
  data_nb <- gnotebook(tab.pos = 3, closebuttons = TRUE, container = bot_g, expand = T)
  
  enabled(sel_vms) <- FALSE
  visible(vms_select_win) <- TRUE
  
  ##################
  
  if(vms_DB$db != "") 
  {
#     svalue(sel_vms_f) <- strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])]
    svalue(sel_vms_f) <- ifelse(.Platform$OS.type == "windows", strsplit(vms_DB$db, "\\\\")[[1]][length(strsplit(vms_DB$db, "\\\\")[[1]])],strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])])
    n_png <- sqldf("select count(*) from ping", dbname = vms_DB$db)
    if (n_png > 0)
    {
      png_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
      svalue(num_vess) <- paste(" N. of Vessels: \n", sqldf("select count(distinct I_NCEE) from ping", dbname = vms_DB$db), sep = "")
      svalue(num_png) <- paste(" N. of Pings: \n", sqldf("select count(*) from ping", dbname = vms_DB$db), sep = "")
      
      min_date <- as.numeric(sqldf("select min(DATE) from ping", dbname = vms_DB$db)[,1])
      max_date <- as.numeric(sqldf("select max(DATE) from ping", dbname = vms_DB$db)[,1])
      
      svalue(min_date_l) <- paste("First: ", days(min_date), "/", months(min_date), "/", years(min_date), sep = "")
      svalue(max_date_l) <- paste("Last: ", days(max_date), "/", months(max_date), "/", years(max_date), sep = "")
      
      enabled(png_b) <- TRUE
    }else{
      png_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
      svalue(num_vess) <- "   ---"
      svalue(num_png) <- "   ---"
      enabled(png_b) <- FALSE
    }
    delete(png_g, png_sta)
    png_sta <- png_sta_n
    add(png_g, png_sta)
    
    n_wrn <- sqldf("select count(*) from warn", dbname = vms_DB$db)
    if (n_wrn > 0)
    {
      wrn_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
      svalue(num_wrn) <- paste(" N. of Warnings: \n", sqldf("select count(*) from warn", dbname = vms_DB$db), sep = "")
      enabled(wrn_b) <- TRUE
    }else{
      wrn_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
      svalue(num_wrn) <- "   ---"
      enabled(wrn_b) <- FALSE
    }
    delete(wrn_g, wrn_sta)
    wrn_sta <- wrn_sta_n
    add(wrn_g, wrn_sta)
    
    n_trk <- sqldf("select count(*) from track", dbname = vms_DB$db)
    if (n_trk > 0)
    {
      trk_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
      svalue(mean_trk) <- paste(" Mean tracks for vessel: \n", round(sqldf("select count(*) from (select distinct I_NCEE, T_NUM from track)", dbname = vms_DB$db)/sqldf("select count(distinct I_NCEE) from track", dbname = vms_DB$db), 1), sep = "")
      svalue(tot_trk) <- paste(" N. of tracks: \n", sqldf("select count(*) from (select distinct I_NCEE, T_NUM from track)", dbname = vms_DB$db), sep = "")
      enabled(trk_b) <- TRUE
      enabled(sel_vms) <- TRUE  
    }else{
      trk_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
      svalue(mean_trk) <- "   ---"
      svalue(tot_trk) <- "   ---"
      enabled(trk_b) <- FALSE
    }
    delete(trk_g, trk_sta)
    trk_sta <- trk_sta_n
    add(trk_g, trk_sta)
    
    n_ntr <- sqldf("select count(*) from intrp", dbname = vms_DB$db)
    if (n_ntr > 0)
    {
      ntr_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
      svalue(num_ntr) <- paste(" N. of Pings: \n", sqldf("select count(*) from intrp", dbname = vms_DB$db), sep = "")
      enabled(ntr_b) <- TRUE
    }else{
      ntr_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
      svalue(num_ntr) <- "   ---"
      enabled(ntr_b) <- FALSE
    }
    delete(ntr_g, ntr_sta)
    ntr_sta <- ntr_sta_n
    add(ntr_g, ntr_sta)
    
    n_dep <- sqldf("select count(*) from p_depth", dbname = vms_DB$db)
    if (n_dep > 0)
    {
      dep_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
      svalue(num_dep) <- paste(" N. of Pings: \n", sqldf("select count(*) from p_depth", dbname = vms_DB$db), sep = "")
      enabled(dep_b) <- TRUE
    }else{
      dep_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
      svalue(num_dep) <- "   ---"
      enabled(dep_b) <- FALSE
    }
    delete(dep_g, dep_sta)
    dep_sta <- dep_sta_n
    add(dep_g, dep_sta)
    
    n_are <- sqldf("select count(*) from p_area", dbname = vms_DB$db)
    if (n_are > 0)
    {
      are_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
      svalue(num_are) <- paste(" N. of Pings: \n", sqldf("select count(*) from p_area", dbname = vms_DB$db), sep = "")
      enabled(are_b) <- TRUE
    }else{
      are_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
      svalue(num_are) <- "   ---"
      enabled(are_b) <- FALSE
    }
    delete(are_g, are_sta)
    are_sta <- are_sta_n
    add(are_g, are_sta)
    
    n_mat <- sqldf("select count(*) from vms_lb", dbname = vms_DB$db)
    if (n_mat > 0)
    {
      mat_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
      svalue(num_mat) <- paste(" N. of match: \n", sqldf("select count(*) from vms_lb", dbname = vms_DB$db), sep = "")
      enabled(mat_b) <- TRUE
    }else{
      mat_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
      svalue(num_mat) <- "   ---"
      enabled(mat_b) <- FALSE
    }
    delete(mat_g, mat_sta)
    mat_sta <- mat_sta_n
    add(mat_g, mat_sta)
    
    n_fis <- sqldf("select count(*) from p_fish", dbname = vms_DB$db)
    if (n_fis > 0)
    {
      fis_sta_n <- gimage(system.file("ico/user-available.png", package="vmsbase"))
      svalue(num_fis) <- paste(" N. of points: \n", sqldf("select count(*) from p_fish where FISH = 1", dbname = vms_DB$db), sep = "")
      enabled(fis_b) <- TRUE
    }else{
      fis_sta_n <- gimage(system.file("ico/user-busy.png", package="vmsbase"))
      svalue(num_fis) <- "   ---"
      enabled(fis_b) <- FALSE
    }
    delete(fis_g, fis_sta)
    fis_sta <- fis_sta_n
    add(fis_g, fis_sta)
    
    enabled(sel_vms) <- TRUE
  }
}
