#' VMSbase main GUI
#' 
#' The \code{gui_main} function implements the main graphical user interface of the
#'  VMSbase package
#' 
#' With this gui the user can access all the functionalities of VMSbase
#' 
#' @return This function does not return a value. 
#' 
#' @usage gui_main()
#' 
#' @export gui_main

gui_main <- function ()
{
  vms_db_file <- ""
  lb_db_file <- ""
  themap_file <- ""
  harb_file <- ""
  bathy_file <- ""
  
  main_win <- gwindow(paste("VMS Base - Version ", packageVersion("vmsbase"), sep = "") , visible = FALSE)
  
  theg <- ggroup(horizontal = TRUE, container = main_win)
  
  left <- ggroup(horizontal = FALSE, container = theg)
  gimage(system.file("img/vms_300_spada_ais.png", package = "vmsbase"), container = left)
  gseparator(container = left)
  
  addSpring(left)
  
  wo_sp_fr <- gframe(text = "Project Management", horizontal = FALSE, container = left)
  ##VMS DataBase file
  vms_db_f <- gframe(text = "VMS DB file", horizontal = TRUE, container = wo_sp_fr)
  addSpring(vms_db_f)
  sel_vms_f <- glabel("Select VMS DB file", container = vms_db_f)
  addSpring(vms_db_f)
  gimage(system.file("ico/folder-blue.png", package="vmsbase"), container = vms_db_f,
         handler = function(h,...){
           vms_db_file <<- gfile(text = "Select VMS DataBase file",
                                 type = "open",
                                 filter = list("VMS DB file" = list(patterns = c("*.vms.sqlite"))))
           svalue(sel_vms_f) <- ifelse(.Platform$OS.type == "windows", strsplit(vms_db_file, "\\\\")[[1]][length(strsplit(vms_db_file, "\\\\")[[1]])],strsplit(vms_db_file, "/")[[1]][length(strsplit(vms_db_file, "/")[[1]])])
         })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = vms_db_f,
         handler = function(h,...){
           vms_db_file <<- ""
           svalue(sel_vms_f) <- "Select VMS DB file"
         })
  
  ##LB DataBase file
  lb_db_f <- gframe(text = "LogBook DB file", horizontal = TRUE, container = wo_sp_fr)
  addSpring(lb_db_f)
  sel_lb_f <- glabel("Select LB DB file", container = lb_db_f)
  addSpring(lb_db_f)
  gimage(system.file("ico/folder-orange.png", package="vmsbase"), container = lb_db_f,
         handler = function(h,...){
           lb_db_file <<- gfile(text = "Select LB DataBase file",
                                type = "open",
                                filter = list("LB DB file" = list(patterns = c("*.lb.sqlite"))))
           svalue(sel_lb_f) <- ifelse(.Platform$OS.type == "windows", strsplit(lb_db_file, "\\\\")[[1]][length(strsplit(lb_db_file, "\\\\")[[1]])],strsplit(lb_db_file, "/")[[1]][length(strsplit(lb_db_file, "/")[[1]])])
         })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = lb_db_f,
         handler = function(h,...){
           lb_db_file <<- ""
           svalue(sel_lb_f) <- "Select LB DB file"
         })
  
  ##Land file
  cus_map_g <- gframe(text = "Land Shape File", horizontal = TRUE, container = wo_sp_fr)
  addSpring(cus_map_g)
  cus_map_lab <- glabel("Select Land Shape File", container = cus_map_g)
  addSpring(cus_map_g)
  gimage(system.file("ico/folder-html.png", package="vmsbase"), container = cus_map_g,
         handler = function(h,...){
           themap_file <<- gfile(text = "Select Land ShapePoly map",
                                 type = "open",
                                 filter = list("shp data" = list(patterns = c("*.shp"))))
           svalue(cus_map_lab) <- paste("Land: ", ifelse(.Platform$OS.type == "windows", strsplit(themap_file, "\\\\")[[1]][length(strsplit(themap_file, "\\\\")[[1]])],strsplit(themap_file, "/")[[1]][length(strsplit(themap_file, "/")[[1]])]), sep = "")         })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = cus_map_g,
         handler = function(h,...){
           themap_file <<- ""
           svalue(cus_map_lab) <- "Select Land Shape File"
         })
  
  ##Harbours file
  cus_har_g <- gframe(text = "Harbours Shape File", horizontal = TRUE, container = wo_sp_fr)
  addSpring(cus_har_g)
  cus_har_lab <- glabel("Select Harbours Shape File", container = cus_har_g)
  addSpring(cus_har_g)
  gimage(system.file("ico/folder-man.png", package="vmsbase"), container = cus_har_g,
         handler = function(h,...){
           harb_file <<- gfile(text = "Select ShapePoints map",
                               type = "open",
                               filter = list("shp data" = list(patterns = c("*.shp"))))
           svalue(cus_har_lab) <- paste("Harbour: ", ifelse(.Platform$OS.type == "windows", strsplit(harb_file, "\\\\")[[1]][length(strsplit(harb_file, "\\\\")[[1]])],strsplit(harb_file, "/")[[1]][length(strsplit(harb_file, "/")[[1]])]), sep = "")          })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = cus_har_g,
         handler = function(h,...){
           harb_file <<- ""
           svalue(cus_har_lab) <- "Select Harbours Shape File"
         })
  
  ##Depth file
  cus_dep_g <- gframe(text = "Bathymetry File", horizontal = TRUE, container = wo_sp_fr)
  addSpring(cus_dep_g)
  cus_dep_lab <- glabel("Select Bathymetry File", container = cus_dep_g)
  addSpring(cus_dep_g)
  gimage(system.file("ico/folder-download.png", package="vmsbase"), container = cus_dep_g,
         handler = function(h,...){
           bathy_file <<- gfile(text = "Select Bathymetry File",
                                type = "open",
                                filter = list("bathy data" = list(patterns = c("*sqlitebathy.rData"))))
           svalue(cus_dep_lab) <- paste("Depth: ", ifelse(.Platform$OS.type == "windows", strsplit(bathy_file, "\\\\")[[1]][length(strsplit(bathy_file, "\\\\")[[1]])],strsplit(bathy_file, "/")[[1]][length(strsplit(bathy_file, "/")[[1]])]), sep = "")
         })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = cus_dep_g,
         handler = function(h,...){
           bathy_file <<- ""
           svalue(cus_dep_lab) <- "Select Bathymetry File"
         })
  
  sa_lo_ws_g <- ggroup(horizontal = TRUE, container = wo_sp_fr)
  addSpring(sa_lo_ws_g)
  save_ws_b <-gbutton(text = "Save\nWorkSpace", container = sa_lo_ws_g, handler = function(h,...)
  {
    ws_file <- ""
    ws_file <- gfile(text = "Enter the WorkSpace name...",
                     initialfilename = "*.ws.rData",
                     type = "save",
                     filter = list("R WorkSpace" = list(patterns = c("*ws.rData"))))
    if(!is.na(ws_file))
    {
      cat(ls(name = parent.env(environment()), pattern = "file"), sep = " - ")
      save(list = ls(name = parent.env(environment()), pattern = "file"), file = ws_file)
      gmessage(paste("WorkSpace Saved!\n\nFile: ", ws_file, sep = ""), title="WorkSpace Saved",
               icon = "info") 
    }else{
      gmessage(paste("WorkSpace not saved!\n\nPlease enter a valid name and\ndestination of the WorkSpace file", sep = ""), title="WorkSpace Save Error!",
               icon = "error")
    }
  })
  addSpring(sa_lo_ws_g)
  load_ws_b <-gbutton(text = "Load\nWorkSpace", container = sa_lo_ws_g, handler = function(h,...)
  {
    ws_file <- ""
    ws_file <- gfile(text = "Select a WorkSpace file...",
                     type = "open",
                     filter = list("R WorkSpace" = list(patterns = c("*ws.rData"))))
    
    if(!is.na(ws_file))
    {
      load(file = ws_file, envir = parent.env(environment()))
      
      ifelse(harb_file != "", svalue(cus_har_lab) <- paste("Harbour: ", ifelse(.Platform$OS.type == "windows", strsplit(harb_file, "\\\\")[[1]][length(strsplit(harb_file, "\\\\")[[1]])],strsplit(harb_file, "/")[[1]][length(strsplit(harb_file, "/")[[1]])]), sep = ""), svalue(cus_har_lab) <- "Select Harbours Shape File")
      ifelse(bathy_file != "",  svalue(cus_dep_lab) <- paste("Depth: ", ifelse(.Platform$OS.type == "windows", strsplit(bathy_file, "\\\\")[[1]][length(strsplit(bathy_file, "\\\\")[[1]])],strsplit(bathy_file, "/")[[1]][length(strsplit(bathy_file, "/")[[1]])]), sep = ""), svalue(cus_dep_lab) <- "Select Bathymetry File")
      ifelse(themap_file != "", svalue(cus_map_lab) <- paste("Land: ", ifelse(.Platform$OS.type == "windows", strsplit(themap_file, "\\\\")[[1]][length(strsplit(themap_file, "\\\\")[[1]])],strsplit(themap_file, "/")[[1]][length(strsplit(themap_file, "/")[[1]])]), sep = ""), svalue(cus_map_lab) <- "Select Land Shape File")
      ifelse(lb_db_file != "", svalue(sel_lb_f) <- ifelse(.Platform$OS.type == "windows", strsplit(lb_db_file, "\\\\")[[1]][length(strsplit(lb_db_file, "\\\\")[[1]])],strsplit(lb_db_file, "/")[[1]][length(strsplit(lb_db_file, "/")[[1]])]), svalue(sel_lb_f) <- "Select LB DB file")
      ifelse(vms_db_file != "", svalue(sel_vms_f) <- ifelse(.Platform$OS.type == "windows", strsplit(vms_db_file, "\\\\")[[1]][length(strsplit(vms_db_file, "\\\\")[[1]])],strsplit(vms_db_file, "/")[[1]][length(strsplit(vms_db_file, "/")[[1]])]), svalue(sel_vms_f) <- "Select VMS DB file")
      
      gmessage(paste("WorkSpace Loaded!\n\nFile: ", ws_file, sep = ""), title="WorkSpace Loaded",
               icon = "info") 
    }else{
      gmessage(paste("WorkSpace not loaded!\n\nPlease select a valid WorkSpace file", sep = ""), title="WorkSpace Load Error!",
               icon = "error")
    }
  })
  addSpring(sa_lo_ws_g)
  
  addSpace(theg, 20)
  right <- ggroup(horizontal = TRUE, container = theg)
  
  ### VMS
  top <- ggroup(horizontal = FALSE, container = right)
  
  # VMS DATA MANAGER
  top_left <- gframe(text = "VMS Data Management", horizontal = FALSE, container = top)
  top_tl <- ggroup(horizontal = TRUE, container = top_left)
  mid_tl <- ggroup(horizontal = TRUE, container = top_left)
  bot_tl <- ggroup(horizontal = TRUE, container = top_left)
  
  # top
  addSpring(top_tl)
  vms_raw <- ggroup(horizontal = FALSE, container = top_tl, spacing = 0)
  gimage(system.file("ico/edit-4.png", package="vmsbase"), container = vms_raw)
  gbutton("Edit Raw\n    Data", container = vms_raw, handler = function(h,...)
  {gui_vms_editraw()})
  addSpring(top_tl)
  vms_db <- ggroup(horizontal = FALSE, container = top_tl, spacing = 0)
  gimage(system.file("ico/db_add.png", package="vmsbase"), container = vms_db)
  gbutton("   Create\nDataBase", container = vms_db, handler = function(h,...)
  {gui_vms_db_create()})
  addSpring(top_tl)
  
  vms_db <- ggroup(horizontal = FALSE, container = top_tl, spacing = 0)
  gimage(system.file("ico/db_mix_sou.ico", package="vmsbase"), container = vms_db)
  gbutton("   Merge\nSources", container = vms_db, handler = function(h,...)
  {gui_vmsdb_mixsou()})
  addSpring(top_tl)
  
  data_vie <- ggroup(horizontal = FALSE, container = top_tl, spacing = 0)
  gimage(system.file("ico/folder-saved-search.png", package="vmsbase"), container = data_vie)
  gbutton("VMS Data\n   Viewer", container = data_vie, handler = function(h,...)
  {gui_vms_db_stat(vms_db_name = vms_db_file)})
  addSpring(top_tl)
  
  # mid
  addSpring(mid_tl)
  vms_clean <- ggroup(horizontal = FALSE, container = mid_tl, spacing = 0)
  gimage(system.file("ico/edit-clear-2.png", package="vmsbase"), container = vms_clean)
  gbutton("Clean DB\n    Data", container = vms_clean, handler = function(h,...)
  {gui_vms_db_clean(vms_db_name = vms_db_file, map_file_name = themap_file, harb_file_name = harb_file)})
  addSpring(mid_tl)
  vms_cut <- ggroup(horizontal = FALSE, container = mid_tl, spacing = 0)
  gimage(system.file("ico/retroshare1.png", package="vmsbase"), container = vms_cut)
  gbutton("  Track\nCutting", container = vms_cut, handler = function(h,...)
  {gui_vms_db_cut(vms_db_name = vms_db_file)})
  addSpring(mid_tl)
  vms_intr <- ggroup(horizontal = FALSE, container = mid_tl, spacing = 0)
  gimage(system.file("ico/draw-bezier-curves.png", package="vmsbase"), container = vms_intr)
  gbutton("        Track\nInterpolation", container = vms_intr, handler = function(h,...)
  {gui_vms_db_intr(vms_db_name = vms_db_file)})
  addSpring(mid_tl)
  
  ##########################
  # bot
  addSpring(bot_tl)
  vms_dep <- ggroup(horizontal = FALSE, container = bot_tl, spacing = 0)
  gimage(system.file("ico/go-down-3.png", package="vmsbase"), container = vms_dep)
  gbutton("      Assign\nBathymetry", container = vms_dep, handler = function(h,...)
  {gui_vms_db_dep(vms_db_name = vms_db_file)})
  addSpring(bot_tl)
  vms_are <- ggroup(horizontal = FALSE, container = bot_tl, spacing = 0)
  gimage(system.file("ico/map-compass.png", package="vmsbase"), container = vms_are)
  gbutton("Assign\n  Area", container = vms_are, handler = function(h,...)
  {gui_vms_db_are(vms_db_name = vms_db_file)})
  addSpring(bot_tl)
  vms_savdep <- ggroup(horizontal = FALSE, container = bot_tl, spacing = 0)
  gimage(system.file("ico/download-3.png", package="vmsbase"), container = vms_savdep)
  gbutton("     Get\nIsoBaths", container = vms_savdep, handler = function(h,...)
  {gui_vms_save_bat(vms_db_name = vms_db_file)})
  addSpring(bot_tl)
  
  ####
  addSpring(top)
  
  # VMS VIEWER
  mid_left <- gframe(text = "VMS Data Visualization", horizontal = FALSE, container = top)
  top_ml <- ggroup(horizontal = TRUE, container = mid_left)
  
  addSpring(top_ml)
  vms_viz_ping <- ggroup(horizontal = FALSE, container = top_ml, spacing = 0)
  gimage(system.file("ico/system-wifi.png", package="vmsbase"), container = vms_viz_ping)
  gbutton("  Ping\nViewer", container = vms_viz_ping, handler = function(h,...)
  {gui_vms_view_ping(vms_db_name = vms_db_file, bathy_file_name = bathy_file)})
  addSpring(top_ml)
  vms_viz_track <- ggroup(horizontal = FALSE, container = top_ml, spacing = 0)
  gimage(system.file("ico/retroshare1.png", package="vmsbase"), container = vms_viz_track)
  gbutton(" Track\nViewer", container = vms_viz_track, handler = function(h,...)
  {gui_vms_view_track(vms_db_name = vms_db_file, bathy_file_name = bathy_file)})
  addSpring(top_ml)
  vms_viz_ntr <- ggroup(horizontal = FALSE, container = top_ml, spacing = 0)
  gimage(system.file("ico/draw-bezier-curves.png", package="vmsbase"), container = vms_viz_ntr)
  gbutton("Interpolation\n      Viewer", container = vms_viz_ntr, handler = function(h,...)
  {gui_vms_view_intrp(vms_db_name = vms_db_file, bathy_file_name = bathy_file)})
  addSpring(top_ml)
  
  vms_viz_adv <- ggroup(horizontal = FALSE, container = top_ml, spacing = 0)
  gimage(system.file("ico/adva-vie.png", package="vmsbase"), container = vms_viz_adv)
  vms_viz_adv <- gbutton("Google\nViewer", container = vms_viz_adv, handler = function(h,...)
  {gui_vms_viz_adv(vms_db_name = vms_db_file)})
  addSpring(top_ml)
  
  addSpring(top)
  
  ### LOGBOOK
  top2 <- ggroup(horizontal = FALSE, container = right)
  
  # LB DATA MANAGER
  top_rig <- gframe(text = "LogBook Data Management", horizontal = TRUE, container = top2)
  addSpring(top_rig)
  lb_raw <- ggroup(horizontal = FALSE, container = top_rig, spacing = 0)
  gimage(system.file("ico/edit-4.png", package="vmsbase"), container = lb_raw)
  gbutton("Edit Raw\n    Data", container = lb_raw, handler = function(h,...)
  {gui_lb_editraw()})
  addSpring(top_rig)
  lb_db <- ggroup(horizontal = FALSE, container = top_rig, spacing = 0)
  gimage(system.file("ico/db_add.png", package="vmsbase"), container = lb_db)
  gbutton("   Create\nDataBase", container = lb_db, handler = function(h,...)
  {gui_lb_db_create()})
  addSpring(top_rig)
  lb_data_vie <- ggroup(horizontal = FALSE, container = top_rig, spacing = 0)
  gimage(system.file("ico/folder-saved-search.png", package="vmsbase"), container = lb_data_vie)
  gbutton("LogBook Data\n   Viewer", container = lb_data_vie, handler = function(h,...)
  {gui_lb_db_stat(lb_db_name = lb_db_file)})

  addSpring(top_rig)
  
  #LB EDIT
  mid_rig <- gframe(text = "LogBook Data Analysis", horizontal = TRUE, container = top2)
  addSpring(mid_rig)
  lb_met_dis <- ggroup(horizontal = FALSE, container = mid_rig, spacing = 0)
  gimage(system.file("ico/office-chart-ring.png", package="vmsbase"), container = lb_met_dis)
  gbutton("    Metier\nDiscovery", container = lb_met_dis, handler = function(h,...)
  {gui_lb_met_dis(lb_db_name = lb_db_file)})
  addSpring(mid_rig)
  lb_met_edi <- ggroup(horizontal = FALSE, container = mid_rig, spacing = 0)
  gimage(system.file("ico/preferences-desktop-2.png", package="vmsbase"), container = lb_met_edi)
  gbutton("Metier\nEditing", container = lb_met_edi, handler = function(h,...)
  {gui_lb_met_edi()})
  addSpring(mid_rig)
  lb_met_cla <- ggroup(horizontal = FALSE, container = mid_rig, spacing = 0)
  gimage(system.file("ico/office-chart-pie.png", package="vmsbase"), container = lb_met_cla)
  gbutton("        Metier\nClassification", container = lb_met_cla, handler = function(h,...)
  {gui_lb_met_cla(lb_db_name = lb_db_file)})
  addSpring(mid_rig)
  
  
  bot_rig <- gframe(text = "VMS-LogBook Analysis",horizontal = TRUE, container = top2)
  addSpring(bot_rig)
  lb_vms_join <- ggroup(horizontal = FALSE, container = bot_rig, spacing = 0)
  gimage(system.file("ico/referencer.png", package="vmsbase"), container = lb_vms_join)
  gbutton("LogBook-VMS\n     Matching", container = lb_vms_join, handler = function(h,...)
  {gui_join_lb_vms(lb_db_name = lb_db_file, vms_db_name = vms_db_file)})
  addSpring(bot_rig)
  
  lb_vms_pred <- ggroup(horizontal = FALSE, container = bot_rig, spacing = 0)
  gimage(system.file("ico/office-chart-area.png", package="vmsbase"), container = lb_vms_pred)
  gbutton("Predict\n Metier", container = lb_vms_pred, handler = function(h,...)
  {gui_vms_met_pred(vms_db_name = vms_db_file)})
  
  addSpring(bot_rig)
  lb_vms_fish <- ggroup(horizontal = FALSE, container = bot_rig, spacing = 0)
  gimage(system.file("ico/office-chart-scatter.png", package="vmsbase"), container = lb_vms_fish)
  gbutton("   Find\nFishing\n Points", container = lb_vms_fish, handler = function(h,...)
  {gui_mark_fis_poi(vms_db_name = vms_db_file, harb_file_name = harb_file)})
  addSpring(bot_rig)
  
  
  ### DATA OUTPUT
  
  bot2_rig <- gframe(text = "Data Output",horizontal = TRUE, container = top2)
  addSpring(bot2_rig)
  vms_grid <- ggroup(horizontal = FALSE, container = bot2_rig, spacing = 0)
  gimage(system.file("ico/view-calendar-timeline.png", package="vmsbase"), container = vms_grid)
  lb_vms_dcf_b <- gbutton("Gridding\n& Mapping", container = vms_grid, handler = function(h,...)
  {gui_out_grid(vms_db_name = vms_db_file)}) 
  addSpring(bot2_rig)
  lb_vms_dcf <- ggroup(horizontal = FALSE, container = bot2_rig, spacing = 0)
  gimage(system.file("ico/blockdevice-2.png", package="vmsbase"), container = lb_vms_dcf)
  lb_vms_dcf_b <- gbutton("DCF\nIndicators", container = lb_vms_dcf, handler = function(h,...)
  {gui_dcf_ind()})
  addSpring(bot2_rig)
  trawl <- ggroup(horizontal = FALSE, container = bot2_rig, spacing = 0)
  gimage(system.file("ico/system-search-5.png", package="vmsbase"), container = trawl)
  trawl_b <- gbutton("Trawled Area\nViewer", container = trawl, handler = function(h,...)
  {})
  addSpring(bot2_rig)
  
  ### About
  addSpring(top2)
  gseparator(container = top2)
  left_b <- ggroup(horizontal = FALSE, container = top2)
  torver_logo <- gimage(system.file("img/torvergata.png", package = "vmsbase"))
  add(left_b, torver_logo)
  addHandlerClicked(torver_logo, handler = function(h,...){
    browseURL(url = "http://web.uniroma2.it", browser = getOption("browser"))
  }) 
  addSpring(left_b)
  gbutton("About...", container = left_b, handler = function(h,...)
  {
    ab_di <- gbasicdialog(title="VMSBASE About...", do.buttons=FALSE)
    size(ab_di) <- c(750, 400)
    ab_lab <- gtext("
                                                                               ---    VMSbase  2.1    ---
                    
                                                                                      List of Authors 
                           =====================================================================
                    
                    VMS BASE represents the implementation of several R routines which have been developed
                    by the \"Tor Vergata\" University of Rome Team involved in the Italian National Program
                    for the Data Collection Framework for Fisheries Data between 2009-2012.                      
                    Namely, the major part of the original routines have been written by:
                    Tommaso Russo, Researcher at the Dept. of Biology
                    Antonio Parisi, Researcher at the Dept. of Economics and Finance 
                    
                    VMS BASE routines, structure, and interface have been developed in 2012 by the following 
                    people (listed with respect to the conceptual contribute):
                    
                    General Idea and Rationale: Tommaso Russo, Antonio Parisi, Stefano Cataudella
                    (Full professor of Ecology, Dept. of Biology, \"Tor Vergata\" University of Rome)
                    
                    Software realization: 
                    Lorenzo D'Andrea
                    (PhD student, Research Fellow at the \"Tor Vergata\" University of Rome),
                    Tommaso Russo, Antonio Parisi.
                    
                    Data Control, Software Testing: 
                    Lorenzo D'Andrea, Jacopo Pulcinella
                    (PhD student, Research Fellow at the \"Tor Vergata\" University of Rome).
                    
                    Other people have contributed to some VMS BASE modules:
                    Metier Discovery and Metier Classification:
                    Paolo Carpentieri (INEA, National Institute of Agricultural Economics) 
                    Fabio Fiorentino
                    [Institute for the Coastal Marine Environment (IAMC) - Italian National Research Council (CNR)]
                    Enrico Arneri (FAO - Adriamed Project)        
                    Trawled Area Viewer: 
                    Antonello Sala
                    [Institute of Marine Sciences (ISMAR) - Italian National Research Council (CNR)]
                    
                    Without any one of them the project would not exist in its current form. 
                    Please allow us to extend our most cordial thanks to all of you. 
                    
                    Appropriate reference for some routines/approaches can be found in:
                    1) Russo T., Parisi A., Cataudella S. (2011). New insights in interpolating 
                    fishing tracks from VMS data for different metiers. Fisheries Research 108: 184-194.
                    2) Russo T., Parisi A., Prorgi M., Boccoli F., Cignini I., Tordoni M.,
                    Cataudella S. (2011). When behaviour reveals activity: inferring fishing 
                    metier from VMS data by artificial neural network. Fisheries Research 111: 53-64.
                    3) Russo T., Parisi A., Cataudella S. (2013). Spatial  indicators  of  fishing
                    pressure:  Preliminary  analyses  and  possible developments. Ecological  Indicators 26, 141-153.
                    
                           =====================================================================

                             Thanks to the \"Open Icon Library Project\" for the icons used in the VMSbase GUIs
                                                                       http://openiconlibrary.sourceforge.net/
                    ",
                    container = ab_di)
    enabled(ab_lab) <- FALSE
    visible(ab_di, set = TRUE)
    
  })
  
  visible(main_win) <- TRUE
}