
#' VMS DB View Track Data GUI
#'  
#' The \code{gui_vms_view_track} function implements the graphical user interface for the
#'  VMS DB routine to view vessel and track data with the Google Maps capabilities.
#' 
#' This function, with a VMS database,
#'  plots vessel and track data on Hybrid Satellite images from Google Maps.
#'   
#' @param vms_db_name The path of a VMS DataBase
#'   
#' @return This function does not return a value. 
#' 
#' @usage gui_vms_viz_adv(vms_db_name = "")
#' 
#' @export gui_vms_viz_adv
#'
#'@seealso \code{\link{gui_vms_view_ping}} \code{\link{gui_vms_view_track}} \code{\link{gui_vms_view_intrp}}

gui_vms_viz_adv <- function (vms_db_name = "")
{
  vms_DB <- vms_DB$new()
  vms_DB$db <- vms_db_name
  
  adv_view_win <- gwindow("Advanced VMS Viewer", visible = FALSE)
  
  big_g <- ggroup(horizontal = T, container = adv_view_win, expand = T)
  
  g_lef <- ggroup(horizontal = F, container = big_g)
  #g_lef_b <- ggroup(horizontal = F, container = big_g)
  
  chk_g3 <- ggroup(horizontal = TRUE, container = g_lef)
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
           svalue(s_bar) <- "Loading DataSet"
           incee <- sqldf("select distinct I_NCEE from intrp order by I_NCEE", dbname = vms_DB$db)
           #colnames(incee) <- "Vessel ID"
           selves[] <- incee
           
           svalue(s_bar) <- "DataSet loaded! Click on a vessel ID to load tracks"
           
           enabled(d_zoom) <- TRUE
           enabled(b_draw) <- TRUE
           enabled(b_save_j) <- TRUE
           enabled(c_opt ) <- TRUE
           enabled(selves) <- TRUE
           enabled(vestra) <- TRUE
         })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = vms_db_f,
         handler = function(h,...){
           vms_DB$db <- ""
           enabled(selves) <- FALSE
           enabled(vestra) <- FALSE
           enabled(d_zoom) <- FALSE
           enabled(b_draw) <- FALSE
           enabled(b_save_j) <- FALSE
           enabled(c_opt ) <- FALSE
           svalue(sel_vms_f) <- "Select VMS DB file"
         })
  addSpring(chk_g3)
  ################
  
  d_zoom <- gdroplist(2:20, selected = 8, container = g_lef, handler = function(h,...)
  {
    svalue(s_bar) <- "Press `Draw Map` button to update view"
  })
  
  c_opt <- gcheckboxgroup(c("Points","Path"), checked = c(TRUE, TRUE), container = g_lef, handler = function(h,...)
  {
    svalue(s_bar) <- "Press `Draw Map` button to update view"
  })
  
  b_draw <- gbutton("Draw\nMap", container = g_lef, handler = function(h,...)
  {
    enabled(selves) <- FALSE
    enabled(vestra) <- FALSE
    enabled(b_draw) <- FALSE
    enabled(b_save_j) <- FALSE
    
    vnum <- svalue(selves)
    if(length(svalue(vestra)) != 0)
    {
      trackn <- svalue(vestra)
      
      s_track  <- fn$sqldf("select * from intrp where I_NCEE = `vnum` and T_NUM = `trackn` order by DATE", dbname = vms_DB$db)
      
      not_intr <- which(s_track[,"P_INT"] == 0)
      data <- as.data.frame(cbind(s_track[,"LON"], s_track[,"LAT"]))
      real_p <- as.data.frame(cbind(s_track[not_intr,"LON"], s_track[not_intr,"LAT"]))
      colnames(data) <- c("lon", "lat")
      colnames(real_p) <- c("lon", "lat")
      
      if(nrow(data) > 2)
      {
        b_box <- make_bbox(s_track[,"LON"], s_track[,"LAT"], f = 0.3)
        
        svalue(s_bar) <- paste("Loading map...", sep = "")
        
        map <- suppressWarnings(get_map(location = b_box,
                                        zoom = svalue(d_zoom), 
                                        scale = "auto",
                                        maptype = "hybrid",
                                        filename = "ggmapTemp",
                                        source = "google"))
        
        if("Path" %in% svalue(c_opt) & "Points" %in% svalue(c_opt))
        {
          print(ggmap(map) +
                  geom_point(data = data, colour = "darkorange") +  
                  geom_path(data = data, colour = "firebrick"))
        }
        
        if("Path" %in% svalue(c_opt) & !("Points" %in% svalue(c_opt)))
        {
          print(ggmap(map) +
                  geom_path(data = data, colour = "darkorange"))
        }
        
        if(!("Path" %in% svalue(c_opt)) & "Points" %in% svalue(c_opt))
        {
          print(ggmap(map) +
                  geom_point(data = data, colour = "firebrick"))
        }
        
        if(!("Path" %in% svalue(c_opt)) & !("Points" %in% svalue(c_opt)))
        {
          print(ggmap(map) +
                  geom_point(data = data, colour = "firebrick",alpha = 1/10)+  
                  geom_path(data = data, colour = "firebrick",alpha = 1/10))
        }
        
        svalue(s_bar) <- paste("Map loading complete!", sep = "")
        
      }
    }else{
      
      svalue(s_bar) <- paste("Loading Vessel ", vnum, sep = "")
      vessel <- fn$sqldf("select * from intrp where I_NCEE = `vnum`", dbname = vms_DB$db)
      data <- as.data.frame(cbind(vessel[,"LON"], vessel[,"LAT"]))
      colnames(data) <- c("lon", "lat")
      vestra[] <- unique(vessel["T_NUM"])
      
      svalue(s_bar) <- paste("Loading Vessel Tracks complete!", sep = "")
      
      if(nrow(vessel) != 0)
      {
        
        num_trk <- nrow(vessel)
        
        b_box <- make_bbox(vessel[,"LON"], vessel[,"LAT"], f = 0.3)
        
        svalue(s_bar) <- paste("Loading map...", sep = "")
        
        map <- suppressWarnings(get_map(location = b_box,
                                        zoom = svalue(d_zoom), 
                                        scale = "auto",
                                        maptype = "hybrid",
                                        filename = "ggmapTemp",
                                        source = "google"))
        
        if("Path" %in% svalue(c_opt) & "Points" %in% svalue(c_opt))
        {
          print(ggmap(map) +
                  geom_point(data = data, colour = "darkorange") +  
                  geom_path(data = data, colour = "firebrick"))
        }
        
        if("Path" %in% svalue(c_opt) & !("Points" %in% svalue(c_opt)))
        {
          print(ggmap(map) +
                  geom_path(data = data, colour = "darkorange"))
        }
        
        if(!("Path" %in% svalue(c_opt)) & "Points" %in% svalue(c_opt))
        {
          print(ggmap(map) +
                  geom_point(data = data, colour = "firebrick"))
        }
        
        if(!("Path" %in% svalue(c_opt)) & !("Points" %in% svalue(c_opt)))
        {
          print(ggmap(map) +
                  geom_point(data = data, colour = "firebrick",alpha = 1/10)+  
                  geom_path(data = data, colour = "firebrick",alpha = 1/10))
        }
        
        svalue(s_bar) <- paste("Map loading complete! Select a track!", sep = "")
        
      }
    }  
    enabled(selves) <- TRUE
    enabled(vestra) <- TRUE
    enabled(b_draw) <- TRUE
    enabled(b_save_j) <- TRUE
  })
  
  b_save_j <- gbutton("Save Map\nto Jpeg", container = g_lef, handler = function(h,...)
  {
    jpeg(filename = gfile(text = "Jpeg file path and name",
                          initialfilename = "*.jpeg",
                          type = "save"),
         width = 800, height = 800)
    
    enabled(selves) <- FALSE
    enabled(vestra) <- FALSE
    enabled(b_draw) <- FALSE
    enabled(b_save_j) <- FALSE
    
    vnum <- svalue(selves)
    if(length(svalue(vestra)) != 0)
    {
      trackn <- svalue(vestra)
      
      s_track  <- fn$sqldf("select * from intrp where I_NCEE = `vnum` and T_NUM = `trackn` order by DATE", dbname = vms_DB$db)
      
      not_intr <- which(s_track[,"P_INT"] == 0)
      data <- as.data.frame(cbind(s_track[,"LON"], s_track[,"LAT"]))
      real_p <- as.data.frame(cbind(s_track[not_intr,"LON"], s_track[not_intr,"LAT"]))
      colnames(data) <- c("lon", "lat")
      colnames(real_p) <- c("lon", "lat")
      
      if(nrow(data) > 2)
      {
        b_box <- make_bbox(s_track[,"LON"], s_track[,"LAT"], f = 0.3)
        
        svalue(s_bar) <- paste("Loading map...", sep = "")
        
        map <- suppressWarnings(get_map(location = b_box,
                                        zoom = svalue(d_zoom), 
                                        scale = "auto",
                                        maptype = "hybrid",
                                        filename = "ggmapTemp",
                                        source = "google"))
        
        if("Path" %in% svalue(c_opt) & "Points" %in% svalue(c_opt))
        {
          print(ggmap(map) +
                  geom_point(data = data, colour = "darkorange") +  
                  geom_path(data = data, colour = "firebrick"))
        }
        
        if("Path" %in% svalue(c_opt) & !("Points" %in% svalue(c_opt)))
        {
          print(ggmap(map) +
                  geom_path(data = data, colour = "darkorange"))
        }
        
        if(!("Path" %in% svalue(c_opt)) & "Points" %in% svalue(c_opt))
        {
          print(ggmap(map) +
                  geom_point(data = data, colour = "firebrick"))
        }
        
        if(!("Path" %in% svalue(c_opt)) & !("Points" %in% svalue(c_opt)))
        {
          print(ggmap(map) +
                  geom_point(data = data, colour = "firebrick",alpha = 1/10)+  
                  geom_path(data = data, colour = "firebrick",alpha = 1/10))
        }
        
        svalue(s_bar) <- paste("Map loading complete!", sep = "")
        
      }
    }else{
      
      svalue(s_bar) <- paste("Loading Vessel ", vnum, sep = "")
      vessel <- fn$sqldf("select * from intrp where I_NCEE = `vnum`", dbname = vms_DB$db)
      data <- as.data.frame(cbind(vessel[,"LON"], vessel[,"LAT"]))
      colnames(data) <- c("lon", "lat")
      vestra[] <- unique(vessel["T_NUM"])
      
      svalue(s_bar) <- paste("Loading Vessel Tracks complete!", sep = "")
      
      if(nrow(vessel) != 0)
      {
        
        num_trk <- nrow(vessel)
        
        b_box <- make_bbox(vessel[,"LON"], vessel[,"LAT"], f = 0.3)
        
        svalue(s_bar) <- paste("Loading map...", sep = "")
        
        map <- suppressWarnings(get_map(location = b_box,
                                        zoom = svalue(d_zoom), 
                                        scale = "auto",
                                        maptype = "hybrid",
                                        filename = "ggmapTemp",
                                        source = "google"))
        
        if("Path" %in% svalue(c_opt) & "Points" %in% svalue(c_opt))
        {
          print(ggmap(map) +
                  geom_point(data = data, colour = "darkorange") +  
                  geom_path(data = data, colour = "firebrick"))
        }
        
        if("Path" %in% svalue(c_opt) & !("Points" %in% svalue(c_opt)))
        {
          print(ggmap(map) +
                  geom_path(data = data, colour = "darkorange"))
        }
        
        if(!("Path" %in% svalue(c_opt)) & "Points" %in% svalue(c_opt))
        {
          print(ggmap(map) +
                  geom_point(data = data, colour = "firebrick"))
        }
        
        if(!("Path" %in% svalue(c_opt)) & !("Points" %in% svalue(c_opt)))
        {
          print(ggmap(map) +
                  geom_point(data = data, colour = "firebrick",alpha = 1/10)+  
                  geom_path(data = data, colour = "firebrick",alpha = 1/10))
        }
        
        svalue(s_bar) <- paste("Map loading complete! Select a track!", sep = "")
        
      }
    }  
    enabled(selves) <- TRUE
    enabled(vestra) <- TRUE
    enabled(b_draw) <- TRUE
    enabled(b_save_j) <- TRUE
    
    dev.off()
  })
  
  
  g_lefm <- ggroup(expand = TRUE, horizontal = T, container = g_lef, use.scrollwindow = TRUE)
  selves <- gtable(data.frame("Vessel" = numeric(0)), chosencol = 1, container = g_lefm, expand = T , handler = function(h,...)
  {
    enabled(selves) <- FALSE
    enabled(vestra) <- FALSE
    enabled(b_draw) <- FALSE
    enabled(b_save_j) <- FALSE
    
    svalue(vestra) <- character(0)
    vnum <- svalue(selves)
    
    svalue(s_bar) <- paste("Loading Vessel ", vnum, sep = "")
    vessel <- fn$sqldf("select * from intrp where I_NCEE = `vnum`", dbname = vms_DB$db)
    svalue(s_bar) <- paste("Loading Vessel Tracks complete!", sep = "")
    
    vestra[] <- unique(vessel["T_NUM"])
    
    if(nrow(vessel) != 0)
    {
      
      num_trk <- nrow(vessel)
      
      b_box <- make_bbox(vessel[,"LON"], vessel[,"LAT"], f = 0.3)
      data <- as.data.frame(cbind(vessel[,"LON"], vessel[,"LAT"]))
      colnames(data) <- c("lon", "lat")
      lon_range <- b_box[c("left", "right")]
      lat_range <- b_box[c("bottom", "top")]
      lonlength <- diff(lon_range)
      latlength <- diff(lat_range)
      zoomlon <- round(log2(360 * 2/lonlength))
      zoomlat <- round(log2(180 * 2/latlength))
      zoom <- max(zoomlon, zoomlat)
      
      svalue(d_zoom) <- zoom
      
      svalue(s_bar) <- paste("Loading map...", sep = "")
      
      map <- suppressWarnings(get_map(location = b_box,
                                      zoom = svalue(d_zoom), 
                                      scale = "auto",
                                      maptype = "hybrid",
                                      filename = "ggmapTemp",
                                      source = "google"))
      
      if("Path" %in% svalue(c_opt) & "Points" %in% svalue(c_opt))
      {
        print(ggmap(map) +
                geom_point(data = data, colour = "darkorange") +  
                geom_path(data = data, colour = "firebrick"))
      }
      
      if("Path" %in% svalue(c_opt) & !("Points" %in% svalue(c_opt)))
      {
        print(ggmap(map) +
                geom_path(data = data, colour = "darkorange"))
      }
      
      if(!("Path" %in% svalue(c_opt)) & "Points" %in% svalue(c_opt))
      {
        print(ggmap(map) +
                geom_point(data = data, colour = "firebrick"))
      }
      
      if(!("Path" %in% svalue(c_opt)) & !("Points" %in% svalue(c_opt)))
      {
        print(ggmap(map) +
                geom_point(data = data, colour = "firebrick",alpha = 1/10)+  
                geom_path(data = data, colour = "firebrick",alpha = 1/10))
      }
      
      svalue(s_bar) <- paste("Map loading complete! Select a track!", sep = "")
      
    }
    enabled(selves) <- TRUE
    enabled(vestra) <- TRUE
    enabled(b_draw) <- TRUE
    enabled(b_save_j) <- TRUE
    
  })
  #g_lefb <- ggroup(expand = TRUE, horizontal = T, container = g_lef, use.scrollwindow = TRUE)
  vestra <- gtable(data.frame("Tracks" = numeric(0)), chosencol = 1, container = g_lefm, expand = T , handler = function(h,...)
  {
    enabled(selves) <- FALSE
    enabled(vestra) <- FALSE
    enabled(b_draw) <- FALSE
    enabled(b_save_j) <- FALSE
    
    vnum <- svalue(selves)
    trackn <- svalue(vestra)
    
    s_track  <- fn$sqldf("select * from intrp where I_NCEE = `vnum` and T_NUM = `trackn` order by DATE", dbname = vms_DB$db)
    
    not_intr <- which(s_track[,"P_INT"] == 0)
    data <- as.data.frame(cbind(s_track[,"LON"], s_track[,"LAT"]))
    real_p <- as.data.frame(cbind(s_track[not_intr,"LON"], s_track[not_intr,"LAT"]))
    colnames(data) <- c("lon", "lat")
    colnames(real_p) <- c("lon", "lat")
    
    if(nrow(data) > 2)
    {
      b_box <- make_bbox(s_track[,"LON"], s_track[,"LAT"], f = 0.3)
      lon_range <- b_box[c("left", "right")]
      lat_range <- b_box[c("bottom", "top")]
      lonlength <- diff(lon_range)
      latlength <- diff(lat_range)
      zoomlon <- round(log2(360 * 2/lonlength))
      zoomlat <- round(log2(180 * 2/latlength))
      zoom <- max(zoomlon, zoomlat)
      
      svalue(d_zoom) <- zoom
      
      svalue(s_bar) <- paste("Loading map...", sep = "")
      
      map <- suppressWarnings(get_map(location = b_box,
                                      zoom = svalue(d_zoom), 
                                      scale = "auto",
                                      maptype = "hybrid",
                                      filename = "ggmapTemp",
                                      source = "google"))
      
      if("Path" %in% svalue(c_opt) & "Points" %in% svalue(c_opt))
      {
        print(ggmap(map) +
                geom_point(data = data, colour = "darkorange") +  
                geom_path(data = data, colour = "firebrick"))
      }
      
      if("Path" %in% svalue(c_opt) & !("Points" %in% svalue(c_opt)))
      {
        print(ggmap(map) +
                geom_path(data = data, colour = "darkorange"))
      }
      
      if(!("Path" %in% svalue(c_opt)) & "Points" %in% svalue(c_opt))
      {
        print(ggmap(map) +
                geom_point(data = data, colour = "firebrick"))
      }
      
      if(!("Path" %in% svalue(c_opt)) & !("Points" %in% svalue(c_opt)))
      {
        print(ggmap(map) +
                geom_point(data = data, colour = "firebrick",alpha = 1/10)+  
                geom_path(data = data, colour = "firebrick",alpha = 1/10))
      }
      
      svalue(s_bar) <- paste("Map loading complete!", sep = "")
      
    }
    
    enabled(selves) <- TRUE
    enabled(vestra) <- TRUE
    enabled(b_draw) <- TRUE
    enabled(b_save_j) <- TRUE
  })
  
  g_mid <- ggroup(horizontal = T, container = big_g, expand = T)
  themap <- ggraphics(container = g_mid)
  
  names(selves) <- "Vessel"
  names(vestra) <- "Track"
  
  s_bar <- gstatusbar(text = "", container = adv_view_win)
  
  enabled(d_zoom) <- FALSE
  enabled(b_draw) <- FALSE
  enabled(b_save_j) <- FALSE
  enabled(c_opt) <- FALSE
  enabled(selves) <- FALSE
  enabled(vestra) <- FALSE
  
  visible(adv_view_win) <- TRUE
  
  svalue(s_bar) <- "Loading Map"
  
  map <- suppressWarnings(get_googlemap(center = c(lon = 12.482778, lat = 41.893056),
                                        zoom = 3,
                                        scale = 1,
                                        maptype = "hybrid"))
  
  print(ggmap(map))
  
  svalue(s_bar) <- "Ready! Load a VMS DB file"
  
  if(vms_DB$db != "")
  {
    #     svalue(sel_vms_f) <- strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])]
    svalue(sel_vms_f) <- ifelse(.Platform$OS.type == "windows", strsplit(vms_DB$db, "\\\\")[[1]][length(strsplit(vms_DB$db, "\\\\")[[1]])],strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])])
    svalue(s_bar) <- "Loading DataSet"
    incee <- sqldf("select distinct I_NCEE from intrp order by I_NCEE", dbname = vms_DB$db)
    #colnames(incee) <- "Vessel ID"
    selves[] <- incee
    
    svalue(s_bar) <- "DataSet loaded! Click on a vessel ID to load tracks"
    
    enabled(d_zoom) <- TRUE
    enabled(b_draw) <- TRUE
    enabled(b_save_j) <- TRUE
    enabled(c_opt ) <- TRUE
    enabled(selves) <- TRUE
    enabled(vestra) <- TRUE
  }  
}

#   b_load_vms <- gbutton(text = "Load\nVMS DB", container = g_lef, handler = function(h,...)
#   {
#     
#     vms_DB$db <- gfile(text = "Select VMS DataBase file",
#                        type = "open",
#                        filter = list("VMS DB file" = list(patterns = c("*.vms.sqlite"))))
#     
#     svalue(s_bar) <- "Loading DataSet"
#     incee <- sqldf("select distinct I_NCEE from intrp order by I_NCEE", dbname = vms_DB$db)
#     #colnames(incee) <- "Vessel ID"
#     selves[] <- incee
#     
#     svalue(s_bar) <- "DataSet loaded! Click on a vessel ID to load tracks"
#     
#     enabled(d_zoom) <- TRUE
#     enabled(b_draw) <- TRUE
#     enabled(c_opt ) <- TRUE
#     enabled(selves) <- TRUE
#     enabled(vestra) <- TRUE
#   })
#   