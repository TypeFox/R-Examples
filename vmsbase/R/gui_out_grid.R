
#' VMS Effort Gridding GUI
#' 
#' The \code{gui_out_grid} function implements the graphical user interface
#'  for the VMS Effort Gridding
#' 
#' This function, with a VMS DB and a Grid Sea Area Map shape file, computes the total
#'  fishing effort (in hours) over each cell of the submitted grid, relative to the
#'  selected metier
#'
#' @param vms_db_name The path of a VMS DataBase
#' 
#' @return This function does not return a value.
#'  The result count will be plotted on the submitted grid. The user can both save the 
#'  result count vector as an r object (necessary for \code{\link{gui_dcf_ind}}),
#'   or the annotated grid shape file.
#'  
#' @usage gui_out_grid(vms_db_name = "")
#' 
#' @export gui_out_grid
#'
#'@seealso \code{\link{gui_dcf_ind}}

gui_out_grid <- function(vms_db_name = "")
{
  vms_DB <- vms_DB$new()
  vms_DB$db <- vms_db_name
  count <- 0
  themap <- polymap$new()
  themap$path <- ""
  
  vms_out_grid_win <- gwindow("VMS Effort Gridding", visible = FALSE)
  
  # VMS DB check for pings in harbour, coherence and on land
  big_g <- ggroup(horizontal = TRUE, container =  vms_out_grid_win)
  
  gri_g <- gframe(horizontal = FALSE, container = big_g)
  rigth_g <- ggroup(horizontal = FALSE, container = big_g)
  gri_g2 <- ggroup(horizontal = TRUE, container = gri_g)
  gri_g3 <- ggroup(horizontal = FALSE, container = gri_g)
  gri_g3f <- ggroup(horizontal = TRUE, container = gri_g3)
  gri_g3f2 <- ggroup(horizontal = TRUE, container = gri_g3)
  
  
  addSpring(gri_g2)
  gimage(system.file("ico/view-calendar-timeline.png", package = "vmsbase"), container = gri_g2)
  proglab <- glabel("Gridding & Mapping" , container = gri_g2)
  addSpring(gri_g2)
  
  #################
  addSpring(gri_g3)
  addSpring(gri_g3f)
  vms_db_f <- gframe(text = "VMS DB file", horizontal = TRUE, container = gri_g3f)
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
           
           if(vms_DB$db != "")
           {
             met_sel_d[] <- sqldf("select distinct(met_des) from vms_lb", dbname = vms_DB$db)[,1]
             svalue(met_sel_d, index = TRUE) <- 1
             enabled(gri_g3f3) <- TRUE
             enabled(gri_g3f4) <- TRUE
             if(themap$path != "")
             {
               enabled(start_b) <- TRUE
             }
           }
           
         })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = vms_db_f,
         handler = function(h,...){
           vms_DB$db <- ""
           enabled(start_b) <- FALSE
           enabled(set_count_g) <- FALSE
           enabled(save_j_b) <- FALSE
           enabled(gri_g3f3) <- FALSE
           enabled(gri_g3f4) <- FALSE
           enabled(save_vesh_g) <- FALSE
           svalue(sel_vms_f) <- "Select VMS DB file"
         })
  addSpring(gri_g3f)
  addSpring(gri_g3)
  ################
  addSpring(gri_g3f2)
  ##Land file
  cus_map_g <- gframe(text = "Grid Shape File", horizontal = TRUE, container = gri_g3f2)
  addSpring(cus_map_g)
  cus_map_lab <- glabel("Select Grid Shape File", container = cus_map_g)
  addSpring(cus_map_g)
  gimage(system.file("ico/folder-html.png", package="vmsbase"), container = cus_map_g,
         handler = function(h,...){
           
           themap$path <- gfile(text = "Select ShapePoly map",
                                type = "open",
                                filter = list("shp data" = list(patterns = c("*.shp"))))
           svalue(cus_map_lab) <- paste("Grid map: ", ifelse(.Platform$OS.type == "windows", strsplit(themap$path, "\\\\")[[1]][length(strsplit(themap$path, "\\\\")[[1]])],strsplit(themap$path, "/")[[1]][length(strsplit(themap$path, "/")[[1]])]), sep = "")
           
           
           if(themap$path != "" & vms_DB$db != "")
           {
             met_sel_d[] <- sqldf("select distinct(met_des) from vms_lb", dbname = vms_DB$db)[,1]
             svalue(met_sel_d, index = TRUE) <- 1
             enabled(gri_g3f3) <- TRUE
             enabled(gri_g3f4) <- TRUE
             enabled(start_b) <- TRUE
           }
         })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = cus_map_g,
         handler = function(h,...){
           themap$path <- ""
           svalue(cus_map_lab) <- "Select ShapePoly map"
           
           enabled(start_b) <- FALSE
           enabled(set_count_g) <- FALSE
           enabled(save_j_b) <- FALSE
           enabled(gri_g3f3) <- FALSE
           enabled(gri_g3f4) <- FALSE
           enabled(save_vesh_g) <- FALSE
         })
  addSpring(gri_g3f2)
  
  gri_g3f3 <- ggroup(horizontal = TRUE, container = gri_g3)
  addSpring(gri_g3f3)
  met_sel_f <- gframe(text = "Metier Selection", horizontal=TRUE, container = gri_g3f3) 
  met_sel_d <- gdroplist(c("Metier"), selected = 1, container = met_sel_f) 
  addSpring(gri_g3f3)
  glabel("Interpolation\nfrequency" , container = gri_g3f3)
  int_fre_cho <- gspinbutton(from = 1, to = 60, by = 1, value = 10, horizontal = TRUE, container = gri_g3f3)
  glabel("minutes" , container = gri_g3f3)
  addSpring(gri_g3f3)
  enabled(gri_g3f3) <- FALSE
  
  gri_g3f4 <- ggroup(horizontal = TRUE, container = gri_g3)
  addSpring(gri_g3f4)
  dat_sel_f <- gframe(text = "Metier Data Source", horizontal=TRUE, container = gri_g3f4) 
  dat_sel_d <- gdroplist(c("VMS-LB Match", "NN Prediction"), selected = 1, container = dat_sel_f, handler = function(h,...)
  {
    met_sel_d[] <- sqldf(paste("select distinct met_des from ",ifelse(svalue(dat_sel_d) == "VMS-LB Match", "vms_lb","nn_clas"),sep = ""), dbname = vms_DB$db)
  }
  ) 
  addSpring(gri_g3f4)
  enabled(gri_g3f4) <- FALSE
  
  addSpring(gri_g)
  infolab <- glabel("" , container = gri_g)
  addSpring(gri_g)
  
  start_b <- gbutton("Start gridding", container = gri_g, handler = function(h,...)
  {
    enabled(cus_map_g) <- FALSE
    enabled(vms_db_f) <- FALSE
    enabled(start_b) <- FALSE
    enabled(set_count_g) <- FALSE
    enabled(save_j_b) <- FALSE
    enabled(gri_g3f3) <- FALSE
    enabled(gri_g3f4) <- FALSE
    enabled(save_vesh_g) <- FALSE
    
    svalue(infolab) <- "Gridding Started"
    count <<- 0
    cat("\n\n   ---   VMS DB Gridding Started   ---\n", sep = "")
    
    cat("\n   -     Loading Data From DB...", sep = "")
    sel_met <- svalue(met_sel_d)
    
    re_pi <- data.frame()
    if(svalue(dat_sel_d) == "VMS-LB Match")
    {
      f_poi <- sqldf("select i_id from p_fish where FISH = 1", dbname = vms_DB$db)
      
      if(nrow(f_poi) > 0)
      {
        p_met_mat <- fn$sqldf("select intrp.ROWID, LON, LAT, met_des from intrp, vms_lb where vms_lb.vessel = intrp.I_NCEE and intrp.T_NUM = vms_lb.track and vms_lb.met_des = '`sel_met`'", dbname = vms_DB$db)
        re_pi <- p_met_mat[which(p_met_mat[,1] %in% f_poi[,1]),2:3]
        cat("Completed!     -", sep = "")
      }else{
        cat("Error! No Fishing points!")
      }
    }else{
      nn_tab <- as.numeric(sqldf("SELECT count(*) FROM sqlite_master WHERE type='table' AND name='nn_clas'", dbname = vms_DB$db))
      if(nn_tab == 1)
      {
        f_poi <- sqldf("select i_id from p_fish_nn where FISH = 1", dbname = vms_DB$db)
        if(nrow(f_poi) > 0)
        {
          p_met_pre <- fn$sqldf("select intrp.ROWID, LON, LAT, met_des from intrp, nn_clas where nn_clas.I_NCEE = intrp.I_NCEE and intrp.T_NUM = nn_clas.T_NUM and nn_clas.met_des = '`sel_met`'", dbname = vms_DB$db)
          re_pi <- p_met_pre[which(p_met_pre[,1] %in% f_poi[,1]),2:3]
          
          cat("Completed!     -", sep = "")
        }else{
          cat("Error! No Fishing points!")
        }
      }else{
        cat("\n\n -     ERROR! No Neural Network Prediction Available     -\n", sep = "")
      }
    }
    
    if(nrow(re_pi) == 0)
    {
      cat("\n   -     ERROR! No Fishing Points found for metier ", sel_met, sep = "")
    }else{
      svalue(infolab) <- "Loading...\nGrid Map"
      cat("\n   -     Loading Grid Map...", sep = "")
      
      themap$data <- readShapePoly(themap$path)
      ps_themap <- SpatialPolygons2PolySet(themap$data)
      
      cat(" Completed!\n   -     Counting...", sep = "")
      count <<- CountMap(re_pi, ps_themap)
      count <<- (count*svalue(int_fre_cho))/60
      cat(" Completed!\n\n", sep = "")
      
      ncolori <- 10
      s.cols <- plotrix::color.scale(1:ncolori, extremes=c("white","firebrick4"))
      
      log_count <- log(count+1)
      s.int <- seq(floor(min(count)), ceiling(max(count)), length=(length(s.cols)))
      s.vals <- findInterval(count, s.int)
      svalue(infolab) <- "Plotting...\nGrid Map"
      map("worldHires", fill=T, col="springgreen4", bg = "deepskyblue4",
          mar = c(0,0,0,0),
          xlim=c(themap$data@bbox[1,1],themap$data@bbox[1,2]*1.11),
          ylim=c(themap$data@bbox[2,1]*0.98,themap$data@bbox[2,2])
      )
      
      if(sum(count) != 0)
      {
        plot(themap$data, lty = "blank", 
             col=s.cols[s.vals], add = T)
        
        map("worldHires", fill=T, col="springgreen4",
            mar = c(0,0,0,0),
            xlim=c(themap$data@bbox[1,1],themap$data@bbox[1,2]*1.11),
            ylim=c(themap$data@bbox[2,1]*0.98,themap$data@bbox[2,2]), add = TRUE)
        
        legend(x = "bottomright",
               legend = round(s.int,1),
               cex = 0.7,
               title = "Fishing Times",
               col = s.cols,
               bg = "aliceblue",
               lty = 1, lwd = 1, pch = 15)
      }else{
        cat("\n -     No Fishing Points found for metier ", sel_met, " in the submitted area     -\n", sep = "")
      }
      map.axes()
      map.scale(cex = 0.75, ratio = FALSE)
      svalue(infolab) <- "Gridding\nCompleted!"
    }
    enabled(cus_map_g) <- TRUE
    enabled(vms_db_f) <- TRUE
    enabled(start_b) <- TRUE
    enabled(set_count_g) <- TRUE
    enabled(save_j_b) <- TRUE
    enabled(gri_g3f3) <- TRUE
    enabled(gri_g3f4) <- TRUE
    enabled(save_vesh_g) <- TRUE
  })
  enabled(start_b) <- FALSE
  
  addSpring(gri_g)
  
  theplot <- ggraphics(container = rigth_g, width = 600, height = 400)
  
  addSpring(gri_g)
  
  set_count_g <- ggroup(horizontal = TRUE, container = gri_g)
  addSpring(set_count_g)
  set_count_yn <- gradio(c("Sum", "Log10(Sum+1)"),
                         selected = 1,
                         container = set_count_g,
                         horizontal = TRUE,
                         handler = function(h,...){
                           
                           enabled(cus_map_g) <- FALSE
                           enabled(vms_db_f) <- FALSE
                           enabled(start_b) <- FALSE
                           enabled(set_count_g) <- FALSE
                           enabled(save_j_b) <- FALSE
                           enabled(gri_g3f3) <- FALSE
                           enabled(gri_g3f4) <- FALSE
                           enabled(save_vesh_g) <- FALSE
                           
                           #####################
                           
                           ncolori <- 10
                           s.cols <- color.scale(1:ncolori, extremes=c("white","firebrick4"))
                           
                           if(svalue(set_count_yn) == "Sum")
                           {
                             s.int <- seq(floor(min(count)), ceiling(max(count)), length=(length(s.cols)))
                             s.vals <- findInterval(count, s.int)
#                              lege <- round(s.int,0)
                           }else{
                             log_count <- log(count+1)
                             s.int <- seq(floor(min(log_count)), ceiling(max(log_count)), length=(length(s.cols)))
                             s.vals <- findInterval(log_count, s.int)
#                              lege <- round(s.int,1)
                           }
#                            s.vals <- findInterval(log_count, s.int)
                           lege <- round(s.int,1)
                           map("worldHires", fill=T, col="springgreen4", bg = "deepskyblue4",
                               mar = c(0,0,0,0),
                               xlim=c(themap$data@bbox[1,1],themap$data@bbox[1,2]*1.11),
                               ylim=c(themap$data@bbox[2,1]*0.98,themap$data@bbox[2,2])
                           )
                           
                           plot(themap$data, lty= "blank", 
                                col=s.cols[s.vals], add = T)
                           
                           map("worldHires", fill=T, col="springgreen4",
                               mar = c(0,0,0,0),
                               xlim=c(themap$data@bbox[1,1],themap$data@bbox[1,2]*1.11),
                               ylim=c(themap$data@bbox[2,1]*0.98,themap$data@bbox[2,2]), add = TRUE
                           )
                           
                           map.axes()
                           map.scale(cex = 0.75, ratio = FALSE)
                           
                           legend(x = "bottomright",
                                  legend = round(s.int,1),
                                  cex = 0.7,
                                  title =  ifelse(svalue(set_count_yn) == "Sum", "Fishing Times","Log Fishing Times"),
                                  col = s.cols,
                                  bg = "aliceblue",
                                  lty = 1, lwd = 1, pch = 15
                           )
                           
                           #####################
                           
                           enabled(cus_map_g) <- TRUE
                           enabled(vms_db_f) <- TRUE
                           enabled(start_b) <- TRUE
                           enabled(set_count_g) <- TRUE
                           enabled(save_j_b) <- TRUE
                           enabled(gri_g3f3) <- TRUE
                           enabled(gri_g3f4) <- TRUE
                           enabled(save_vesh_g) <- TRUE
                           
                         })
  addSpring(set_count_g)
  enabled(set_count_g) <- FALSE
  
  addSpring(gri_g)
  
  save_vesh_g <- ggroup(horizontal = TRUE, container = gri_g)
  addSpring(save_vesh_g)
  save_j_b <- gbutton("Save Jpeg", container = save_vesh_g, handler = function(h,...)
  {
    enabled(cus_map_g) <- FALSE
    enabled(vms_db_f) <- FALSE
    enabled(start_b) <- FALSE
    enabled(set_count_g) <- FALSE
    enabled(save_j_b) <- FALSE
    enabled(gri_g3f3) <- FALSE
    enabled(gri_g3f4) <- FALSE
    enabled(save_vesh_g) <- FALSE
    
    jpeg(filename = gfile(text = "Save Jpeg File",
                          type = "save"),
         width = 800, height = 800,
         quality = 100)
    
    ncolori <- 10
    s.cols <- color.scale(1:ncolori, extremes=c("white","firebrick4"))
    
    if(svalue(set_count_yn) == "Sum")
    {
      s.int <- seq(floor(min(count)), ceiling(max(count)), length=(length(s.cols)))
                                   s.vals <- findInterval(count, s.int)
      #                              lege <- round(s.int,0)
    }else{
      log_count <- log(count+1)
      s.int <- seq(floor(min(log_count)), ceiling(max(log_count)), length=(length(s.cols)))
                                   s.vals <- findInterval(log_count, s.int)
      #                              lege <- round(s.int,1)
    }
#     s.vals <- findInterval(log_count, s.int)
    lege <- round(s.int,1)
    map("worldHires", fill=T, col="springgreen4", bg = "deepskyblue4",
        mar = c(0,0,0,0),
        xlim=c(themap$data@bbox[1,1],themap$data@bbox[1,2]*1.11),
        ylim=c(themap$data@bbox[2,1]*0.98,themap$data@bbox[2,2])
    )
    
    plot(themap$data, lty= "blank", 
         col=s.cols[s.vals], add = T)
    
    map("worldHires", fill=T, col="springgreen4",
        mar = c(0,0,0,0),
        xlim=c(themap$data@bbox[1,1],themap$data@bbox[1,2]*1.11),
        ylim=c(themap$data@bbox[2,1]*0.98,themap$data@bbox[2,2]), add = TRUE
    )
    
    map.axes()
    map.scale(cex = 0.75, ratio = FALSE)
    
    legend(x = "bottomright",
           legend = round(s.int,1),
           cex = 0.7,
           title =  ifelse(svalue(set_count_yn) == "Sum", "Fishing Times","Log Fishing Times"),
           col = s.cols,
           bg = "aliceblue",
           lty = 1, lwd = 1, pch = 15
    )
    
    dev.off()
    
    enabled(cus_map_g) <- TRUE
    enabled(vms_db_f) <- TRUE
    enabled(start_b) <- TRUE
    enabled(set_count_g) <- TRUE
    enabled(save_j_b) <- TRUE
    enabled(gri_g3f3) <- TRUE
    enabled(gri_g3f4) <- TRUE
    enabled(save_vesh_g) <- TRUE
  })
  enabled(save_j_b) <- FALSE
  addSpring(save_vesh_g)
  save_vec <- gbutton("    Save\nCount Vector", container = save_vesh_g, handler = function(h,...)
  {
    saveRDS(count, paste(gfile(text = "Select dir and Vector Name", type = "save", filter = list("All files" = list(patterns = c("*")), "R grid Vector" = list(patterns = c(".grivec.rData"))) ), ".grivec.rData", sep = ""))
  })
  addSpring(save_vesh_g)
  save_shp <- gbutton("    Save\nShape File", container = save_vesh_g, handler = function(h,...)
  {
    cho_des <- gfile(text = "Select Shapefile Destination Directory", type = "selectdir")
    Join2shp(shpfile = themap$path ,datavector = count, dirdest = cho_des)
    gconfirm("VMS Grid Shapefile Saving Completed!", title = "Confirm", icon = "info")
  })
  enabled(save_vesh_g) <- FALSE
  addSpring(save_vesh_g)
  if(vms_DB$db != "")
  {
    #     svalue(sel_vms_f) <- strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])]
    svalue(sel_vms_f) <- ifelse(.Platform$OS.type == "windows", strsplit(vms_DB$db, "\\\\")[[1]][length(strsplit(vms_DB$db, "\\\\")[[1]])],strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])])
    
    met_sel_d[] <- sqldf("select distinct(met_des) from vms_lb", dbname = vms_DB$db)[,1]
    svalue(met_sel_d, index = TRUE) <- 1
    enabled(gri_g3f3) <- TRUE
    enabled(gri_g3f4) <- TRUE
    if(themap$path != "")
    {
      enabled(start_b) <- TRUE
    }
  }  
  visible(vms_out_grid_win) <- TRUE
}

