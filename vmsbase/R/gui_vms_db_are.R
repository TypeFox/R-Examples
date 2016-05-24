
#' VMS DB Area Assignment GUI
#' 
#' The \code{gui_vms_db_are} function implements the graphical user interface for the
#'  VMS Assign Area routine.
#' 
#' This function, with a VMS database,
#'  assigns areas to VMS tracks according to median positions.
#' 
#' @param vms_db_name The path of a VMS DataBase
#' 
#' @return This function does not return a value. 
#' 
#' @usage gui_vms_db_are(vms_db_name = "")
#' 
#' @export gui_vms_db_are
#'

gui_vms_db_are <- function(vms_db_name = "")
{
  vms_DB <- vms_DB$new()
  vms_DB$db <- vms_db_name
  themap <- polymap$new()
  themap$path = ""
  
  vms_db_are_win <- gwindow("VMS Data Enhancer Utility - Area", visible = FALSE)
  
  # Assign area
  
  are_g <- gframe(horizontal = FALSE, container = vms_db_are_win)
  are_g2 <- ggroup(horizontal = TRUE, container = are_g)
  are_g3 <- ggroup(horizontal = FALSE, container = are_g)
  
  addSpring(are_g2)
  gimage(system.file("ico/map-compass.png", package="vmsbase"), container = are_g2)
  proglab_area <- glabel("Assign Area" , container = are_g2)
  addSpring(are_g2)
  
  #################
  addSpring(are_g3)
  ##VMS DataBase file
  vms_db_f <- gframe(text = "VMS DB file", horizontal = TRUE, container = are_g3)
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
           if(themap$path != "" & vms_DB$db != "")
           {
             enabled(start_b) <- TRUE
           }
         })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = vms_db_f,
         handler = function(h,...){
           vms_DB$db <- ""
           enabled(start_b) <- FALSE
           svalue(sel_vms_f) <- "Select VMS DB file"
         })
  addSpring(are_g3)
  ################
  
  ##Land file
  cus_map_g <- gframe(text = "Land Shape File", horizontal = TRUE, container = are_g3)
  addSpring(cus_map_g)
  cus_map_lab <- glabel("Select Land Shape File", container = cus_map_g)
  addSpring(cus_map_g)
  gimage(system.file("ico/folder-html.png", package="vmsbase"), container = cus_map_g,
         handler = function(h,...){
           themap$path <- gfile(text = "Select Area ShapePoly map",
                                type = "open",
                                filter = list("shp data" = list(patterns = c("*.shp"))))
           #            svalue(cus_map_lab) <- paste("Area Map: ", strsplit(themap$path, "/")[[1]][length(strsplit(themap$path, "/")[[1]])], sep = "")
           svalue(cus_map_lab) <- paste("Area Map: ", ifelse(.Platform$OS.type == "windows", strsplit(themap$path, "\\\\")[[1]][length(strsplit(themap$path, "\\\\")[[1]])],strsplit(themap$path, "/")[[1]][length(strsplit(themap$path, "/")[[1]])]), sep = "")
           if(themap$path != "" & vms_DB$db != "")
           {
             enabled(start_b) <- TRUE
           }
         })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = cus_map_g,
         handler = function(h,...){
           themap$path <- ""
           svalue(cus_map_lab) <- "Select Land Shape File"
         })
  
  addSpring(are_g)
  
  infolab_area <- glabel("" , container = are_g)
  addSpring(are_g)
  start_b <- gbutton("Start\nArea Annotation", container = are_g, handler = function(h,...)
  {
    enabled(cus_map_g) <- FALSE
    enabled(vms_db_f) <- FALSE
    enabled(start_b) <- FALSE
    
    if(sqldf("select count(*) from track", dbname = vms_DB$db) > 0)
    {
      
      svalue(infolab_area) <- "Upadating...\nVMS DataBase"
      sqldf("drop table if exists p_area", dbname = vms_DB$db)
      sqldf("CREATE TABLE p_area(vess_id INT, t_num INT, AREA INT)", dbname = vms_DB$db)
      
      svalue(infolab_area) <- "Loading...\nArea Map"
      themap$data <- readShapePoly(themap$path)
      
      svalue(infolab_area) <- "Generating...\nArea Box"
      Area_box <- SpatialPolygons2PolySet(themap$data)
      
      incee <- sqldf("select distinct I_NCEE from ping", dbname = vms_DB$db)
      for ( v in 1:nrow(incee) )
      {
        cat("\nAnnotating vessel ", v," of ", nrow(incee), "\n", spe = "")
        svalue(infolab_area) <- paste("Loading...\nVessel: ", v," of ", nrow(incee), spe = "")
        
        vessel <- fn$sqldf("select * from track where I_NCEE = `incee[v,1]` order by DATE", dbname = vms_DB$db)
        
        for(i in unique(vessel[,"T_NUM"]))
        {
          svalue(infolab_area) <- paste("\nVessel: ", v," of ", nrow(incee),
                                        "\nTrack: ", i, " of ", length(unique(vessel[,"T_NUM"])), spe = "")
          
          cat(".", sep = "")
          area <- Assign_Area(vessel[which(vessel[,"T_NUM"] == i),c("LON","LAT")], Area_box)
          
          trk_area <- as.data.frame(cbind(incee[v,1], i, area))
          colnames(trk_area) <- c("I_NCEE", "T_NUM", "AREA")
          
          sqldf("insert into p_area select * from `trk_area`", dbname = vms_DB$db)    
        }
        cat("Completed!\n", sep = "")
        
      }
      cat("\n\n   ---   End Assign Area   ---\n")
      
      gconfirm("VMS DB Area Annotation Completed!",
               title = "Confirm",
               icon = "info",
               parent = vms_db_are_win,
               handler = function(h,...){dispose(vms_db_are_win)})
    }else{
      
      gconfirm("Track data not available\n\nExecute Track Cut first!",
               title = "Error",
               icon = "error",
               parent = vms_db_are_win,
               handler = function(h,...){dispose(vms_db_are_win)})
      
    }
    enabled(cus_map_g) <- TRUE
    enabled(vms_db_f) <- TRUE
    enabled(start_b) <- TRUE
  })
  enabled(start_b) <- FALSE
  
  if(vms_DB$db != "")
  {
    #     svalue(sel_vms_f) <- strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])]
    svalue(sel_vms_f) <- ifelse(.Platform$OS.type == "windows", strsplit(vms_DB$db, "\\\\")[[1]][length(strsplit(vms_DB$db, "\\\\")[[1]])],strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])])
  } 
  if(themap$path != "")
  {
    #   svalue(cus_map_lab) <- paste("Area Map: ", strsplit(themap$path, "/")[[1]][length(strsplit(themap$path, "/")[[1]])], sep = "")
    svalue(cus_map_lab) <- paste("Area Map: ", ifelse(.Platform$OS.type == "windows", strsplit(themap$path, "\\\\")[[1]][length(strsplit(themap$path, "\\\\")[[1]])],strsplit(themap$path, "/")[[1]][length(strsplit(themap$path, "/")[[1]])]), sep = "")
  }
  if(themap$path != "" & vms_DB$db != "")
  {
    enabled(start_b) <- TRUE
  }
  visible(vms_db_are_win) <- TRUE
  
}