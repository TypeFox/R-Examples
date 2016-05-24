
#' VMS DB Save Bathymetry GUI
#'  
#' The \code{gui_vms_save_bat} function implements the graphical user interface for the
#'  VMS Save Bathymetry routine.
#' 
#' This function, with a VMS database,
#'  downloads an XYZ bathymetry file from NOAA servers.
#'   
#' @param vms_db_name The path of a VMS DataBase
#' 
#' @return This function does not return a value. 
#' 
#' @usage gui_vms_save_bat(vms_db_name = "")
#' 
#' @export gui_vms_save_bat
#'
#'@seealso \code{\link{gui_vms_view_ping}} \code{\link{gui_vms_view_track}} \code{\link{gui_vms_view_intrp}}

gui_vms_save_bat <- function(vms_db_name = "")
{
  vms_DB <- vms_DB$new()
  vms_DB$db <- vms_db_name
  
  vms_db_bat_win <- gwindow("Bathymetry and Isobaths Download Tool", visible = FALSE)
  
  bat_g <- gframe(horizontal = FALSE, container = vms_db_bat_win)
  bat_g2 <- ggroup(horizontal = TRUE, container = bat_g)
  addSpring(bat_g2)
  gimage(system.file("ico/download-3.png", package = "vmsbase"), container = bat_g2)
  proglab_bat <- glabel("Get Isobaths" , container = bat_g2)
  addSpring(bat_g2)
  
  
  bat_g3 <- ggroup(horizontal = TRUE, container = bat_g)
  
  #################
  addSpring(bat_g3)
  vms_db_f <- gframe(text = "VMS DB file", horizontal = TRUE, container = bat_g3)
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
           enabled(start_b) <- TRUE
         })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = vms_db_f,
         handler = function(h,...){
           vms_DB$db <- ""
           enabled(start_b) <- FALSE
           svalue(sel_vms_f) <- "Select VMS DB file"
         })
  addSpring(bat_g3)
  ################
  
  res_g <- ggroup(horizontal = TRUE, container = bat_g3)
  addSpring(res_g)
  glabel("Set Resolution", container = res_g)
  use_res <- gdroplist(1:4,
                       selected = 2, container = res_g) 
  addSpring(res_g)
  
  
  addSpring(bat_g)
  
  
  theplot <- ggraphics(container = bat_g, width = 400, height = 400)
  
  addSpring(bat_g)
  start_b <- gbutton("Start\nBathymetry Download", container = bat_g, handler = function(h,...)
  {
    enabled(vms_db_bat_win) <- FALSE
    
#     enabled(res_g) <- FALSE
#     enabled(vms_db_f) <- FALSE
#     enabled(start_b) <- FALSE
    
    pings <- sqldf("select * from ping", dbname = vms_DB$db)
    
    xmax <- as.numeric(sqldf("select max(LON) from track", dbname = vms_DB$db))
    xmin <- as.numeric(sqldf("select min(LON) from track", dbname = vms_DB$db))
    ymax <- as.numeric(sqldf("select max(LAT) from track", dbname = vms_DB$db))
    ymin <- as.numeric(sqldf("select min(LAT) from track", dbname = vms_DB$db))
    
    xrange = extendrange(c(xmax, xmin), f = 0.05)
    yrange = extendrange(c(ymax, ymin), f = 0.05)
    
#     bat_blo <- getNOAA.bathy(xrange[1],
#                              xrange[2],
#                              yrange[1],
#                              yrange[2], resolution = svalue(use_res))
    
    bat_blo <- getNOAA.bathy(xmin-0.5,
                             xmax+0.5,
                             ymin-0.5,
                             ymax+0.5, resolution = svalue(use_res))
    
    isob <- c(-20, -50, -100,-200)
    lon <- unique(as.numeric(rownames(bat_blo)))
    lat <- unique(as.numeric(colnames(bat_blo)))

    icol <- rgb(0,0, seq(255,100, len = length(isob)), maxColorValue = 255)

    par(mar = c(1,1,1,1))
    blues<-c("lightsteelblue4","lightsteelblue3","lightsteelblue2","lightsteelblue1")
    plot(bat_blo,image=TRUE,land=TRUE,lwd=0.1,bpal=list(c(0,max(bat_blo),"grey"),c(min(bat_blo),0,blues)))
    plot(bat_blo,deep=0,shallow=0,step=0,lwd=0.4,add=TRUE)
    points(pings[,"LON"], pings[,"LAT"], col = "firebrick", cex = 0.1)
#     plot(pings[,"LON"], pings[,"LAT"],
#          xlab = "Longitude",
#          ylab = "Latitude",
#          mar = c(1,1,0,0))
#     contour(lon, lat, bat_blo, add=TRUE,
#             lwd=0.3,
#             col=icol, 
#             levels = c(-20, -50, -100,-200),
#             drawlabel=FALSE,
#             xlab="a")

    gconfirm("Bathymetry downloas Complete!\n\nSave current Bathymetry and Isobaths?",
             title="Confirm",
             icon = "question",
             parent = vms_db_bat_win,
             handler = function(h,...){
               saveRDS(bat_blo, file = paste(vms_DB$db, "bathy.rData", sep = ""))
               dispose(vms_db_bat_win)
             })
  })
  enabled(start_b) <- FALSE
  
  if(vms_DB$db != "")
  {
#     svalue(sel_vms_f) <- strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])]
    svalue(sel_vms_f) <- ifelse(.Platform$OS.type == "windows", strsplit(vms_DB$db, "\\\\")[[1]][length(strsplit(vms_DB$db, "\\\\")[[1]])],strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])])
    enabled(start_b) <- TRUE
  } 
  visible(vms_db_bat_win) <- TRUE
}