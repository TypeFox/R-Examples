
#' VMS DB View Raw Ping Data GUI
#'  
#' The \code{gui_vms_view_ping} function implements the graphical user interface for the
#'  VMS DB routine to view raw ping data.
#' 
#' This function, with a VMS database,
#'  plots vessel data.
#'   
#' @param vms_db_name The path of a VMS DataBase
#' @param bathy_file_name The path of a Bathymetry file
#'   
#' @return This function does not return a value. 
#' 
#' @usage gui_vms_view_ping(vms_db_name = "", bathy_file_name = "")
#' 
#' @export gui_vms_view_ping
#'
#'@seealso \code{\link{gui_vms_save_bat}} \code{\link{gui_vms_view_track}} \code{\link{gui_vms_view_intrp}}

gui_vms_view_ping <- function (vms_db_name = "", bathy_file_name = "")
{
  vms_DB <- vms_DB$new()
  vms_DB$db <- vms_db_name
  bathy <- bathymetry$new()
  bathy$path <- bathy_file_name
  
  #Avvio Interfaccia
  
  ping_view_win <- gwindow("Ping View Device", visible = FALSE)
  
  big_g <- ggroup(horizontal = TRUE, container = ping_view_win)
  left_g <- ggroup(horizontal = FALSE, container = big_g)
  chk_g3 <- ggroup(horizontal = TRUE, container = left_g)
  expo_gr <- ggroup(horizontal = TRUE, container = left_g)
  
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
           
           incee <- sqldf("select distinct I_NCEE from ping order by I_NCEE", dbname = vms_DB$db)
           selves[] <- incee
           
           enabled(selves) <- TRUE
           enabled(expo_gr) <- TRUE
           enabled(bbox_exp) <- TRUE
         })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = vms_db_f,
         handler = function(h,...){
           vms_DB$db <- ""
           enabled(selves) <- FALSE
           enabled(expo_gr) <- FALSE
           enabled(bbox_exp) <- FALSE
           svalue(sel_vms_f) <- "Select VMS DB file"
         })
  addSpring(chk_g3)
  ################
  addSpring(expo_gr)
  save_jpeg <- gbutton(text = "save jpeg", container = expo_gr, handler = function(h,...)
  {
    enabled(ping_view_win) <- FALSE
    
    xle <- par()$pin[1]
    yle <- par()$pin[2]
    
    j_fil_na <- gfile(text = "Jpeg file path and name",
                      initialfilename = "*.jpeg",
                      type = "save")
    if(length(unlist(strsplit(j_fil_na, "[.]"))) == 1){j_fil_na <- paste(j_fil_na, ".jpeg", sep = "")}
    
    
    jpeg(filename = j_fil_na,
         width = round(1000*xle), height = round(1000*yle),
         quality = 100, bg = "aliceblue", pointsize = 80)
    par(lwd = 10)
    vnum <- svalue(selves)
    vessel <- fn$sqldf("select * from ping where I_NCEE = `vnum`", dbname = vms_DB$db)
    if(nrow(vessel) != 0)
    {
      xrange <- c(svalue(mi_lo), svalue(ma_lo))
      yrange <- c(svalue(mi_la), svalue(ma_la))
      pingview(vessel, bathy, xrange, yrange)
    }
    dev.off()
    enabled(ping_view_win) <- TRUE
  })
  addSpring(expo_gr)
  export_csv <- gbutton(text = "export csv", container = expo_gr, handler = function(h,...)
  {
    enabled(ping_view_win) <- FALSE
    vnum <- svalue(selves)
    vessel <- fn$sqldf("select * from ping where I_NCEE = `vnum` order by DATE", dbname = vms_DB$db)
    if(nrow(vessel) != 0)
    {
      csv_fil_na <- gfile(text = "Save vessel route as CSV file", type = "save", initialfilename = "*.csv", 
                          filter = list("All files" = list(patterns = c("*")), "CSV files" =
                                          list(patterns = "*.csv")))
      
      if(length(unlist(strsplit(csv_fil_na, "[.]"))) == 1){csv_fil_na <- paste(csv_fil_na, ".csv", sep = "")}
      
      write.table(vessel,
                  file = csv_fil_na,
                  append = FALSE,
                  sep = ";",
                  dec = ".",
                  row.names = FALSE,
                  col.names = TRUE)
    }
    enabled(ping_view_win) <- TRUE
  })
  addSpring(expo_gr)
  
  
  bbox_exp <- gexpandgroup(text = "Custom B-Box", container = left_g, horizontal = TRUE)
  addSpring(bbox_exp)
  bbox_lay <- glayout(container = bbox_exp)
  bbox_lay[1,2] <- ggroup(horizontal = FALSE)
  glabel("Max Lat", container = bbox_lay[1,2])
  ma_la <- gspinbutton(from = -90, to = 90, by = 0.5, value = 0, container = bbox_lay[1,2])
  bbox_lay[2,1] <- ggroup(horizontal = FALSE)
  glabel("Min Lon", container = bbox_lay[2,1])
  mi_lo <- gspinbutton(from = -180, to = 180, by = 0.5, value = 0, container = bbox_lay[2,1])
  bbox_lay[2,3] <- ggroup(horizontal = FALSE)
  glabel("Max Lon", container = bbox_lay[2,3])
  ma_lo <- gspinbutton(from = -180, to = 180, by = 0.5, value = 0, container = bbox_lay[2,3])
  bbox_lay[3,2] <- ggroup(horizontal = FALSE)
  glabel("Min Lat", container = bbox_lay[3,2])
  mi_la <- gspinbutton(from = -90, to = 90, by = 0.5, value = 0, container = bbox_lay[3,2])
  addSpring(bbox_exp)
  re_plot <- gbutton(text = "Custom Plot", container = bbox_exp, handler = function(h,...){
    
    vnum <- svalue(selves)
    vessel <- fn$sqldf("select * from ping where I_NCEE = `vnum`", dbname = vms_DB$db)
    
    xrange <- c(svalue(mi_lo), svalue(ma_lo))
    yrange <- c(svalue(mi_la), svalue(ma_la))
    
    pingview(vessel, bathy, xrange, yrange)
  })
  addSpring(bbox_exp)
  enabled(bbox_exp) <- FALSE
  
  chk_g4 <- ggroup(horizontal = TRUE, expand = TRUE, container = left_g)
  selves <- gtable(items = data.frame("Vessel" = numeric(0)), chosencol = 1, container = chk_g4, expand = TRUE, handler = function(h,...)
  {
    enabled(ping_view_win) <- FALSE
    enabled(expo_gr) <- FALSE
    
    vnum <- svalue(selves)
    vessel <- fn$sqldf("select * from ping where I_NCEE = `vnum`", dbname = vms_DB$db)
    if(nrow(vessel) != 0)
    {
      span <- 0.25
      xrange <- extendrange(x = vessel["LON"], f = span)
      yrange <- extendrange(x = vessel["LAT"], f = span)
      svalue(ma_la) <- yrange[2]
      svalue(mi_lo) <- xrange[1]
      svalue(ma_lo) <- xrange[2]
      svalue(mi_la) <- yrange[1]
      pingview(vessel, bathy, xrange, yrange)
    }
    enabled(ping_view_win) <- TRUE
    enabled(expo_gr) <- TRUE
  })
  enabled(selves) <- FALSE
  enabled(expo_gr) <- FALSE
  ###################
  cus_dep_g <- gframe(text = "Bathymetry File", horizontal = TRUE, container = left_g)
  addSpring(cus_dep_g)
  cus_dep_lab <- glabel("Select Bathymetry File", container = cus_dep_g)
  addSpring(cus_dep_g)
  gimage(system.file("ico/folder-download.png", package="vmsbase"), container = cus_dep_g,
         handler = function(h,...){
           bathy$path <- gfile(text = "Select Bathymetry File",
                               type = "open",
                               filter = list("bathy data" = list(patterns = c("*sqlitebathy.rData"))))
           #            svalue(cus_dep_lab) <- paste("File: ", strsplit(bathy$path, "/")[[1]][length(strsplit(bathy$path, "/")[[1]])], sep = "")
           
           svalue(cus_dep_lab) <- paste("File: ", ifelse(.Platform$OS.type == "windows", strsplit(bathy$path, "\\\\")[[1]][length(strsplit(bathy$path, "\\\\")[[1]])],strsplit(bathy$path, "/")[[1]][length(strsplit(bathy$path, "/")[[1]])]), sep = "")
           
           
           bathy$data <- readRDS(bathy$path)
         })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = cus_dep_g,
         handler = function(h,...){
           bathy$path <- ""
           svalue(cus_dep_lab) <- "Select Bathymetry File"
         })
  ###################
  right_g <- gframe(horizontal = F, container = big_g, expand = T) 
  theplot <- ggraphics(container = right_g, width = 600, height = 450)
  visible(ping_view_win) <- TRUE
  maps::map("world",col="black", bg = "lightsteelblue1", mar=c(6,6,0,0), fill = TRUE, interior = FALSE)
  maps::map.axes()
  title(main = "Ping Viewer", line = 0.3)
  title(xlab = "Lon", ylab = "Lat", line = 2)
  if(vms_DB$db != "")
  {
    enabled(ping_view_win) <- FALSE
    
    #     svalue(sel_vms_f) <- strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])]
    svalue(sel_vms_f) <- ifelse(.Platform$OS.type == "windows", strsplit(vms_DB$db, "\\\\")[[1]][length(strsplit(vms_DB$db, "\\\\")[[1]])],strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])])
    
    incee <- sqldf("select distinct I_NCEE from ping order by I_NCEE", dbname = vms_DB$db)
    selves[] <- incee
    enabled(ping_view_win) <- TRUE
#     enabled(selves) <- TRUE
#     enabled(bbox_exp) <- TRUE
  }
  if(bathy$path != "")
  {
    #   svalue(cus_dep_lab) <- paste("File: ", strsplit(bathy$path, "/")[[1]][length(strsplit(bathy$path, "/")[[1]])], sep = "")
    svalue(cus_dep_lab) <- paste("File: ", ifelse(.Platform$OS.type == "windows", strsplit(bathy$path, "\\\\")[[1]][length(strsplit(bathy$path, "\\\\")[[1]])],strsplit(bathy$path, "/")[[1]][length(strsplit(bathy$path, "/")[[1]])]), sep = "")
    
    bathy$data <- readRDS(bathy$path)
  }
}
