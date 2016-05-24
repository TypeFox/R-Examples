

#' VMS DB Clean GUI
#'  
#' The \code{gui_vms_db_clean} function implements the graphical user interface for the
#'  VMS DataBase Cleaning routine.
#' 
#' This function, with a VMS database and two shape files with land polygon and harbours points,
#'  performs a filtered search over the whole db assigning warning status to the vms raw data.
#'  
#' @param vms_db_name The path of a VMS DataBase
#' @param map_file_name The path of a shape file with land polygon data
#' @param harb_file_name The path of a shape file with harbours point data
#' 
#' @return This function does not return a value. 
#' 
#' @usage gui_vms_db_clean(vms_db_name = "", map_file_name = "", harb_file_name = "")
#' 
#' @export gui_vms_db_clean
#'

gui_vms_db_clean <- function(vms_db_name = "", map_file_name = "", harb_file_name = "")
{
  vms_DB <- vms_DB$new()
  vms_DB$db <- vms_db_name
  
  themap <- polymap$new()
  themap$path = map_file_name
  
  harb <- harbCoo$new()
  harb$path = harb_file_name
  
  vms_db_clean_win <- gwindow("VMS Data Cleaning Utility", visible = FALSE)
  
  # VMS DB check for pings in harbour, coherence and on land
  chk_g <- gframe(horizontal = FALSE, container = vms_db_clean_win)
  chk_g2 <- ggroup(horizontal = TRUE, container = chk_g)
  chk_g3 <- ggroup(horizontal = FALSE, container = chk_g)
  chk_g3a <- ggroup(horizontal = TRUE, container = chk_g3)
  chk_g3b <- ggroup(horizontal = TRUE, container = chk_g3)
  chk_g3c <- ggroup(horizontal = TRUE, container = chk_g3)  
  
  addSpring(chk_g2)
  gimage(system.file("ico/edit-clear-2.png", package = "vmsbase"), container = chk_g2)
  proglab <- glabel("Ping check" , container = chk_g2)
  addSpring(chk_g2)
  
  #################
  addSpring(chk_g3)
  addSpring(chk_g3a)
  vms_db_f <- gframe(text = "VMS DB file", horizontal = TRUE, container = chk_g3a)
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
           if(themap$path != "" & harb$path != "" & vms_DB$db != "")
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
  addSpring(chk_g3a)
  addSpring(chk_g3)
  ################
  addSpring(chk_g3b)
  cus_map_g <- gframe(text = "Land Shape File", horizontal = TRUE, container = chk_g3b)
  addSpring(cus_map_g)
  cus_map_lab <- glabel("Select Land Shape File", container = cus_map_g)
  addSpring(cus_map_g)
  gimage(system.file("ico/folder-html.png", package="vmsbase"), container = cus_map_g,
         handler = function(h,...){
           themap$path <- gfile(text = "Select Land ShapePoly map",
                                type = "open",
                                filter = list("shp data" = list(patterns = c("*.shp"))))
           #            svalue(cus_map_lab) <- paste("Land: ", strsplit(themap$path, "/")[[1]][length(strsplit(themap$path, "/")[[1]])], sep = "")
           svalue(cus_map_lab) <- paste("Land: ", ifelse(.Platform$OS.type == "windows", strsplit(themap$path, "\\\\")[[1]][length(strsplit(themap$path, "\\\\")[[1]])],strsplit(themap$path, "/")[[1]][length(strsplit(themap$path, "/")[[1]])]), sep = "")
           if(themap$path != "" & harb$path != "" & vms_DB$db != "")
           {
             enabled(start_b) <- TRUE
           }
         })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = cus_map_g,
         handler = function(h,...){
           themap$path <- ""
           svalue(cus_map_lab) <- "Select Land Shape File"
         })
  addSpring(chk_g3b)
  
  ##Harbours file
  addSpring(chk_g3c)
  cus_har_g <- gframe(text = "Harbours Shape File", horizontal = TRUE, container = chk_g3c)
  addSpring(cus_har_g)
  cus_har_lab <- glabel("Select Harbours Shape File", container = cus_har_g)
  addSpring(cus_har_g)
  gimage(system.file("ico/folder-man.png", package="vmsbase"), container = cus_har_g,
         handler = function(h,...){
           harb$path <- gfile(text = "Select ShapePoints map",
                              type = "open",
                              filter = list("shp data" = list(patterns = c("*.shp"))))
           #            svalue(cus_har_lab) <- paste("Harbour: ", strsplit(harb$path, "/")[[1]][length(strsplit(harb$path, "/")[[1]])], sep = "")
           svalue(cus_har_lab) <- paste("Harbour: ", ifelse(.Platform$OS.type == "windows", strsplit(harb$path, "\\\\")[[1]][length(strsplit(harb$path, "\\\\")[[1]])],strsplit(harb$path, "/")[[1]][length(strsplit(harb$path, "/")[[1]])]), sep = "")
           if(themap$path != "" & harb$path != "" & vms_DB$db != "")
           {
             enabled(start_b) <- TRUE
           }
         })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = cus_har_g,
         handler = function(h,...){
           harb$path <- ""
           svalue(cus_har_lab) <- "Select Land Shape File"
         })
  addSpring(chk_g3c)
  addSpring(chk_g)
  infolab <- glabel("" , container = chk_g)
  addSpring(chk_g)
  start_b <- gbutton("Start cleaning", container = chk_g, handler = function(h,...)
  {
    enabled(cus_map_g) <- FALSE
    enabled(cus_har_g) <- FALSE
    enabled(vms_db_f) <- FALSE
    enabled(start_b) <- FALSE
    
    svalue(infolab) <- "Loading...\n Standard Harbours"
    #     harb$path <- system.file("shp/harb_it.shp", package="vmsbase")
    harb$data <- readShapePoints(harb$path)
    XCOORD <- harb$data@coords[,1]
    YCOORD <- harb$data@coords[,2]
    
    svalue(infolab) <- "Loading...\nLand Map"
    #     themap$path <- system.file("shp/Med_Poly.shp", package="vmsbase")
    themap$data <- readShapePoly(themap$path)
    
    svalue(infolab) <- "Updating...\nVMS DataBase"
    sqldf("drop table if exists warn", dbname = vms_DB$db)
    sqldf("CREATE TABLE warn(p_id INT, W_DUPL INT, W_HARB INT, W_LAND INT, W_COHE INT)", dbname = vms_DB$db)
    
    svalue(infolab) <- "Ping Cleaning Started"
    cat("\n\n   ---   Ping Cleaning Started   ---\n")
    incee <- sqldf("select distinct I_NCEE from ping", dbname = vms_DB$db)
    for ( i in 1:nrow(incee) )
    {
      svalue(infolab) <- paste("Processing\nVessel: ", i," of ", nrow(incee), spe = "")
      cat("\nVessel: ", i," of ", nrow(incee), spe = "")
      
      vessel <- fn$sqldf("select ROWID, * from ping where I_NCEE = `incee[i,1]` order by DATE", dbname = vms_DB$db)
      
      numlines <- nrow(vessel)
      cat(" with ", numlines, " pings", sep = "")
      if(numlines == 0)
      {
        cat(" - Skipped!", sep = "")
        next
      }else{
        ann_data <- data.frame("ROWID" = numeric(numlines), 
                               "W_DUPL" = NA, 
                               "W_HARB" = NA, 
                               "W_LAND" = NA, 
                               "W_COHE" = integer(numlines))
        
        ann_data["ROWID"] <- vessel["rowid"]
        
        # Check duplicated pings
        cat("\nChecking duplicated pings", sep = "")
        ann_data["W_DUPL"] <- duplicated(vessel[,2:7])
        cat(" - Completed!", sep = "")
        
        cat("\nChecking pings in harbour and coherence", sep = "")
        for (j in 1:numlines)
        {
          #         svalue(proglab) <- paste("Vessel: ", i," of ", nrow(incee), spe = "")      
          
          #Check pings in harbour
          hdist <- min(spDistsN1(cbind(XCOORD, YCOORD), as.matrix(c(vessel[j,"LON"],vessel[j,"LAT"])), longlat = TRUE))
          ifelse(hdist < 2, ann_data[j,"W_HARB"] <- T, ann_data[j,"W_HARB"] <- F)
          
          #Check pings coherence
          if(j == 1)
          {
            if(nrow(ann_data) > 1)
            {
              if(ann_data[j+1,"W_DUPL"] == F)
              {
                succdist <- spDists(matrix(c(vessel[j,"LON"],vessel[j,"LAT"]), ncol = 2), matrix(c(vessel[j+1,"LON"],vessel[j+1,"LAT"]), ncol = 2), longlat = TRUE)
                succlag <- (vessel[j+1,"DATE"]-vessel[j,"DATE"])*24
                if(succlag == 0) {succlag <- 0.01}
                succvel <- succdist/succlag
                ifelse(succvel < 50, ann_data[j,"W_COHE"] <- 2, ann_data[j,"W_COHE"] <- 0)
              }
              else{
                ann_data[j,"W_COHE"] <- 2
              }
            }
            else{
              ann_data[j,"W_COHE"] <- 3
            }
            next
          }
          
          if(j == numlines)
          {
            if(ann_data[j,"W_DUPL"] == F)
            {
              predist <- spDists(matrix(c(vessel[j,"LON"],vessel[j,"LAT"]), ncol = 2), matrix(c(vessel[j-1,"LON"],vessel[j-1,"LAT"]), ncol = 2), longlat = TRUE)
              prelag <- (vessel[j,"DATE"]-vessel[j-1,"DATE"])*24
              if(prelag == 0) {prelag <- 0.01}
              prevel <- predist/prelag
              ifelse(prevel < 50, ann_data[j,"W_COHE"] <- 1, ann_data[j,"W_COHE"] <- 0)
            }
            else{
              ann_data[j,"W_COHE"] <- 1
            }
            next
          }
          
          if (j > 1 & j < numlines)
          {
            if(ann_data[j,"W_DUPL"] == T)
            {
              ann_data[j,"W_COHE"] <- ann_data[j-1,"W_COHE"]
              next
            }
            
            if(ann_data[j-1,"W_DUPL"] == T)
            {
              succdist <- spDists(matrix(c(vessel[j,"LON"],vessel[j,"LAT"]), ncol = 2), matrix(c(vessel[j+1,"LON"],vessel[j+1,"LAT"]), ncol = 2), longlat = TRUE)
              succlag <- (vessel[j+1,"DATE"]-vessel[j,"DATE"])*24
              if(succlag == 0) {succlag <- 0.01}
              succvel <- succdist/succlag
              ifelse(succvel < 50, ann_data[j,"W_COHE"] <- 2, ann_data[j,"W_COHE"] <- 0)
              next
            }
            
            if(ann_data[j+1,"W_DUPL"] == T)
            {
              predist <- spDists(matrix(c(vessel[j,"LON"],vessel[j,"LAT"]), ncol = 2), matrix(c(vessel[j-1,"LON"],vessel[j-1,"LAT"]), ncol = 2), longlat = TRUE)
              prelag <- (vessel[j,"DATE"]-vessel[j-1,"DATE"])*24
              if(prelag == 0) {prelag <- 0.01}
              prevel <- predist/prelag
              ifelse(prevel < 50, ann_data[j,"W_COHE"] <- 1, ann_data[j,"W_COHE"] <- 0)
              next
            }
            
            succdist <- spDists(matrix(c(vessel[j,"LON"],vessel[j,"LAT"]), ncol = 2), matrix(c(vessel[j+1,"LON"],vessel[j+1,"LAT"]), ncol = 2), longlat = TRUE)
            succlag <- (vessel[j+1,"DATE"]-vessel[j,"DATE"])*24
            if(succlag == 0) {succlag <- 0.01}
            succvel <- succdist/succlag
            
            predist <- spDists(matrix(c(vessel[j,"LON"],vessel[j,"LAT"]), ncol = 2), matrix(c(vessel[j-1,"LON"],vessel[j-1,"LAT"]), ncol = 2), longlat = TRUE)
            prelag <- (vessel[j,"DATE"]-vessel[j-1,"DATE"])*24
            if(prelag == 0) {prelag <- 0.01}
            prevel <- predist/prelag
            
            if(prevel < 50)
            {
              ann_data[j,"W_COHE"] <- 1
            }
            if(succvel < 50)
            {
              ann_data[j,"W_COHE"] <- 2
            }
            if(prevel < 50 & succvel < 50)
            {
              ann_data[j,"W_COHE"] <- 3
            }
            if(prevel > 50 & succvel > 50)
            {
              ann_data[j,"W_COHE"] <- 0
            }
            
          }
        }
        cat(" - Completed!", sep = "")
        
        #Check pings on land
        
        cat("\nChecking pings on land", sep = "")
        
        onland <- over(SpatialPoints(c(vessel["LON"], vessel["LAT"])), themap$data)
        olpts <- which(!is.na(onland[,1]))
        nlpts <- which(is.na(onland[,1]))
        
        ann_data[olpts,"W_LAND"] <- T
        ann_data[nlpts,"W_LAND"] <- F
        
        cat(" - Completed!\n", sep = "")
        
        sqldf("insert into warn select * from `ann_data`", dbname = vms_DB$db)
        
      }
    }
    
    cat("\n\n   ---   End Ping Cleaning   ---\n", sep = "")
    
    tot_p <- 100/sqldf("select count(*) from warn", dbname = vms_DB$db)[1,]
    no_p1 <- sqldf("select count(*) from warn where W_DUPL = 1", dbname = vms_DB$db)[1,]
    no_p2 <- sqldf("select count(*) from warn where W_HARB = 1", dbname = vms_DB$db)[1,]
    no_p3 <- sqldf("select count(*) from warn where W_LAND = 1", dbname = vms_DB$db)[1,]
    no_p4 <- sqldf("select count(*) from warn where W_COHE = 0", dbname = vms_DB$db)[1,]
    
    cat("\n - Found ", no_p1, " (", round(tot_p * no_p1, 2), "% of total)", " duplicated pings!", sep = "")
    cat("\n - Found ", no_p2, " (", round(tot_p * no_p2, 2), "% of total)", " pings in harbour!", sep = "")
    cat("\n - Found ", no_p3, " (", round(tot_p * no_p3, 2), "% of total)", " pings on land!", sep = "")
    cat("\n - Found ", no_p4, " (", round(tot_p * no_p4, 2), "% of total)", " not coherent pings!\n\n", sep = "")
    
    gconfirm("VMS DB Ping check Completed!",
             title = "Confirm",
             icon = "info",
             parent = vms_db_clean_win,
             handler = function(h,...){dispose(vms_db_clean_win)})
  })
  enabled(start_b) <- FALSE
  
  if(vms_DB$db != "")
  {
    #     svalue(sel_vms_f) <- strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])]    
    svalue(sel_vms_f) <- ifelse(.Platform$OS.type == "windows", strsplit(vms_DB$db, "\\\\")[[1]][length(strsplit(vms_DB$db, "\\\\")[[1]])],strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])])
  }
  if(themap$path != "")
  {
    #     svalue(cus_map_lab) <- paste("Land: ", strsplit(themap$path, "/")[[1]][length(strsplit(themap$path, "/")[[1]])], sep = "")
    svalue(cus_map_lab) <- paste("Land: ", ifelse(.Platform$OS.type == "windows", strsplit(themap$path, "\\\\")[[1]][length(strsplit(themap$path, "\\\\")[[1]])],strsplit(themap$path, "/")[[1]][length(strsplit(themap$path, "/")[[1]])]), sep = "")
    
  }
  if(harb$path != "")
  {
    #     svalue(cus_har_lab) <- paste("Harbour: ", strsplit(harb$path, "/")[[1]][length(strsplit(harb$path, "/")[[1]])], sep = "")
    svalue(cus_har_lab) <- paste("Harbour: ", ifelse(.Platform$OS.type == "windows", strsplit(harb$path, "\\\\")[[1]][length(strsplit(harb$path, "\\\\")[[1]])],strsplit(harb$path, "/")[[1]][length(strsplit(harb$path, "/")[[1]])]), sep = "")
  }
  if(themap$path != "" & harb$path != "" & vms_DB$db != "")
  {
    enabled(start_b) <- TRUE
  }
  
  visible(vms_db_clean_win) <- TRUE
}