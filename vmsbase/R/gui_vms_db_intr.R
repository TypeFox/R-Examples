
#' VMS DB Interpolation GUI
#'  
#' The \code{gui_vms_db_intr} function implements the graphical user interface for the
#'  VMS Track Interpolation routine.
#' 
#' This function, with a VMS cleaned and cutted database (see \code{\link{gui_vms_db_cut}} and \code{\link{gui_vms_db_clean}}),
#'  interpolates track data of each vessel in the VMS DB.
#'  
#' @param vms_db_name The path of a VMS DataBase
#' 
#' @return This function does not return a value. 
#' 
#' @usage gui_vms_db_intr(vms_db_name = "")
#' 
#' @export gui_vms_db_intr
#'
#' @seealso \code{\link{gui_vms_db_cut}} \code{\link{gui_vms_db_clean}}
#'
#' @references 
#' Russo, T., Parisi, A. and Cataudella, S. (2011) New insights in interpolating fishing tracks from VMS data for different metiers. \emph{Fisheries Research}, \bold{108(1)}, 184--194.
#' \url{http://www.sciencedirect.com/science/article/pii/S0165783610003450}


gui_vms_db_intr <- function(vms_db_name = "")
{
  vms_DB <- vms_DB$new()
  vms_DB$db <- vms_db_name
  
  vms_db_intr_win <- gwindow("VMS Data Interpolation Utility", visible = FALSE)
  
  # TRACK InterpolatioN
  int_g <- gframe(horizontal = FALSE, container = vms_db_intr_win)
  int_g2 <- ggroup(horizontal = TRUE, container = int_g)
  int_g3 <- ggroup(horizontal = TRUE, container = int_g)
  addSpring(int_g2)
  gimage(system.file("ico/draw-bezier-curves.png", package="vmsbase"), container = int_g2)
  proglab_int <- glabel("Track Interpolation" , container = int_g2)
  addSpring(int_g2)
  
  
  #################
  addSpring(int_g3)
  vms_db_f <- gframe(text = "VMS DB file", horizontal = TRUE, container = int_g3)
  addSpring(vms_db_f)
  sel_vms_f <- glabel("Select VMS DB file", container = vms_db_f)
  addSpring(vms_db_f)
  gimage(system.file("ico/folder-blue.png", package="vmsbase"), container = vms_db_f,
         handler = function(h,...){
           vms_DB$db <- gfile(text = "Select VMS DataBase file",
                              type = "open",
                              filter = list("VMS DB file" = list(patterns = c("*.vms.sqlite"))))
           svalue(sel_vms_f) <- ifelse(.Platform$OS.type == "windows", strsplit(vms_DB$db, "\\\\")[[1]][length(strsplit(vms_DB$db, "\\\\")[[1]])],strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])])
           enabled(start_b) <- TRUE
         })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = vms_db_f,
         handler = function(h,...){
           vms_DB$db <- ""
           enabled(start_b) <- FALSE
           svalue(sel_vms_f) <- "Select VMS DB file"
         })
  addSpring(int_g3)
  ################
  
  addSpring(int_g)
  infolab_int <- glabel("" , container = int_g)
  addSpring(int_g)
  int_g3 <- ggroup(horizontal = TRUE, container = int_g)
  addSpring(int_g3)
  glabel("Interpolate every" , container = int_g3)
  int_fre_cho <- gspinbutton(from = 1, to = 60, by = 1, value = 10, horizontal = TRUE, container = int_g3)
  glabel("minutes" , container = int_g3)
  addSpring(int_g3)
  start_b <- gbutton("Start Interpolation", container = int_g, handler = function(h,...)
  {
    enabled(int_g3) <- FALSE
    enabled(start_b) <- FALSE
    enabled(vms_db_f) <- FALSE
    
    if(sqldf("select count(*) from track", dbname = vms_DB$db) > 0)
    {
      
      int_fre <- svalue(int_fre_cho)
      
      svalue(infolab_int) <- "Updating...\nVMS DataBase"
      
      sqldf("drop table if exists intrp", dbname = vms_DB$db)
      
      sqldf("CREATE TABLE intrp(I_NCEE INT, LAT REAL, LON REAL, DATE REAL, SPE REAL, HEA REAL, W_HARB INT, T_NUM INT, P_ID INT, P_INT INT, T_ID INT)", dbname = vms_DB$db)
      
      cat("\n   ---   Interpolation Started   ---\n", sep = "")
      
      incee <- sqldf("select distinct I_NCEE from track", dbname = vms_DB$db)
      
      num_incee <- nrow(incee)
      
      for ( v in 1:num_incee )
      {
        vessel <- fn$sqldf("select ROWID, * from track where I_NCEE = `incee[v,1]` order by DATE", dbname = vms_DB$db)
        
        if(nrow(vessel) == 0)
        {
          
          cat(" - Skipped, no pings", sep = "")
          
        }else{
          maxtrk_n <- max(vessel["T_NUM"])
          svalue(infolab_int ) <- paste("\nVessel: ", v, " of ", num_incee, sep = "")
          cat("\nVessel: ", incee[v,1], " - ",  v, " of ", num_incee, " with ", maxtrk_n, " tracks\nInterpolating...", sep = "")
          
          for(ind in 1:maxtrk_n)
          {
            #        svalue(infolab3 ) <- paste("\nVessel: ", v," of ", num_incee, "\nTrack: ", ind, " of ", maxtrk_n, spe = "")
            cat(".", sep = "")
            tracks <- vessel[which(vessel["T_NUM"] == ind),]
            tracks <-tracks[!is.na(tracks["DATE"]),]
            tracki <- tracks
            
            numrow <- nrow(tracki)
            
            if(numrow > 4)
            {
              
              COO_LL <- as.data.frame(cbind(tracki[,"LON"],tracki[,"LAT"]))
              sign_chk <- FALSE
              if(abs(sum(sign(COO_LL[,1]))) != nrow(COO_LL)){
                COO_LL[which(COO_LL[,1] < 0),1] <- COO_LL[which(COO_LL[,1] < 0),1] + 360
                sign_chk <- TRUE
              }
              
              
              colnames(COO_LL) <- c("X","Y")
              attr(COO_LL, "projection") <- "LL"
              # utm_zone <- attr(suppressMessages(convUL(COO_LL)), "zone")
              # utm_zone <- min(long2UTM(COO_LL[,1]))
              mid_poi <- min(COO_LL[,1])+((max(COO_LL[,1])-min(COO_LL[,1]))/2)
              utm_zone <- long2UTM(mid_poi)
              SP_UTM <- spTransform(SpatialPoints(COO_LL,
                                                  proj4string=CRS("+proj=longlat +datum=WGS84")),
                                    CRS(paste("+proj=utm +zone=",
                                              utm_zone,
                                              "+datum=WGS84", sep = "")))
              num_ro <- nrow(SP_UTM@coords)
              COO_UTM <- data.frame("X" = numeric(num_ro),
                                    "Y" = numeric(num_ro),
                                    "zone" = numeric(num_ro))
              COO_UTM$X <- coordinates(SP_UTM)[,1]/1000
              COO_UTM$Y <- coordinates(SP_UTM)[,2]/1000
              
              tracki["LON"] <- COO_UTM[,"X"]
              tracki["LAT"] <- COO_UTM[,"Y"]
              hp <- data.frame("LON" = numeric(numrow), "LAT" = numeric(numrow))
              hcrm <- data.frame("LON" = numeric(numrow), "LAT" = numeric(numrow))
              hdri <- data.frame("LON" = numeric(numrow), "LAT" = numeric(numrow))
              hdrim <- data.frame("LON" = numeric(1), "LAT" = numeric(1))
              tmai <- data.frame("LON" = numeric(numrow), "LAT" = numeric(numrow))
              elle <- numeric(numrow)
              NAXNA = tracki[,"LON"]
              NAYNA = tracki[,"LAT"]
              dts = diff(tracki[,"DATE"])
              dts[which(dts == 0)] <- dts[which(dts == 0)] + 0.0000001
              dts <- dts * 24
              ltrack = nrow(tracki)-2
              temp = 0.5*(cbind(diff(NAXNA)[-1],diff(NAYNA)[-1])/dts[-1] + cbind(diff(NAXNA)[-(ltrack+1)],diff(NAYNA)[-(ltrack+1)])/dts[-(ltrack+1)])    
              hcr <- as.data.frame(rbind(rep(0,2),temp,rep(0,2)))
              colnames(hcr) <- c("LON","LAT")
              oldnorms = sqrt(apply(hcr^2,1,sum))
              distances = sqrt((NAXNA[-1]-NAXNA[-length(NAXNA)])^2 + (NAYNA[-1]-NAYNA[-length(NAYNA)])^2)
              newnorms = c(NA,(distances[-length(distances)] + distances[-1])/(dts[-length(dts)]+dts[-1]),NA)
              hcrm = as.data.frame(matrix(rep(newnorms/oldnorms,2),ncol=2)*hcr)
              colnames(hcrm) <- c("LON","LAT")
              hp["LON"] <- tracki["SPE"] * cos((450-tracki["HEA"])*(2*pi)/360)
              hp["LAT"] <- tracki["SPE"] * sin((450-tracki["HEA"])*(2*pi)/360)
              hdri["LON"] <- hcrm["LON"] - hp["LON"]
              hdri["LAT"] <- hcrm["LAT"] - hp["LAT"]
              hdrim["LON"] <- median(hdri[which(!is.na(hdri["LON"])),"LON"])
              hdrim["LAT"] <- median(hdri[which(!is.na(hdri["LAT"])),"LAT"])
              tmai["LON"] <- (hdrim[1,"LON"] + hp["LON"])
              tmai["LAT"] <- (hdrim[1,"LAT"] + hp["LAT"])
              temp_inte <- 0.0006944445*int_fre
              tms <- 1:(ceiling(1/temp_inte) * (ceiling(tracks[nrow(tracks),"DATE"]) - floor(tracks[1,"DATE"]))) * (temp_inte)
              floo <- floor(tracks[1,"DATE"] )
              dec <- tracks["DATE"] - floo
              int <- sort(c(tms[which(tms > min(dec[,"DATE"]) & tms < max(dec[,"DATE"]))], dec[,"DATE"]))
              numpoi <- length(int)
              newpoi <- data.frame("I_NCEE" = numeric(numpoi),
                                   "LAT" = numeric(numpoi),
                                   "LON" = numeric(numpoi),
                                   "DATE" = numeric(numpoi),
                                   "SPE" = numeric(numpoi),
                                   "HEA" = numeric(numpoi),
                                   "W_HARB" = integer(numpoi),
                                   "T_NUM" = numeric(numpoi),
                                   "P_ID" = numeric(numpoi),
                                   "P_INT" = numeric(numpoi),
                                   "T_ID" = numeric(numpoi))
              newpoi["T_ID"] <- NA
              newpoi["P_ID"] <- NA
              newpoi["I_NCEE"] <- tracki[1,"I_NCEE"]
              newpoi["T_NUM"] <- ind
              newpoi["DATE"] <- int + floo
              newpoi["P_INT"] <- 0
              newpoi[1,"LON"] <- tracki[1,"LON"]
              newpoi[1,"LAT"] <- tracki[1,"LAT"]
              newpoi[1,"SPE"] <- tracki[1,"SPE"]
              newpoi[1,"HEA"] <- tracki[1,"HEA"]
              newpoi[1,"W_HARB"] <- tracki[1,"W_HARB"]
              newpoi[1,"T_ID"] <- tracki[1,"rowid"]
              newpoi[1,"P_ID"] <- tracki[1,"P_ID"]
              newpoi[nrow(newpoi),"LON"] <- tracki[nrow(tracki),"LON"]
              newpoi[nrow(newpoi),"LAT"] <- tracki[nrow(tracki),"LAT"]
              newpoi[nrow(newpoi),"SPE"] <- tracki[nrow(tracki),"SPE"]
              newpoi[nrow(newpoi),"HEA"] <- tracki[nrow(tracki),"HEA"]
              newpoi[nrow(newpoi),"W_HARB"] <- tracki[nrow(tracki),"W_HARB"]
              newpoi[nrow(newpoi),"T_ID"] <- tracki[nrow(newpoi),"rowid"]
              newpoi[nrow(newpoi),"P_ID"] <- tracki[nrow(newpoi),"P_ID"]
              dup <- which(newpoi[,"DATE"] %in% tracks[,"DATE"])
              newpoi[dup,"LON"] <-  tracki[which(tracks["DATE"] == newpoi[dup,"DATE"]),"LON"]
              newpoi[dup,"LAT"] <-  tracki[which(tracks["DATE"] == newpoi[dup,"DATE"]),"LAT"]
              newpoi[dup,"SPE"] <-  tracki[which(tracks["DATE"] == newpoi[dup,"DATE"]),"SPE"]
              newpoi[dup,"HEA"] <-  tracki[which(tracks["DATE"] == newpoi[dup,"DATE"]),"HEA"]
              newpoi[dup,"T_ID"] <- tracki[which(tracks["DATE"] == newpoi[dup,"DATE"]),"rowid"]
              newpoi[dup,"P_ID"] <- tracki[which(tracks["DATE"] == newpoi[dup,"DATE"]),"P_ID"]
              
              for(que in 1:(nrow(dec)-1))
              {
                sus <- which(int > dec[que,"DATE"] & int < dec[que+1,"DATE"])           
                dimen <- length(sus)
                if(dimen > 0)
                {
                  if(tracki[que + 1,"LON"]  == tracki[que,"LON"] & tracki[que + 1,"LAT"]  == tracki[que,"LAT"])
                  {
                    newpoi[sus, "LON"] <- tracki[que,"LON"]
                    newpoi[sus, "LAT"] <- tracki[que,"LAT"]
                    newpoi[sus,"P_INT"] <- 1
                  }
                  else
                  {
                    if(que == 1 | que == (nrow(dec)-1))
                    {
                      for(repe in 1:dimen)
                      {
                        essei <- data.frame("LON" = numeric(1), "LAT" = numeric(1))
                        deltati <- (dec[que + 1,"DATE"] - dec[que,"DATE"])*24
                        coe <- (int[sus[repe]] - dec[que,"DATE"])*24/deltati
                        xi <- cbind(tracki[que,"LON"], tracki[que,"LAT"])
                        xip1 <- cbind(tracki[que+1,"LON"], tracki[que+1,"LAT"])
                        essei["LON"] <- xi[1,1] + (coe * (xip1[1,1]-xi[1,1]))
                        essei["LAT"] <- xi[1,2] + (coe * (xip1[1,2]-xi[1,2]))
                        newpoi[sus[repe],"LON"] <- essei["LON"]
                        newpoi[sus[repe],"LAT"] <- essei["LAT"]
                        newpoi[sus[repe],"P_INT"] <- 1
                      }
                    }else{
                      for(repe in 1:dimen)
                      {
                        essei <- data.frame("LON" = numeric(1), "LAT" = numeric(1))
                        difft = (tracki[que + 1,"DATE"] - tracki[que,"DATE"])*24
                        invdifft = 1/(difft)
                        normTimes = ((newpoi[sus[repe],"DATE"]-tracki[que,"DATE"])*24)/difft
                        xi <- cbind(tracki[que,"LON"], tracki[que,"LAT"])
                        xip1 <- cbind(tracki[que+1,"LON"], tracki[que+1,"LAT"])
                        essei["LON"] = H0(normTimes)*xi[1] + H1(normTimes)*xip1[1] + H2(normTimes)*difft*tmai[que,"LON"]+ H3(normTimes)*difft*tmai[que+1,"LON"]
                        essei["LAT"] = H0(normTimes)*xi[2] + H1(normTimes)*xip1[2] + H2(normTimes)*difft*tmai[que,"LAT"] + H3(normTimes)*difft*tmai[que+1,"LAT"]
                        dX = invdifft*dH0(normTimes)*xi[1] + invdifft*dH1(normTimes)*xip1[1] + dH2(normTimes)*tmai[que,"LON"] + dH3(normTimes)*tmai[que+1,"LON"]
                        dY = invdifft*dH0(normTimes)*xi[2] + invdifft*dH1(normTimes)*xip1[2] + dH2(normTimes)*tmai[que,"LAT"] + dH3(normTimes)*tmai[que+1,"LAT"]
                        norms = sqrt(dX^2 + dY^2)
                        thetas=numeric(1)
                        thetas[dX>=0]=atan(dY[dX>=0]/dX[dX>=0])
                        thetas[dX<0]= -atan(dY[dX<0]/dX[dX<0])+pi
                        thetas[(dX==0)&(dY==0)]=0
                        newpoi[sus[repe],"LON"] <- essei["LON"]
                        newpoi[sus[repe],"LAT"] <- essei["LAT"]
                        newpoi[sus[repe],"SPE"] <- norms
                        new_hea <- 450-((thetas)*180/pi)
                        newpoi[sus[repe],"HEA"] <- ifelse(new_hea < 360, new_hea, new_hea - 360)
                        newpoi[sus[repe],"P_INT"] <- 1
                      }
                    }
                  }
                }
              }
              
              utm <- as.data.frame(cbind(newpoi["LON"], newpoi["LAT"]))
              utm_poi <- SpatialPoints(utm*1000, CRS(paste("+proj=utm +zone=",
                                                           utm_zone,
                                                           "+datum=WGS84", sep = "")))
              ll <- spTransform(utm_poi, CRS("+proj=longlat +datum=WGS84"))
              
              newpoi[,c("LON","LAT")] <- coordinates(ll)
              
              if(sign_chk){
                newpoi[which(newpoi[,1] > 180),1] <- newpoi[which(newpoi[,1] > 180),1] - 360
              }
              
              sqldf("insert into intrp select * from `newpoi`", dbname = vms_DB$db)
            }
            
          }
          cat(" Completed!\n")
        }
      }
      
      cat("\n   ---   End Interpolation   --- \n\n")
      gconfirm("VMS DB Track interpolation completed!",
               title = "Confirm",
               icon = "info",
               parent = vms_db_intr_win,
               handler = function(h,...){dispose(vms_db_intr_win)})
      
    }else{
      
      gconfirm("Track data not available\n\nExecute Track Cut first!",
               title = "Error",
               icon = "error",
               parent = vms_db_intr_win,
               handler = function(h,...){dispose(vms_db_intr_win)})
      
    }
    
  })
  enabled(start_b) <- FALSE
  
  if(vms_DB$db != "")
  {
    svalue(sel_vms_f) <- ifelse(.Platform$OS.type == "windows", strsplit(vms_DB$db, "\\\\")[[1]][length(strsplit(vms_DB$db, "\\\\")[[1]])],strsplit(vms_DB$db, "/")[[1]][length(strsplit(vms_DB$db, "/")[[1]])])
    enabled(start_b) <- TRUE
  }
  
  visible(vms_db_intr_win) <- TRUE
}