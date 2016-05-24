
#' LogBook Metier Discovery GUI
#'  
#' 
#' The \code{gui_lb_met_dis} function implements the graphical user interface for the
#'  LogBook Metier Discovery
#' 
#' This function,  with a logbook database, performs a clustering over the whole
#'  db calculating the distance or dissimilarity between every observation
#'   in the sample.
#'   
#' @param lb_db_name The path of a LogBook DataBase
#' 
#' @return This function does not return a value. 
#' After the execution the user is asked where to save the result file.
#' 
#' @usage gui_lb_met_dis(lb_db_name = "")
#' 
#' @export gui_lb_met_dis  
#'

gui_lb_met_dis <- function(lb_db_name = "")
{
  
  lb_DB <- log_DB$new()
  lb_DB$db <- lb_db_name
  lb_CLA <- log_Cla$new()
  
  met_dis_win <- gwindow("Metier Discovery Tool", visible = FALSE)
  big_g <- ggroup(horizontal = FALSE, container = met_dis_win)
  
  new_g <- ggroup(horizontal = TRUE, container = big_g)
  addSpring(new_g)
  lb_db_f <- gframe(text = "LogBook DB file", horizontal = TRUE, container = new_g)
  addSpring(lb_db_f)
  sel_lb_f <- glabel("Select LB DB file", container = lb_db_f)
  addSpring(lb_db_f)
  gimage(system.file("ico/folder-orange.png", package="vmsbase"), container = lb_db_f,
         handler = function(h,...){
           lb_DB$db <<- gfile(text = "Select LB DataBase file",
                              type = "open",
                              filter = list("LB DB file" = list(patterns = c("*.lb.sqlite"))))
           svalue(sel_lb_f) <- ifelse(.Platform$OS.type == "windows", strsplit(lb_DB$db, "\\\\")[[1]][length(strsplit(lb_DB$db, "\\\\")[[1]])],strsplit(lb_DB$db, "/")[[1]][length(strsplit(lb_DB$db, "/")[[1]])])
           enabled(start_b) <- TRUE
           enabled(g_sup) <- TRUE
         })
  gimage(system.file("ico/application-exit-5.png", package="vmsbase"), container = lb_db_f,
         handler = function(h,...){
           lb_DB$db <<- ""
           enabled(start_b) <- FALSE
           enabled(g_sup) <- FALSE
           svalue(sel_lb_f) <- "Select LB DB file"
         })
  addSpring(new_g)
  
  g_sup <- ggroup(horizontal = TRUE, container = big_g)
  addSpring(g_sup)
  
  g_s_b <- gframe("Dataset Cleaning", horizontal = F, container = g_sup)
  cles_ex <- gcheckbox("Remove Exotics", checked = TRUE, container = g_s_b)
  
  cles_ef_g <- ggroup(horizontal = TRUE, container = g_s_b)
  cles_ef <- gcheckbox("Species f <", checked = TRUE, handler = function(h,...){enabled(cles_ef_q) <- !enabled(cles_ef_q)}, container = cles_ef_g)
  cles_ef_q <- gspinbutton (from = 0, to = 3, by = 0.001, value = 0, container = cles_ef_g)
  
  cles_ou <- gcheckbox("Species Outliers", checked = TRUE, container = g_s_b)
  cles_oul <- gcheckbox("Logs Outliers", checked = TRUE, container = g_s_b)
  cles_lo <- gcheckbox("logarithm", checked = FALSE, container = g_s_b)
  cles_st <- gcheckbox("Stand by time", checked = FALSE, container = g_s_b)
  
  g_s_c <- gframe("Clustering", horizontal = F, container = g_sup)
  g_k_ran <- ggroup(horizontal = T, container = g_s_c)
  l_k_ran <- glabel("Clusters to test Min", container = g_k_ran)
  addSpring(g_k_ran)
  min_k_ran <- gspinbutton(from = 2, to = 60, by = 1, value = 1, digits = 0,
                           horizontal = TRUE, container = g_k_ran)
  g_k_ran2 <- ggroup(horizontal = T, container = g_s_c)
  
  l_k_ran2 <- glabel("Clusters to test Max", container = g_k_ran2)
  addSpring(g_k_ran2)
  max_k_ran <- gspinbutton(from = 2, to = 60, by = 1, value = 60, digits = 0,
                           horizontal = TRUE, container = g_k_ran2)
  
  g_samp <- ggroup(horizontal = T, container = g_s_c)
  l_samp <- glabel("Number of samples", container = g_samp)
  addSpring(g_samp)
  s_samp <- gspinbutton(from = 10, to = 1000, by = 10, value = 10, digits = 0,
                        horizontal = TRUE, container = g_samp)
  
  g_size <- ggroup(horizontal = T, container = g_s_c)
  l_size <- glabel("Size of samples", container = g_size)
  addSpring(g_size)
  s_size <- gspinbutton(from = 10, to = 1000, by = 10, value = 10, digits = 0,
                        horizontal = TRUE, container = g_size)
  
  g_metr <- gframe("Metric:", horizontal = FALSE, container = g_sup)
  s_metr <- gradio(c("euclidean", "manhattan"), selected = 2,
                   container = g_metr)
  g_subs <- gframe("Sub-Sampling:", horizontal = FALSE, container = g_sup)
  s_subs <- gcheckbox("N. SubSample", checked = TRUE, handler = function(h,...){enabled(subs_q) <- !enabled(subs_q)}, container = g_subs)
  subs_q <- gspinbutton (from = 100, to = 10000, by = 100, value = 1000, container = g_subs)
  
  g_go <- ggroup(horizontal = TRUE, container = g_sup)
  addSpring(g_go)
  start_b <- gbutton(text = "Start\nMetier Discovery", container = g_go, handler = function(h,...){
    enabled(g_sup) <- FALSE
    enabled(save_bu) <- FALSE
#     
#     if(exists("x") == TRUE){rm(x)}
#     
    cat("\n\nGetting the DataSet...")
    svalue(sb) <- "Getting the DataSet..."
    x <- sqldf("select * from elobo", dbname = lb_DB$db)
    
    if(svalue(cles_st) == TRUE)
    {
      cat("\nStandardizing by time...")
      svalue(sb) <- "Standardizing by time..."
      tms=ceiling(x[,3]-x[,2])
      tms.0 = which(tms==0)
      
      if(length(tms.0)>0){
        x = x[-tms.0,]
        tms = tms[-tms.0]
      }
      tms.outlier <- which(tms > 21)
      if(length(tms.outlier)>0){
        tms = tms[-tms.outlier]
        x=x[-tms.outlier,]
      }
      
      divid <- tms
      for(l in 1:(ncol(x)-6)){divid <- cbind(divid, tms)}
      x[,6:ncol(x)] <- x[,6:ncol(x)]/divid
      rm(divid)
    }
    
    cat("\nClearing Header...")
    svalue(sb) <- "Clearing Header..."
    x <- x[,-c(1,2,3,4,5)]
    

    if(svalue(cles_ex) == TRUE)
    {
      #Elimina le specie esotiche
      cat("\nExcluding exotic species")
      svalue(sb) <- "Excluding exotic species..."
      bad_str <- c("SQI", "CRQ", "GBO", "GGY", "OFJ", "COD", "SQL", "POI", "MAW")
      bad_str_FAO <- paste("FAO_", bad_str, sep = "")
      tole_ba <- which(colnames(x) %in% bad_str_FAO)
      if(length(tole_ba)>0){
        cat(" - ", length(tole_ba), " species removed", sep = " ")
        bad_str_col <- which(colnames(x) %in% bad_str_FAO)
        x <- x[,-bad_str_col]  
      } 
    }
    
    #Elimina le specie con f < 1x100
    if(svalue(cles_ef) == TRUE)
    {
      cat("\nExcluding species with f < ", svalue(cles_ef_q), sep = "")
      svalue(sb) <- paste("Removing species with f < ", svalue(cles_ef_q), "% ", sep = "")
      nposit <- function(vec) return(length(which(vec>0)))
      nocc <- 100*apply(x,2,nposit)/nrow(x)
      nocc_bru <- which(nocc <= svalue(cles_ef_q))
      if(length(nocc_bru)>0)
      {
        cat(" - ", length(nocc_bru)," species removed", sep = "")
        x <- x[,-nocc_bru]
      }
    }
    
    
    #Rimuove gli outliers
    if(svalue(cles_ou) == TRUE)
    {
      cat("\nExcluding Species Outliers")
      svalue(sb) <- "Excluding Species Outliers..."
      
      #################
#       x <- rm.outlier(x, fill = TRUE, median = TRUE, opposite = FALSE)
      
      tot_sp <- apply(x, 2, sum)
      sp_out <- boxplot(tot_sp, plot = FALSE)$out
      tole_sp <- which(colnames(x) %in% names(sp_out))
      if(length(tole_sp) > 0){
        cat(" - ", length(tole_sp), " species removed", sep = " ")
        x <- x[,-tole_sp]
      } 
    }
    
    if(svalue(cles_oul) == TRUE)
    {
      cat("\nExcluding Logs Outliers")
      svalue(sb) <- "Excluding Logs Outliers..."
      tot_ca <- apply(x, 1, sum)
      ca_out <- as.numeric(names(boxplot(tot_ca, plot = FALSE)$out))
      if(length(ca_out) > 0)
      {
        cat(" - ", length(ca_out), " logs removed", sep = "")
        x <- x[-ca_out,]
      }
    }
    
    if(svalue(s_subs) == TRUE & svalue(subs_q) < nrow(x))
    {
      totak <- sample(1:nrow(x), svalue(subs_q), replace = FALSE)
      x <- x[totak,]
    }
    
    if(svalue(cles_lo) == TRUE)
    {
      cat("\nCalculating logarithm", sep = "")
      svalue(sb) <- paste("Calculating logarithm", sep = "")
      x <- log(x+1)
    }
    
    #Ultimi controlli x scrupolo
    col_sum <- apply(x,2,sum)
    col_sum_zero <- which(col_sum == 0)
    if(length(col_sum_zero)>0) x <- x[,-col_sum_zero]
    row_sum <- apply(x,1,sum)
    row_sum_na <- which(is.na(row_sum)==TRUE)
    if(length(row_sum_na) > 0) x <- x[-row_sum_na,]
    
    if(nrow(x) > 0 & ncol(x) > 0)
    {    
    #Set k 
    krange = c(svalue(min_k_ran):min(nrow(x)-1, svalue(max_k_ran)))
    silhvec=numeric(length(krange))      
    #Perform the CLARA
    i=thr=0
    max_siz <- max(svalue(s_size), max(krange))
    cat("\n\n   -     Beginning Metier Discovery!     -\n\n - Testing partioning:\n")
    svalue(sb) <- "Metier Discovery..."
    tonorm <- ifelse((svalue(s_metr) == "Bray-Curtis"),FALSE,TRUE)
    S1=Sys.time()
    for(k in krange){
      if(max_siz <= k){
        cat("\n\n - STOP - Sample Size < number of clusters -")
        break}
      cat("\n\nTesting ", k," clusters - ", sep = "")
      svalue(sb) <- paste("Testing partioning with ",k," clusters", sep = "")
      S2=Sys.time()
      i=i+1

      xclara <- clara(x,
                        k,
                        metric = svalue(s_metr),
                        stand = tonorm,
                        samples = svalue(s_samp),
                        sampsize = min(nrow(x), max_siz),
                        trace = 0,
                        medoids.x = TRUE, keep.data = TRUE, rngR = TRUE, pamLike = TRUE)
      
      silhvec[i]=xclara$silinfo$avg.width
      
      cat("Average width: ", round(silhvec[i], 3), sep = "")
      svalue(sb) <- paste("Average width: ", silhvec[i], sep = "")
      
      if(silhvec[i] > thr){
        thr = silhvec[i]
        lb_CLA$data <- xclara
        
        cat(" - New Best Result!", sep = "")
        svalue(sb) <- paste("Cluster ", i," with average width: ", round(silhvec[i], 2), " - NEW BEST!", sep = "")
      }
      
      par(las=1)
      plot(silhvec, type="b", lwd=2, axes=F, cex=1.5, 
           xlab = "N. of Clusters", ylab = "Avg. Width",
           mar = c(2,2,0,0),
           ylim=c(0.9*min(silhvec),1.1*max(silhvec)))
      axis(1, 1:length(krange),krange)
      abline(v = krange-1, lty = 2, col = "gray", lwd = 1)
#       abline(h = 0.5, lty = 2, col = 2)
      points(krange[which.max(silhvec)]-(svalue(min_k_ran)-1), max(silhvec), pch=19, col=2, cex=1.3)      
      S3=Sys.time()
      cat(" - Required time: ",round(S3-S2, 3), " ",  attr(S3-S2, "units"), sep ="")
      svalue(sb) <- paste("Required time: ", S3-S2, attr(S3-S2, "units"), " with k =", i, sep =" ")
    }
    S4 = Sys.time()
    cat("\n\n   ---   Metier Discovery Completed!   ---\n\nBest Avg Width ", lb_CLA$data$silinfo$avg.width, " - with ", nrow(lb_CLA$data$medoids), " clusters in ", round(S4-S1, 3), " ", attr(S4-S1, "units"), "\n", sep = "")
    svalue(sb) <- paste("END\nCompleted in ", round(S4-S1, 3), " ", attr(S4-S1, "units"), sep = "")
    
    gconfirm("Metier discovery completed", title = "Confirm")
    enabled(save_bu) <- TRUE
    enabled(g_sup) <- TRUE
    }else{
      cat("\nNot enough data to perform Metier Discovery!", sep = "")
      enabled(g_sup) <- TRUE
    }
  })
  addSpring(g_go)
  addSpring(g_sup)
  enabled(g_sup) <- FALSE
  
  g_med <- gframe(text = "Results", horizontal = FALSE, container = big_g)
  theplot <- ggraphics(width = 700, height = 350, container = g_med, expand = TRUE)
  
  g_inf <- ggroup(horizontal = TRUE, container = big_g)
  addSpring(g_inf)
  save_bu <- gbutton(text = "Save\nMetier Discovery\nResults", container = g_inf, handler = function(h,...)
  {
    res_file <- gfile(text = "Save Metier Discovery result",
                      type = "save")
    saveRDS(lb_CLA, file = paste(res_file, ".rData", sep = ""))
    
    gconfirm("Metier discovery results saved", title = "Confirm")
  })
  enabled(save_bu) <- FALSE
  addSpring(g_inf)
  
  sb <- gstatusbar("", container = met_dis_win)
  
  enabled(start_b) <- FALSE
  
  if(lb_DB$db != "")
  {
    svalue(sel_lb_f) <- ifelse(.Platform$OS.type == "windows", strsplit(lb_DB$db, "\\\\")[[1]][length(strsplit(lb_DB$db, "\\\\")[[1]])],strsplit(lb_DB$db, "/")[[1]][length(strsplit(lb_DB$db, "/")[[1]])])
    enabled(start_b) <- TRUE
    enabled(g_sup) <- TRUE
  }
  
  visible(met_dis_win) <- TRUE
  
}