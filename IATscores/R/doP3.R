doP3 <- function(IATdata, 
                 P3 = c("dscore", "gscore", "wpr90", "minid", "minid_t10",
                        "minid_w10", "minid_i10"), maxMemory = 1000,
                 verbose = TRUE)
{
  Scores <- data.frame("subject" = unique(IATdata$subject))
  # convert data to long format
  id <- c("subject", "latency", "correct", "blockcode",
          "index", "praccrit")
  
  IATlong <- melt(IATdata, id = id,
                  value.name = "RT")
  
  # P3 = dscore
  if("dscore" %in% P3)
  {
    if(verbose) print(paste0(Sys.time(), ": Applying parameter P3 = dscore"))
    Dsc <- Dscores(IATlong)
    Scores <- left_join(Scores, Dsc, by = "subject")
    rm(list = "Dsc") 
  }
  
  # P3 = gscore
  if("gscore" %in% P3)
  {
    if(verbose) print(paste0(Sys.time(), ": Applying parameter P3 = gscore"))
    Gsc <- Gscores(IATlong)
    Scores <- left_join(Scores, Gsc, by = "subject")
    rm(list = "Gsc") 
  }
  
  # P3 = WPR
  if("wpr90" %in% P3)
  {
    if(verbose) print(paste0(Sys.time(), ": Applying parameter P3 = wpr"))
    WPRsc <- WPRscores(IATlong)
    Scores <- left_join(Scores, WPRsc, by = "subject")
    rm(list = "WPRsc")
  }
  
  # P3 = minid
  if("minid" %in% P3)
  {
    if(verbose) print(paste0(Sys.time(), ": Applying parameter P3 = minid and variants"))
    MDsc <- MiniDscores(IATlong, P3, maxMemory = maxMemory)
    Scores <- left_join(Scores, MDsc, by = "subject")
    rm(list = "MDsc")      
  }
  Scores
}