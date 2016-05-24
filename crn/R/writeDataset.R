writeDataset <- function(filename, cnames = colnamesDaily, varname = "T_DAILY_MEAN" ){
  HOURLY <- FALSE
  if ("LST_TIME" %in% cnames) HOURLY <- TRUE
  if (!(varname %in% cnames)){
    cat(varname, " is not in ", cnames, "\n")
    stop("check variable name")
    
  }
  
  
  
  cClass <- rep("NULL", length(cnames))
  cClass[1] <- "numeric"
  dex <- which(cnames == varname)
  cClass[dex] <- "numeric"
  
  
  if (HOURLY == TRUE){
    dirname        <- "HourlyData"
    cClass[c(3,4)] <- "numeric"
    
  } else {
    dirname   <- "DailyData"
    cClass[2] <- "numeric"
  }
  
  X <- read.table(filename, col.names = cnames, colClasses = cClass)
  if (!file.exists(dirname)) dir.create(dirname)   
  fname <- paste(varname, as.character(Sys.Date()), ".dat", sep = "")
  fname <- file.path(dirname, fname, fsep = .Platform$file.sep)          
  write.table(X, file = fname, row.names = FALSE)
  
}