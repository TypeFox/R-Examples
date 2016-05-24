rename_column <- function(x, old_name, new_name, required){
  is_present <- old_name %in% colnames(x)
  if (required & !is_present){
    stop(paste("Expected column", old_name, "not found!"))
  }
  if (is_present){
    colnames(x)[colnames(x) == old_name] <- new_name
  }
  x
}

rename_cond_column <- function(x, old_name1, old_name2, new_name1, new_name2){
  is_present1 <- old_name1 %in% colnames(x)
  is_present2 <- old_name2 %in% colnames(x)
  if (!(is_present1 | is_present2)){
    stop(paste("Expected column", paste(old_name1, 'or', old_name2), "not found!"))
  }
  if (is_present1 | is_present2){
    try(colnames(x)[colnames(x) == old_name1] <- new_name1)
    try(colnames(x)[colnames(x) == old_name2] <- new_name2)
  }
  if (is_present1 & is_present2){
    x$CLOCK <- NULL
  }
  x
}

myPremodeling <- function(myData){
  myData$fNTAFD <- as.factor(myData$NTAFD)
  myData$ID <- as.factor(myData$ID)
  myData$RESPONSE <- as.numeric(myData$RESPONSE)
  myData$EXPOSURE <- as.numeric(myData$EXPOSURE)
  
  if("CLOCK" %in% colnames(myData)){
    # Parse clock time
    myData$CLOCK <- parse_date_time(myData$TOD, "hm") 
    # Transform clock time to 0-24 hours interval
    start.time <- parse_date_time("00:00", "hm")
    tmp <- myData$TOD-start.time
    tmp <- as.duration(tmp)
    myTOD <- as.numeric(tmp/3600)
    myData$TOD <- myTOD
    myData$CLOCK <- NULL
  }
  
  if("PERIOD" %in% colnames(myData)){
    myData$PERIOD <- as.factor(myData$PERIOD)
  }
  return(myData)
}
