fread.OTU <- function(file, sep="auto") {
  
  input <- data.table::fread(file, na.strings=c("",".","NA"))
  df <- as.data.frame(input)
  rownames(df) <- df[,1]
  df <- df[, -1]
  valid.OTU(df)
  
  df$taxonomy <- as.character(df$taxonomy)
  
  return(df)
}

read.OTU <- function(file, sep=",") {
  
  df <- read.table(file, sep=sep, header=T, row.names=1,
                   na.strings=c("",".","NA"))
 
  valid.OTU(df)
  
  df$taxonomy <- as.character(df$taxonomy)
  
  return(df)
}

fread.meta <- function(file, sep="auto") {
  
  input <- data.table::fread(file, na.strings=c("",".","NA"))
  df <- as.data.frame(input)
  rownames(df) <- df[,1]
  df <- df[, -1]
 
  return(df)
}

read.meta <- function(file, sep=",") {
  
  input <- read.table(file, header=TRUE, row.names=1, sep=sep,
                      na.strings=c("",".","NA"))
  ### PUT META DATA VALIDATION HERE
    
  return(input)
}

write.data <- function(data, file) {
  
  file <- .ensure.filepath(file, "csv")
  
  write.csv(data, file=file, quote=FALSE)
}


