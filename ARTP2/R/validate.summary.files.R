
validate.summary.files <- function(summary.files){
  
  # validate summary.files
  if(!is.vector(summary.files)){
    msg <- "Invalid summary.files"
    stop(msg)
    
    tmp <- !file.exists(summary.files)
    if(any(tmp)){
      msg <- paste(c("File(s) below were not found: ", summary.files[tmp]), collapse = "\n")
      stop(msg)
    }
  }
  
}
