load_test <-
function(filename, dir = getwd(), column_separators = c("\t", " ", "", ",", ";"), test_nrows = 1000, ... ) {
  file_dir <- paste(dir, filename, sep = "/")
  if(!file.exists(file_dir)) {
    print(" - - ERROR: file not found!")
    return(list(success = FALSE, error = "file not found"))
  }
  load_success <- FALSE
  load_error	<- "no errors"
  file_type <- tolower(substr(filename, nchar(filename) - 2L, nchar(filename)))
  if(file_type == "zip" | file_type == ".gz") {
    file_cal <- if(file_type == "zip") parse(text = "unz(file_dir, substr(filename, 1L, nchar(filename) - 4L))") else parse(text = "gzfile(file_dir)")
  } else {	file_cal <- parse(text = "file_dir") }
  
  for(countI in 1:length(column_separators)) {
    tryI <- try(ncol(read.table(eval(file_cal), sep = column_separators[countI], nrows = test_nrows, ...) ), silent = TRUE)
    if(class(tryI) == "try-error") { load_error <- tryI
    } else {
      if(tryI < 5L) { load_error <- "insufficient columns"
      } else {
        load_error <- NA
        load_success <- TRUE
        break
      }
    }
  }
  if(file_type == "zip" | file_type == ".gz") close(eval(file_cal))
  return(list(success = load_success, error = load_error, type = file_type, sep = column_separators[countI]))	
}
