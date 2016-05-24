load_GWAS <-
function(filename, dir = getwd(), column_separators = c("\t", " ", "", ",", ";"), test_nrows = 1000,
                      header = TRUE, nrows = -1, comment.char = "", na.strings = c("NA", "."), stringsAsFactors = FALSE, ...) {
  if(!is.character(filename) | !is.character(dir)) {
    stop("argument 'filename' or 'dir' is not a character string") }
  
  FL <- load_test(filename, dir, column_separators = column_separators, test_nrows = test_nrows,
                  header = header, comment.char = comment.char, na.strings = na.strings, stringsAsFactors = stringsAsFactors, ... )
  if(!FL$success) {	stop(paste("unable to load data -", FL$error) ) }
  if(FL$type == "zip" | FL$type == ".gz") {
    if(FL$type == "zip") {
      output <- read.table(unz(paste(dir, filename, sep = "/"), substr(filename, 1L, nchar(filename) - 4L)), sep = FL$sep,
                           header = header, nrows = nrows, comment.char = comment.char, na.strings = na.strings, stringsAsFactors = stringsAsFactors, ... )
      close(unz(paste(dir, filename, sep = "/"), substr(filename, 1L, nchar(filename) - 4L)))
    } else {
      output <- read.table(gzfile(paste(dir, filename, sep = "/")), sep = FL$sep,
                           header = header, nrows = nrows, comment.char = comment.char, na.strings = na.strings, stringsAsFactors = stringsAsFactors, ... )
      close(gzfile(paste(dir, filename, sep = "/")))
    }
    return(output)
  } else {
    return(read.table(paste(dir, filename, sep = "/"), sep = FL$sep,
                         header = header, nrows = nrows, comment.char = comment.char, na.strings = na.strings, stringsAsFactors = stringsAsFactors, ... ))
  }
}
