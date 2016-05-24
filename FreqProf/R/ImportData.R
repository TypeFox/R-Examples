#' Import Data Pop Up
#' 
#' This function reads a file, whose extension is either csv, bin or fpw, and 
#' imports it as a data.frame.
#' 
#' @param filename a string indicating the path of the file containing the data.
#'   By default, will open a pop-up so that the user can choose a file with the 
#'   GUI.
#' @return A data.frame ready to be converted into freqprof class (see function 
#'   \code{freqprof}).
#' @importFrom utils read.csv
#' @export
#' @examples
#'  \donttest{
#'  # Select a file
#'  import_data(filename = system.file("extdata", "S58-1-1.bin", 
#'                                     package = "FreqProf"))
#'  }
import_data = function(filename = file.choose()) {
  file.extension = tolower(substr(filename, 
                                  nchar(filename) - 2,
                                  nchar(filename)))
  
  data = switch(file.extension,
                csv = read.csv(filename),
                bin = read.bin(filename),
                fpw = read.fpw(filename))
  
  if(is.null(data)) {
    stop("file extension must be either csv, fpw, or bin")  
  }
  
  return(data)
}

#' Reads the data in the file "filename", which is supposed to be a .bin file
#' 
#' @param filename a string indicating the path of the file.
#' @return A data.frame giving the raw data.
#' @export
#' @examples
#'  \donttest{
#'  read.bin(filename = system.file("extdata", "S58-1-1.bin",
#'                                  package = "FreqProf"))
#'  }
read.bin = function(filename) {
  
  # scanning the file, line by line
  file.scan = scan(file  = filename,
                   what  = "character",
                   sep   = "\n",
                   quiet = TRUE)
  
  # looking for the line that starts with an asterisk
  N = which(sapply(1:length(file.scan),
                   function(x) {
                     substr(file.scan[x], 1, 1) == "*"
                   })
            )
  
  file.data = file.scan[(N + 2):length(file.scan)]
  
  # removing "<" at the end of line
  for(j in 1:length(file.data)) {
    file.data[j] = substr(file.data[j], 1, nchar(file.data[j]) - 1)
  }
  
  # converting these strings into a matrix
  result = matrix(data = NA, 
                  nrow = nchar(file.data[1]),
                  ncol = length(file.data))
  
  for(j in 1:dim(result)[2]) {
    result[, j] = sapply(1:nchar(file.data[j]),
                         function(x) {
                           as.numeric(substr(file.data[j], x, x))
                         })
  }
  
  return(as.data.frame(result))
}

#' Reads the data in the file "filename", which is supposed to be a .fpw file
#' 
#' @param filename a string indicating the path of the file.
#' @return A data.frame of the .fpw data.
#' @importFrom utils read.table
#' 
read.fpw = function(filename) {
           
  file.scan = scan(file  = filename,
                   what  = "character",
                   sep   = "\n",
                   quiet = TRUE)
  
  # looking for the line that indicates "[DATA]"
  N = which(sapply(1:length(file.scan),
                   function(x) {
                     file.scan[x] == "[DATA]"
                   })
            )
  
  # converting these strings into a data.frame
  result = read.table(file = filename, 
                      skip = N + 2, 
                      sep  = " ")
  
  # as the lines begin with a " ", the first column is read as NA
  result = result[, 2:ncol(result)]
  names(result) = paste("V", 1:ncol(result), sep = "")
  
  return(result)
}
