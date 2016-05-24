# Function to read only portions of MP3 files, using mp3splt to locate frames
# Can lead to a lot of disk writes if used to loop through long spectrogram files!
# Adapted from tuneR::readMP3 by J. Katz 02 Feb 2013
# 2014 MAR 29

readMP3 <- function (
   filename,                       # Name of full-length MP3 file
   from,                       # Start time to cut from, seconds, decimals rounded to the hundredths
   to                         # End time to cut to, seconds, decimals rounded to the hundredths
   ) { 
      # First four if's straight from tuneR::readMP3 
      if (!is.character(filename)) 
         stop("'filename' must be of type character.")
      if (length(filename) != 1) 
         stop("Please specify exactly one 'filename'.")
      if (!file.exists(filename)) 
         stop("File '", filename, "' does not exist.")
      if (file.access(filename, 4)) 
         stop("No read permission for file ", filename)
      if (!missing(from) && !missing(to)) { # Begin 'cut' option
          # Use a terminal/cmd to see if mp3splt is installed
         if (length(system("mp3splt -h", intern=TRUE))<10) {
            cat("'mp3splt' not found; use tuneR::readMP3 equivalent? (Y/n)\n")
			x <- tolower(readLines(n=1))
			if(x == 'y' || x == ""){
               wave <- tuneR::readMP3(filename)
               return(wave)
            } else stop('Download mp3splt at\nhttp://mp3splt.sourceforge.net/mp3splt_page/home.php')
         }
         if(from>=to) # Make sure cut points are logical
            stop("Invalid cut points: from>=to.")
         if(from>=60) { # Format seconds to mm.ss.xx, where xx are hundredths 
            min <- from%/%60
            sec <- floor(from%%60)+round(from%%floor(from), 2)            
            if(nchar(sec)<5 && round(from%%floor(from), 2)>0) sec <- paste0('0', sec)
            from <- paste0(min, '.', sec)
         } else if(from<10) {
            from <- paste0('00.0', round(from, 2))
         } else from <- paste0('00.', round(from, 2))
         if (to>=60) {
            min <- to%/%60
            sec <- floor(to%%60)+round(to%%floor(to), 2)
            if(nchar(sec)<5 && round(to%%floor(to), 2)>0) sec <- paste0('0', sec)
            to <- paste0(min, '.', sec)
         } else if (to<10) {
            to <- paste0('00.0', round(to, 2))
         } else to <- paste0('00.', round(to, 2))
         # Call mp3splt to cut up the mp3 file and make a temporary exerpt
#         filenameQ <- paste0("'", filename, "'")
         system(paste('mp3splt -Q -d', getwd(), filename, from, to, '-o temp16jp86mz2'))
         # Use tuneR::readMP3 process on the temporary file
         return(tuneR::readMP3('temp16jp86mz2.mp3'))
      } else if(any(all(missing(from) && !missing(to)), all(!missing(from) && missing(to)))) { 
         stop('Single cut point provided, must provide both from and to.')
      } else {
         # Continue standard tuneR::readMP3 process w/out 'cut' option
         return(tuneR::readMP3(filename=filename))
      }
   }
   
