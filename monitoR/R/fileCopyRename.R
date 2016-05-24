# Modified: 2014 MAR 29

fileCopyRename <-
function(
   files,                     # Vector or list of files to transfer, or specify 'from'
   from='.',                  # Entire SD card or other directory with field recordings, or specify 'files'
   to,                   # Destination directory
   csv.dir=to,                # Directory where table should be saved
   csv.name,               # Name of csv file
   loc.prefix,                # Location for files
   ext,                       # Extension
   rec.tz=NA,                 # Time zone setting for recorders
   hours.offset=0,            # optional offset to modification time
   CardRecorderID=NA,         # ID of SD card
   kaleidoscope=TRUE,         # If using Wildlife Acoustics' Kaleidoscope for converting WAC to WAV
   split.channels=FALSE,      # If using Kaleidoscope's split channel option only the left channel will be used
   metadata.only=FALSE,       # If moving new WAV files from 'Recordings' to 'Surveys' set to TRUE
   full.survey.names=FALSE,   # TRUE will use full file paths as survey names
   rename=TRUE, 	              # Set to FALSE to move files without renaming
   copy=TRUE				  # Set to FALSE to rename files without moving
   ){
   
   if(missing(csv.name) & length(loc.prefix) == 1) csv.name <- paste(loc.prefix, '_', as.character(Sys.time(), format='%Y-%m-%d'), '.csv', sep='')
   else csv.name <- paste('metadata_', as.character(Sys.time(), format='%Y-%m-%d'), '.csv', sep='')
   
   if (all(nchar(loc.prefix) != 6)) stop(paste('loc.prefix must be 6 characters, got', paste0(loc.prefix, collapse=", ")))

   # Create "to" directory if it doesn't exist
   if(!missing(to) && !file.exists(to)) dir.create(to)
   # Determine whether to use full survey names
   if(full.survey.names) fpath <- paste0(to, '/')
   else fpath <- ""

   if(metadata.only){
         if(missing(files)){ # if moving a whole directory
            files <- list.files(from, full.names=TRUE, pattern=paste0(".", ext))
            survey.names <- list.files(from, pattern=paste0(".", ext))
         } else { # if only moving a few files, not a whole directory
            survey.full.names <- strsplit(as.character(files), "/")
            survey.names <- unlist(lapply(survey.full.names, function(x) x[length(x)]))
         }
   
         # Info from within files
         if(ext == 'wav') {
            headers <- lapply(files, function(x) tuneR::readWave(x, header=TRUE)) 
            duration <- sapply(headers, function(x) round(x$samples/x$sample.rate, 2))
            samp.rate <- sapply(headers, function(x) x$sample.rate)
            bits <- sapply(headers, function(x) x$bits)
            channels <- sapply(headers, function(x) x$channels)
         } else if (ext == 'mp3') {
            waves <- lapply(files, readMP3) 
            duration <- sapply(waves, function(x) round(length(x@left)/x@samp.rate, 2))
            samp.rate <- sapply(waves, function(x) x@samp.rate)
            bits <- sapply(waves, function(x) x@bits)
            channels <- sapply(waves, function(x) x@stereo + 1)
         } else if (ext == 'wac') {
            duration <- samp.rate <- bits <- channels <- NA
         } else stop('Expect ext of wav, wac, or mp3, but got ', ext)
             
         # Put together table
         survey.dat <- data.frame(fldSurveyName=survey.names, fldSurveyLength=duration, fldSampleRate=samp.rate, fldBitsperSample=bits, fldChannels=channels)
         # Create a destination path for file.copy
         newfiles <- paste(to, '/', survey.names, sep="")
   } else { # Full functionality
         # Create vectors of data that will become data frame columns
         if(missing(files)){ # if moving a whole directory:
            files <- list.files(from, full.names=TRUE, pattern=paste0("\\.", ext))
            survey.names <- list.files(from, pattern=paste0("\\.", ext))
         } else { # if only moving a few files:
            survey.full.names <- strsplit(as.character(files), "/")
            survey.names <- unlist(lapply(survey.full.names, function(x) x[length(x)]))
         }
         # Manually correct for timezone changes
         mdates <- file.info(files)$mtime+hours.offset*3600
         if(is.na(rec.tz)) rec.tz <- format(mdates, format='%Z')
         mdates <- format(mdates, tz=rec.tz, format='%Y-%m-%d_%H%M%S_%Z')
         newfiles <- paste(to, '/', loc.prefix, '_', mdates, '.', ext, sep="")
         if(!kaleidoscope && rename) {newfilenames <- paste0(fpath, loc.prefix, '_', mdates, '.', ext)
         } else if(kaleidoscope && !split.channels && rename) {newfilenames <- paste0(fpath, loc.prefix, '_', mdates, '_00000_000.wav')
         } else if(kaleidoscope && split.channels && rename) {newfilenames <- paste0(fpath, loc.prefix, '_', mdates, '_1_00000.wav')
         } else if(!rename) {newfilenames <- files}
      
         # Info from within files
         if(ext == 'wav') {
            headers <- lapply(files, function(x) tuneR::readWave(x, header=TRUE)) 
            duration <- sapply(headers, function(x) round(x$samples/x$sample.rate, 2))
            samp.rate <- sapply(headers, function(x) x$sample.rate)
            bits <- sapply(headers, function(x) x$bits)
            channels <- sapply(headers, function(x) x$channels)
         } else if (ext == 'mp3') {
            waves <- lapply(files, readMP3) 
            duration <- sapply(waves, function(x) round(length(x@left)/x@samp.rate, 2))
            samp.rate <- sapply(waves, function(x) x@samp.rate)
            bits <- sapply(waves, function(x) x@bits)
            channels <- sapply(waves, function(x) x@stereo + 1)
         } else if (ext == 'wac') {
            duration <- samp.rate <- bits <- channels <- NA
         } else stop('Expect ext of wav, wac, or mp3, but got ', ext)
             
         # Put together table
         survey.dat <- data.frame(fldOriginalDateModified=mdates, fldOriginalRecordingName=survey.names, fldSurveyName=newfilenames, fldRecordingFormat=ext, fkCardRecorderID=CardRecorderID, fldSurveyLength=duration, fldSampleRate=samp.rate, fldBitsperSample=bits, fldChannels=channels)
         }
   # Move files
   if(copy) file.copy(files, newfiles)

   # Return results
   write.csv(x=survey.dat, file=paste(csv.dir, '/', csv.name, sep=''), quote=FALSE, row.names=FALSE)
   invisible(survey.dat)
}
