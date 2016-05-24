# Functions for sorting out sound clips
# It can accept (1) Wave objects or file paths to (2) wav or (3) mp3 objects, and return either 1 or 2
# Modified: 6 Sept 2015

getClip <- function(
  clip, 
  name="clip", 
  output="file", 
  write.wav=FALSE
) {

  if(class(clip) == "list" | (class(clip) == "character" && length(clip)>1)) {
    clist <- list()
    if(grepl(", ",name)) {
      name <- gsub(".*\\(", "", name)
      name <- gsub("\\)", "", name)
      name <- strsplit(name, ",")[[1]]
    }
    for(i in seq(length(clip))) {
      clist[[i]] <- getOneClip(clip[[i]], paste0(name[i], i), output, write.wav) 
    }
    return(clist)
  } 

  return(getOneClip(clip, name, output, write.wav))

}

getOneClip <- function(
  clip, 
  name, 
  output, 
  write.wav
) {

  if(output == "file") {
    if(class(clip) == "Wave") {
      fname <- paste0(name, ".wav")
      if(!write.wav) {
	stop("output argument is \"file\" but write.wav argument is FALSE so this function will not create a file. Set write.wav=TRUE to create a file, or else specify a wav file instead of a Wave object.")
      }
      if(file.exists(fname)) stop("Will not create a wav file from this clip because a file with name ", fname, " already exists.")
      else tuneR::writeWave(clip, fname) 
      return(fname)
    } else 
    if(class(clip) == "character") {
      if(!file.exists(clip)) stop("clip argument seems to be a file name but no file with the name ", clip, " exists.")
      return(clip)
    } else 
    stop("Can\'t figure out what to do with this clip:", clip, "with class:", class(clip))
  }

  if(output == "Wave") {
    if(class(clip) == "Wave") {
      return(clip)
    } else 
    if(class(clip) == "character") {
      if(!file.exists(clip)) stop("clip argument seems to be a file name but no file with the name ", clip, " exists!")
      file.ext <- tolower(gsub(".*\\.", "", clip))
      if(file.ext == "wav") 
        clip <- tuneR::readWave(clip) else 
      if(file.ext == "mp3") 
        clip <- readMP3(clip) else stop("File extension must be wav or mp3, but got ", file.ext)
      return(clip)
    } else 
    stop("Can\'t figure out what to do with this clip:", clip, "with class:", class(clip))
  }

}   

# Reads a single wav or mp3 file 
readClip <- function(clip) {

  if(class(clip) != "character" | length(clip) != 1) stop("Expected a length-one character vector for clip, but got a length ", length(clip), " ", class(clip), " object.")
  if(!file.exists(clip)) stop("clip argument seems to be a file name but no file with the name ", clip, " exists!")
 
  file.ext <- tolower(gsub(".*\\.", "", clip))
  if(file.ext == "wav") return(tuneR::readWave(filename=clip))
  if(file.ext == "mp3") return(readMP3(filename=clip)) 
  stop("File extension must be wav or mp3, but got ", file.ext)

}

