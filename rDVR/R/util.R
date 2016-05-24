#' Check for Video Server binary
#' 
#' \code{checkForVServer}
#' A utility function to check if the Video Server standalone binary is present.
#' @param jarloc A directory in which the binary is to be placed. Defaults to the /bin of rDVR package.
#' @param update A boolean indicating whether to update the binary if it is present.
#' @export
#' @section Detail: The Video Server java binary can be found at https://github.com/johndharrison/rDVR. If users would like to create their own please refe to documentation. This convience function downloads the standalone server and places it in the rDVR package directory bin folder by default.
#' @examples
#' \dontrun{
#' checkForVServer()
#' }

checkForVServer <- function (jarloc = NULL, update = FALSE) 
{
  vsURL <- "http://dl.dropboxusercontent.com/u/38391057/"
  jarLoc <- ifelse(is.null(jarloc), file.path(find.package("rDVR"), "bin"), jarloc)
  dvrFILE <- file.path(jarLoc, "videorecorderservice-2.0.jar")
  if (update || !file.exists(dvrFILE)) {
    cat("        PLEASE NOTE THIS FUNCTION WILL DOWNLOAD A STANDALONE BINARY JAVA
        JAR FROM https://dl.dropboxusercontent.com/u/38391057/videorecorderservice-2.0.jar. THIS JAR
        HAS BEEN COMPILED BY THE AUTHOR OF THIS PACKAGE> IF YOU WOULD 
        PREFER TO COMPILE YOUR OWN PLEASE REFER TO THE DOCUMENTATION.\n")
    ans <- readline("PLEASE REPLY (Y) yes TO PROCEED:")
    if(!identical(ans, "Y")){stop("Please agree to download or read documentation on rolling your own binary.")}
    dir.create(jarLoc, showWarnings=FALSE)
    print("DOWNLOADING STANDALONE VIDEO SERVER. THIS MAY TAKE SEVERAL MINUTES")
    download.file(paste0( vsURL, "videorecorderservice-2.0.jar"), dvrFILE, mode = "wb")
  }else{
    stop("FILE ALREADY EXISTS.")
  }
}

#' Start the Video Server.
#' 
#' \code{startVideoServer}
#' A utility function to start the standalone video server. 
#' @param jarloc A directory in which the standalone video server binary is located. Defaults to the /bin of rDVR package.
#' @param savedir A directory where the user would like videos saved to.  If not declared it defaults to the temp folder (which varies depending on the OS).
#' @param port The port on which the video server will listen. Defaults to 9998.
#' @param distmode You can enable a "distribution" mode for the storage of recorded videos that will use the last two characters of the filename requested to save video to place it in a subfolder. By default, if you want to save a video with the name, say, videofile20987 you will end up with the file stored at: /path/to/dest/folder/videofile20987.mov.
#' If you had enabled the distribution mode with distmode = TRUE the video would be stored at: /path/to/dest/folder/87/videofile20987.mov
#' @param invisible Windows Only Show the video server in a shell. By default it is invisible with setting TRUE. 
#' @export
#' @section Detail: By default the binary is assumed to be in
#' the rDVR package /bin directory. 
#' @examples
#' \dontrun{
#' startVideoServer()
#' }

startVideoServer <- function (jarloc = NULL, savedir = NULL, port = NULL, distmode = FALSE, invisible = TRUE) 
{
  jarLoc <- ifelse(is.null(jarloc), file.path(find.package("rDVR"), "bin"), jarloc)
  saveDIR <- savedir
  if(distmode){distMode <- "true"}else{distMode <- NULL}
  dvrFILE <- file.path(jarLoc, "videorecorderservice-2.0.jar")
  if (!file.exists(dvrFILE)) {
    stop("No Video Recorder binary exists. Run checkForVServer or start video server manually.")
  }
  else {
    jarOpt <- list('-DdestFolder=' = saveDIR
                   , '-Dport=' = port, '-DdistributeFiles=' = distMode, '-jar ' = shQuote(dvrFILE))
    jarOpt <- jarOpt[!sapply(jarOpt, is.null)]
    if (.Platform$OS.type == "unix") {
#      system(paste0('nohup java ', paste0(names(jarOpt), jarOpt, collapse = ' '), ">>/dev/null 2>>/dev/null &")
#             ,ignore.stdout = TRUE, ignore.stderr = TRUE, wait = FALSE)
      system(paste0('java ', paste0(names(jarOpt), jarOpt, collapse = ' '))
             ,ignore.stdout = TRUE, ignore.stderr = TRUE, wait = FALSE)
    }
    else {
      system(paste0('java ', paste0(names(jarOpt), jarOpt, collapse = ' ')), wait = FALSE, 
             invisible = invisible)
    }
  }
}

#' @export .DollarNames.rDVR
#' 

.DollarNames.rDVR <- function(x, pattern){
  grep(pattern, getRefClass(class(x))$methods(), value=TRUE)
}


