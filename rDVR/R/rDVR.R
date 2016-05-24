#' CLASS rDVR
#'
#' class to communicate with the Video Server.
#'
#' This is class composing of methods to communicate with the Java Video Server.
#' The methods start a video, save a running video and stop a running video without saving. Please note videos are currently limited to 10 minutes in length. If you want to change this limit you will need to compile a custom binary.
#' Videos are encoded using the Apple QuickTime RLE codec. Currently to change the codec requires compiling a custom binary.
#'
#' @import RCurl
#' @import methods
#' @section Fields:
#' \describe{
#' \item{\code{remoteServerAddr}:}{Object of class \code{"character"}, giving the ip of the remote server. Defaults to localhost}
#' \item{\code{port}:}{Object of class \code{"numeric"}, the port of the remote server on which to connect to the Video Server.}
#' \item{\code{saveDir}:}{Object of class \code{"character"}, the location which the videos are saved. This is set by the Video Server. By default it is the OS temp directory. This option should be left NULL unless you set the directory when starting the Video Server.}
#' }
#'
#' @section Methods:
#' \describe{
#' \item{\code{new(...)}:}{ Create a new \code{rDVR} object. ... is used to define the appropriate slots.}
#' \item{\code{start(fileName, silent)}:}{ Start a new video recording. 
#' \describe{
#' \item{\code{fileName: }}{The filename by which to save your video. Defaults to Rtemp. This will be saved as RTemp.mov. YourNAME will be saved as YourNAME.mov.}
#' \item{\code{silent: }}{A boolean. If TRUE the method will run silently}
#' }
#' }
#' \item{\code{save(silent, replace)}:}{ Save a currently running video to file. The file is given by the fileName used in the \code{start} method. The directory that the file is written to is the save directory stipulated when the Video Server was started.
#' \describe{
#' \item{\code{silent: }}{A boolean. If TRUE the method will run silently}
#' \item{\code{replace: }}{A boolean. If TRUE then rDVR will replace a file if it already exists. If false rDVR will append 'copy' to the saveFile name.}
#' }
#' }
#' \item{\code{stop(silent)}:}{ Stops a currently running video. Using stop rather then save will result in the video being stopped but not saved.
#' \describe{
#' \item{\code{silent: }}{A boolean. If TRUE the method will run silently}
#' }
#' }
#' \item{\code{closeServer(silent)}:}{ Stops a currently running Video Server.
#' \describe{
#' \item{\code{silent: }}{A boolean. If TRUE the method will run silently}
#' }
#' }
#' }
#'
#' @export rDVR
#' @exportClass rDVR
#' @aliases rDVR
#' 

rDVR <- setRefClass("rDVR",
                    fields = list(
                      remoteServerAddr = "character",
                      port = "integer",
                      serverURL = "character",
                      saveDir = "character",
                      saveFile = "character"
                    ),
                    methods = list(
                      initialize = function(remoteServerAddr = "localhost",
                                            port = 9998L,
                                            saveDir = sub("(.*)/.*", "\\1", normalizePath(tempdir(), winslash='/'))
                      ){
                        remoteServerAddr <<- remoteServerAddr
                        port <<- port
                        saveDir <<- saveDir
                      },
                      
                      start = function(fileName = 'Rtemp', silent = FALSE){
                        if(!silent){print("Starting recording of video:")}
                        #if(file.exists(file.path(saveDir,'Rtemp')))
                        serverURL <<- paste0("http://",remoteServerAddr,":",port)
                        saveFile <<- fileName
                        ipAddr <- paste0(serverURL, '/rec/start')
                        res <- getURLContent(ipAddr, customrequest = "GET", isHTTP = FALSE)
                      },
                      
                      save = function(silent = FALSE, replace = TRUE){
                        if(!silent){print("Saving video ")}
                        if(file.exists(file.path(saveDir, paste0(saveFile, '.mov')))){
                          if(replace){
                            file.remove(file.path(saveDir, paste0(saveFile, '.mov')))
                          }else{
                            saveFile <<- paste0(saveFile, "copy")
                          }
                        }
                        ipAddr <- paste0(serverURL, '/rec/save/', saveFile)
                        res <- getURLContent(ipAddr, customrequest = "GET", isHTTP = FALSE)
                      },
                      
                      stop = function(silent = FALSE){
                        if(!silent){print("Stopped recording. Video not saved.")}
                        ipAddr <- paste0(serverURL, '/rec/stop')
                        res <- getURLContent(ipAddr, customrequest = "GET", isHTTP = FALSE)
                      },
                      
                      closeServer = function(silent = FALSE){
                        if(!silent){print("Close the Video Server.")}
                        ipAddr <- paste0(serverURL, '/rec/closeserver')
                        tryCatch(res <- getURLContent(ipAddr, customrequest = "GET", isHTTP = FALSE), error = function(e){})
                      }
                    )
)

