#' @title Read ImageJ zip file containing several ROI files and extract coordinate information.
#' 
#' @description A wrapper function, which reads an ImageJ zip file containing a collection of ROI files and outputs a list of data frames ready for \code{\link{convert.ijdata}} function.
#' 
#' @param X character string defining the name (including extension) or file path of an ImageJ zip file. Alternatively an \code{\link[=read.ijzip]{ijzip}} object.
#' @param spots optional. A character argument specifying the type of ROI objects that should be considered as sample spot sequences. Alternatively a numeric vector specifying the order of elements or a character vector specifying the names of ROI objects that should be assigned as sampling spot sequences. Defaults to \code{"point"} (See "Details" for further information).
#' @param gbs optional. A character argument, numeric vector or character vector specifying the type of ROI objects that should be considered as growth bands. Defauls to \code{"polyline"}. f left empty the remaining elements that are not defined in holes and main are assumed to be growth lines. For further information see \code{spots} and "Details". At the moment \code{spots} must be specified for this option to work.
#' @param main optional. A character argument specifying the type of ROI object that should be considered as the measurement axis. Only one measurement axis per ImageJ .zip file is allowed. Defaults to \code{"line"}.  At the moment \code{spots} must be specified for this option to work.
#' @param names optional. A character argument specifying how the names of \code{spots} and \code{gbs} should be generated. These names will be used in further functions (\code{\link{convert.ijdata}}, \code{\link{spot.dist}}). In general, it is adviced to use simple ROI names without special characters (for example \code{-} is not allowed in a ROI name; see 'Details'). Possible \code{names} options are:
#' \itemize{
#' \item \code{"generate.invalid"} (default). Uses the ROI object names, except when they are not valid \code{\link{data.frame}} column names. In the latter case sequential names will be generated.  
#' \item \code{"generate"}. Generates sequential names for all elements.
#' \item \code{"keep"}. Uses the ROI object names, except when they are not valid \code{\link{data.frame}} column names. In the latter case \code{\link{make.names}} function will be used to generate \code{data.frame} combatible column names.
#' \item \code{"force.keep"}. Uses the ROI object names as they are in the .zip file. Using this option might cause problems in consequent functions and is not recommended.
#' \item \code{"manual"}. Names for both \code{spots} and \code{gbs} are searched from \code{spot.names} and \code{gbs.names} arguments, respectively. 
#' \item \code{"manual.spots"}. Names for \code{spots} are searched from \code{spot.names} argument. Names for \code{gbs} are generated following \code{"generate.invalid"}.
#'  \item \code{"manual.gbs"}. Names for \code{gbs} are searched from \code{gbs.names} argument. Names for \code{spots} are generated following \code{"generate.invalid"}.
#' }
#' @param spot.names optional. A character vector of equal length to \code{spots} defining the names of sample spot sequences. Required if \code{names = "manual"} or \code{"manual.spots"}. Ignored otherwise.
#' @param gbs.names optional. A character vector of equal length to \code{gbs} defining the names of growth bands. Required if \code{names = "manual"} or \code{"manual.gbs"}. Ignored otherwise.
#' @param main.name optional. A character vector of lenght 1 defining the name of the measurement axis (\code{main}). If \code{main.name = "keep"}, the ROI object name will be used (not recommended, see "Details"). Otherwise the name will be taken from the argument. Defaults to \code{"main"}. 
#' @param sample.name optional. A character vector of length 1 defining the name of the sample. File name without the extension or alternatively \code{ijzip} object name is used as a default (\code{sample.name = "file"}).
#' @param scale optional. A numeric value defining the scale of photograph in pixels / \code{unit}. Defaults to 1.
#' @param unit optional. A charater vector of length 1 defining the unit of measurements. See \code{scale}.
#' 
#' @details In order to minimize the amount of text to be typed by a user, ROI objects of type "point" (this includes the "Multi-point Tool" points) are considered as sample spot sequences (\code{spots}) by default. Further, all "polyline" objects are assumed as growth bands (\code{gbs}) and "line" objects as the measurement axis (\code{main}) resulting to that only one "line" object is allowed per .zip file using the default settings. Alternatively, the user can specify the \code{spots}, \code{gbs}, and \code{main} objects manually using the order of the ImageJ .zip file with the exception that \bold{only one measurement axis is allowed} per \code{\link[=convert.ijdata]{rawDist}} or \code{\link[=spot.dist]{spotDist}} object. 
#' 
#'Punctuation characters other than \code{_} or \code{.} should not be used as names of \code{spots} or \code{gbs}, because they tend to confuse the internal \code{\link[base]{grep}} functions in \code{\link{spot.dist}} function. Hence it is adviced to use one of the options renaming invalid names of \code{spots} and \code{gbs} (\code{"generate.invalid"}, \code{"generate"}, \code{"keep"}). 
#' 
#' @return Returns an "IJDATA" object, which is a list of data frames containing the x and y coordinates for sampling spot sequences (\code{spots.x} and \code{spots.y}), growth bands (\code{gbs.x} and \code{gbs.y}), and measurement axis (\code{main.x} and \code{main.y}) together with sample name, scaling factor and unit of measurement.
#' @author Mikko Vihtakari 
#' @seealso \code{\link{order.ijdata}} for ordering and subsetting \code{read.ijdata} output.
#' 
#' \code{\link{convert.ijdata}} for converting the coordinate information to \link[spatstat]{spatstat} point patterns. 
#' 
#' \code{\link{spot.dist}} for aligning sample spot sequences. 
#' 
#' \code{\link[RImageJROI]{read.ijroi}} and \code{\link[RImageJROI]{read.ijzip}} for reading ImageJ ROI and .zip files.
#' 
#' @examples 
#' # Locate the example zip file
#' path <- file.path(system.file("extdata", package = "sclero"), "shellspots.zip") 
#' 
#' # You can replace 'path' by 'Your_file_name.zip'
#' dat <- read.ijdata(path) 
#' summary(dat)
#' 
#' ## Works also for IJZIP objects
#' dat2 <- read.ijzip(path)
#' dat2 <- read.ijdata(dat2)
#' dat[!(dat %in% dat2)] # Only the sample name differs
#' @import RImageJROI
#' @export

read.ijdata <- function(X, spots = 'point', gbs = 'polyline', main = 'line', names = 'generate.invalid', spot.names = NULL, gbs.names = NULL, main.name = "main", sample.name = 'file', scale = 1, unit = NULL){

## Debugging parameters, remove when ready
#X = file; spots = "point"; gbs = "polyline"; main = "line"; names = "generate.invalid"; spot.names = NULL; gbs.names = NULL; main.name = "main"; sample.name = "file"; scale = 1; unit = NULL

## Read zip file or load IJZIP object
if(class(X) == "character") {
  tmp <- grep(".zip", X, value = T)
  if(length(tmp) != 1) stop("Something wrong with the file name")
  if(length(tmp) == 1 & tmp==X){
  dat <- read.ijzip(X, names = TRUE)
    } 
  } else {
    if(class(X) == "ijzip") {
      dat <- X} else {
        stop("X is not path to a ImageJ .zip file nor a ijzip object")
      }}

## Define parameters  

types <- as.data.frame(do.call("rbind", lapply(dat, function(x) x$strType)))

## 1. Spots

if(class(spots) == "character" & length(spots) == 1 & any(spots %in% types[,1])) hole.seqs <- which(types[,1] %in% spots)
if(class(spots) == "character" & !any(spots %in% types[,1])) hole.seqs <- which(rownames(types) %in% spots)
if(class(spots) == "integer" | class(spots) == "numeric") hole.seqs <- spots

## Find coordinates for spot sequences, x-axis

tmp <-  lapply(hole.seqs, function(i) dat[[i]]$coords[,1]/scale) # Use scale argument to scale coordinates
n <- max(unlist(lapply(tmp, length)))
spots.x <- as.data.frame(do.call("cbind", lapply(tmp, function(x) {length(x) <- n; x})))

## X names
if(names == "force.keep"){
  colnames(spots.x) <- rownames(types)[hole.seqs]} else {
    if(names == "keep"){
      colnames(spots.x) <- make.names(row.names(types))[hole.seqs]} else { 
      if(names == "generate.invalid" | names == "manual.gbs"){
        change <- make.names(row.names(types))[hole.seqs] == row.names(types)[hole.seqs]
        colnames(spots.x) <- ifelse(change, row.names(types)[hole.seqs], paste0("s", 1:ncol(spots.x)))} else {
          if(names == "generate"){
            colnames(spots.x) <- paste0("s", 1:ncol(spots.x))} else {
              if(names == "manual" | names == "manual.spots") {
                if((ncol(spots.x) == length(spot.names)) == FALSE) stop("number of spots and spot.names differ")
                colnames(spots.x) <- spot.names} 
              }
            }
          }
        }

## Find coordinates for spot sequences, y-axis
  
tmp <-  lapply(hole.seqs, function(i) dat[[i]]$coords[,2]/scale) # Use scale argument to scale coordinates
n <- max(unlist(lapply(tmp, length)))
spots.y <- as.data.frame(do.call("cbind", lapply(tmp, function(x) {length(x) <- n; x})))

## Y names
if(names == "force.keep"){
  colnames(spots.y) <- rownames(types)[hole.seqs]} else {
    if(names == "keep"){
      colnames(spots.y) <- make.names(row.names(types))[hole.seqs]} else { 
      if(names == "generate.invalid" | names == "manual.gbs"){
        change <- make.names(row.names(types))[hole.seqs] == row.names(types)[hole.seqs]
        colnames(spots.y) <- ifelse(change, row.names(types)[hole.seqs], paste0("s", 1:ncol(spots.y)))} else {
          if(names == "generate"){
            colnames(spots.y) <- paste0("s", 1:ncol(spots.y))} else {
              if(names == "manual" | names == "manual.spots") {
                if((ncol(spots.y) == length(spot.names)) == FALSE) stop("number of spots and spot.names differ")
                colnames(spots.y) <- spot.names} 
              }
            }
          }
        }
  
## 2. Main axis

if(length(main) != 1) stop("Only one main axis allowed")
if(class(main) == "character" & any(main %in% types[,1])) main.l <- which(types[,1] %in% main)
if(class(main) == "character" & !any(main %in% types[,1])) main.l <- which(rownames(types) %in% main)
if(class(main) == "integer" | class(main) == "numeric") main.l <- main

## Find coordinates
  
main.x <- data.frame(main = dat[[main.l]]$coords[,1]/scale) # Use scale argument to scale coordinates
main.y <- data.frame(main = dat[[main.l]]$coords[,2]/scale)

## Names
  
if(length(main.name) != 1) {
  stop("The required length of main.name is 1. Use either 'keep' or a custom name.")} else {
    if(main.name != "keep") {
      colnames(main.x) <- main.name
      colnames(main.y) <- main.name
    }
  }  

## 3. Growth bands  
  
if(is.null(gbs)) gbss <- seq_along(row.names(types))[-c(hole.seqs, main.l)] else {
  if(class(gbs) == "character" & length(gbs) == 1 & any(gbs %in% types[,1])) {
    gbss <- which(types[,1] %in% gbs)} else {
      if(class(gbs) == "integer" | class(gbs) == "numeric") gbss <- gbs}}

## Find coordinates, x-axis  

tmp <-  lapply(gbss, function(i) dat[[i]]$coords[,1]/scale) # Use scale argument to scale coordinates
n <- max(unlist(lapply(tmp, length)))
gbs.x <- as.data.frame(do.call("cbind", lapply(tmp, function(x) {length(x) <- n; x})))

## X names  

if(names == "force.keep"){
  colnames(gbs.x) <- rownames(types)[gbss]} else {
    if(names == "keep"){
      colnames(gbs.x) <- make.names(row.names(types))[gbss]} else { 
      if(names == "generate.invalid" | names == "manual.spots"){
        change <- make.names(row.names(types))[gbss] == row.names(types)[gbss]
        colnames(gbs.x) <- ifelse(change, row.names(types)[gbss], paste0("l", 1:ncol(gbs.x)))} else {
          if(names == "generate"){
            colnames(gbs.x) <- paste0("l", 1:ncol(gbs.x))} else {
              if(names == "manual" | names == "manual.gbs") {
                if((ncol(gbs.x) == length(gbs.names)) == FALSE) stop("number of gbs and gbs.names differ")
                colnames(gbs.x) <- gbs.names} 
              }
            }
          }
        }

## Find coordinates, y-axis
  
tmp <-  lapply(gbss, function(i) dat[[i]]$coords[,2]/scale) # Use scale argument to scale coordinates
n <- max(unlist(lapply(tmp, length)))
gbs.y <- as.data.frame(do.call("cbind", lapply(tmp, function(x) {length(x) <- n; x})))

## Y names  

if(names == "force.keep"){
  colnames(gbs.y) <- rownames(types)[gbss]} else {
    if(names == "keep"){
      colnames(gbs.y) <- make.names(row.names(types))[gbss]} else { 
      if(names == "generate.invalid" | names == "manual.spots"){
        change <- make.names(row.names(types))[gbss] == row.names(types)[gbss]
        colnames(gbs.y) <- ifelse(change, row.names(types)[gbss], paste0("l", 1:ncol(gbs.y)))} else {
          if(names == "generate"){
            colnames(gbs.y) <- paste0("l", 1:ncol(gbs.y))} else {
              if(names == "manual" | names == "manual.gbs") {
                if((ncol(gbs.y) == length(gbs.names)) == FALSE) stop("number of gbs and gbs.names differ")
                colnames(gbs.y) <- gbs.names} 
              }
            }
          }
        }

## 4. Sample name 

if(length(sample.name) > 1) stop("Length of sample name > 1")

if(class(X) == "ijzip") {
  deparse(substitute(X))}

name <- if(sample.name == "file") {
  if(class(X) == "ijzip") {
  deparse(substitute(X))} else {
    if(grepl(".zip", X)) {
      if(grepl('/', X)) {
      tmp <- unlist(strsplit(X, '/'))
      sub(".zip", "", tmp[length(tmp)])} else {
     sub(".zip", "", X)}
  }}} else sample.name

## 5. Compile to a list

dat <- list(spots.x = spots.x, spots.y = spots.y, gbs.x = gbs.x, gbs.y = gbs.y, main.x = main.x, main.y = main.y, sample.name = name, scaling.factor = scale, unit = unit)

class(dat) <- "IJDATA"
  
return(dat)}

