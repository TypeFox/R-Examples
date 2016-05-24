#' @title Order IJDATA spot sequences and growth lines
#' @description Reorders spot sequences and growth lines within object of class \code{\link[=read.ijdata]{IJDATA}}.
#' @details Reorders IJDATA \code{spot.x} and \code{spot.y} and/or \code{gbs.x} and \code{gbs.y} coordinate data.frames. Useful when order of ROIs does not correspond with the desired order of \code{\link{convert.ijdata}} or spot.dist output. Can also be used to print the order of spot sequences and growth lines within IJDATA object (see 'print.order'). In addition the function can also be used to drop spot sequences or growth lines from the data set by leaving out ROI names. In this case a warning is produced to confirm that the user has not forgotten something. 
#' @param IJDATA an \code{\link[=read.ijdata]{IJDATA}} object.
#' @param spots a character or numeric vector specifying the desired order of sample spot sequences. 
#' @param gbs a character or numeric vector specifying the desired order of growth lines. 
#' @param print.order logical. Should the current order of spot sequences and growth lines be printed instead of changing the order? 
#' @author Mikko Vihtakari
#' @seealso \code{\link{read.ijdata}} for reading zip files containing ImageJ ROIs.
#' 
#' \code{\link{convert.ijdata}} for converting the coordinate information to \link[spatstat]{spatstat} point patterns. 
#' 
#' @examples data(shellspots)
#' order.ijdata(shellspots, print.order = TRUE) # Prints the current order. Does not change anything
#' dat <- order.ijdata(shellspots, gbs = c(1,3,6:14,4,5,2)) # Changes order of growth bands
#' order.ijdata(dat, print.order = TRUE)
#' 
#' ## Subset the first sample spot sequence
#' dat2 <- order.ijdata(shellspots, gbs = 1:13)
#' ## Warning message:
#' ## In order.ijdata(shellspots, gbs = 1:13) :
#' ## Length of gbs does not correspond the number of columns. Data removed.
#' order.ijdata(dat2, print.order = TRUE)
#' @export

order.ijdata <- function(IJDATA, spots = "", gbs = "", print.order = FALSE){
  
## Debugging parameters, remove when ready
#IJDATA <- dat; spots = NULL; spots = c("Laser", "SIMS1", "SIMS2", "SIMS3"); gbs = ""; print.order = T 

## Print order

if(print.order) order.list <- list(spots = coln(IJDATA$spots.x), gbs = coln(IJDATA$gbs.x))

spots.ncol <- ncol(IJDATA$spots.x)
gbs.ncol <- ncol(IJDATA$gbs.x)

### Spots
empty <- all(length(spots) == 1 & spots == "")
## If spots is empty, do not change the order
if(empty) {
  } else {
## If spots is a character vector, order by the character vector
if(class(spots) == "character") {
  IJDATA$spots.x <- IJDATA$spots.x[spots]
  IJDATA$spots.y <- IJDATA$spots.y[spots]
}
## If spots is numeric, order by the numeric values
if(class(spots) == "integer" | class(spots) == "numeric") {
  IJDATA$spots.x <- IJDATA$spots.x[spots]
  IJDATA$spots.y <- IJDATA$spots.y[spots]
}}

if(all(!empty, length(spots) != spots.ncol)) warning("Length of spots does not correspond the number of columns. Data removed.") 

### Gbs
empty <- all(length(gbs) == 1 & gbs == "")
## If gbs is empty, do not change the order
if(empty) {
  } else {
## If gbs is a character vector, order by the character vector
if(class(gbs) == "character") {
  IJDATA$gbs.x <- IJDATA$gbs.x[gbs]
  IJDATA$gbs.y <- IJDATA$gbs.y[gbs]
}
## If gbs is numeric, order by the numeric values
if(class(gbs) == "integer" | class(gbs) == "numeric") {
  IJDATA$gbs.x <- IJDATA$gbs.x[gbs]
  IJDATA$gbs.y <- IJDATA$gbs.y[gbs]
}}

if(all(!empty, length(gbs) != gbs.ncol)) warning("Length of gbs does not correspond the number of columns. Data removed.") 

if(print.order) return(order.list)
if(!print.order) return(IJDATA)}
