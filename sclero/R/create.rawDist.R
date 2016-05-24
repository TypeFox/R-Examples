#' @title Create rawDist data from arbitary coordinates
#' 
#' @description Creates a \code{rawDist} data object from arbitary coordinates ready for \link[=plot.rawDist]{plotting} or \link[=spot.dist]{sample spot alignment}.  
#' 
#' @param spots A list of \code{\link[spatstat]{ppp}} objects or a single \code{ppp} object defining the sample spot sequences. If \code{\link[spatstat]{marks}} are not specified, sequencial names will be used.
#' @param gbs \code{\link[spatstat]{psp}} object defining the growth lines. If \code{\link[spatstat]{marks}} are not specified, sequencial names will be used.
#' @param main \code{\link[spatstat]{psp}} object defining the measurement axis. If \code{\link[spatstat]{marks}} are not specified, sequencial 'main' will assigned as the marks. Only on `main` axis is allowed per `rawDist` object.
#' @param spot.seq.names optional. A character vector of equal length to number of sample \code{spots} sequences defining the name for each sequence. If left empty sequencial names will be generated.
#' @param sample.name optional. A character vector (\code{length == 1}) defining the name of the sample.
#' @param scaling.factor optional. A numeric value defining the scale of photograph in pixels / \code{unit}. Defaults to 1.
#' @param unit optional. A charater vector (\code{length == 1}) defining the unit of measurements. See \code{scale}.
#' @details This function can be used to create arbitary test data, which can be passed further on in \code{sclero} package function hierarchy. 
#' @return Returns a list of class \code{rawDist}, which contains \link[spatstat]{spatstat} point patterns. The returned \code{rawDist} can be plotted using the generic plotting command.
#' @author Mikko Vihtakari
#' @seealso \code{\link{read.ijdata}} for generating \code{IJDATA} objects.
#' 
#' \code{\link{convert.ijdata}} for converting \code{IJDATA} objects to \code{rawDist} objects.
#' 
#' \code{\link{plot.rawDist}} for plotting.
#' 
#' \code{\link{spot.dist}} for aligning sample spots.
#' 
#' @examples 
#' dev.off()
#' W <- square(10)
#' S <- ppp(x = c(7, 5, 3), y = rep(5,3), window = W)
#' G <- psp(x0 = c(8,6,4,2), y0 = rep(2,4), x1 = c(8,6,4,2), y1 = rep(8,4), window = W)
#' M <- psp(x0 = 0, x1 = 8, y0 = 5, y1 = 5, window = W)
#' x <- create.rawDist(spots = S, gbs = G, main = M)
#' plot(x)
#' 
#' ## Generate random points for alignment
#' set.seed(1)
#' S <- rpoint(n = 5, win = owin(xrange = c(2,7), yrange = c(5,7)))
#' S$window <- W
#' G <- psp(x0 = c(7,5,3,1), y0 = rep(2,4), x1 = c(9,7,5,3), y1 = rep(8,4), window = W)
#' M <- psp(x0 = 0, x1 = 8, y0 = 1, y1 = 1, window = W)
#' x <- create.rawDist(spots = S, gbs = G, main = M)
#' plot(x)
#' y <- spot.dist(x)
#' plot(y)
#' @import spatstat
#' @export

create.rawDist <- function(spots, gbs, main, spot.seq.names = NULL, sample.name = NULL, scaling.factor = 1, unit = NULL){

if(!class(gbs)[1] %in% c("psp")) stop("'gbs' is not a spatstat 'psp' object")
if(!class(main)[1] %in% c("psp")) stop("'main' is not a  spatstat 'psp' object")

if(class(spots) == "list") {
  if(!all(unlist(lapply(spots, function(k) class(k) == "ppp")))) stop("list of 'spots' contains elements that are not spatstat 'ppp' objects") 
} else {
  if(!class(spots)[1] %in% "ppp") stop("'spots' is not a spatstat 'ppp' object")
}
  
if(main$n != 1) stop("Only one main axis allowed")
if(is.null(marks(main))) marks(main) <- "main"

if(is.null(marks(gbs))) marks(gbs) <- paste0("l", 1:gbs$n)
  
if(class(spots) == "list") {spots <- spots} else {
if(class(spots) == "ppp") spots <- list(spots)} 
 
if(is.null(spot.seq.names)) {
  names(spots) <- paste0("s", 1:length(spots))} else {
    if(length(spot.seq.names) == length(spots)) {names(spots) <- spot.seq.names} else
  stop("length of 'spot.seq.names' and number of 'spots' sequences differ")}

spots <- lapply(spots, function(k){
  if(is.null(marks(k))) marks(k) <- 1:k$n
  k})
  
range.x <- range(c(spots$window$xrange, gbs$window$xrange, main$window$xrange))
range.y <- range(c(spots$window$yrange, gbs$window$yrange, main$window$yrange))

Win <- owin(xrange = range.x, yrange = range.y)

start <- ppp(x = main$ends$x0, y = main$ends$y0, marks = c("start"), window = Win)
end <- ppp(x = main$ends$x1, y = main$ends$y1, marks = c("end"), window = Win)
  
df <- list(spots = spots, gbs = gbs, main = main, window = Win, start.main = start, end.main = end, sample.name = sample.name, scaling.factor =  scaling.factor, unit =  unit)
  
class(df) <- "rawDist"
  
return(df)}