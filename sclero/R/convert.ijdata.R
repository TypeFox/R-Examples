#' @title Convert IJDATA object to a list of spatstat objects
#' @description Converts an \code{\link[=read.ijdata]{IJDATA}} to a list of \link[spatstat]{spatstat} patterns ready for \link[=plot.rawDist]{plotting} or \link[=spot.dist]{sample spot alignment}. 
#' @param X an \code{\link[=read.ijdata]{IJDATA}} object to be converted.
#' @return Returns a list of class \code{rawDist}, which contains \link[spatstat]{spatstat} point patterns. The returned \code{rawDist} can be plotted using the generic plotting command.
#' @author Mikko Vihtakari
#' @seealso \code{\link{read.ijdata}} for generating \code{IJDATA} objects.
#' 
#' \code{\link{plot.rawDist}} for plotting.
#' 
#' \code{\link{spot.dist}} for aligning sample spots.
#' 
#' @examples data(shellspots)
#' shell_map <- convert.ijdata(shellspots)
#' plot(shell_map)
#' @import spatstat plyr
#' @export

convert.ijdata <- function(X){
  
tmp <- lapply(c("spots.x", "gbs.x", "main.x"), function(i) range(X[[i]], na.rm = TRUE))
  
range.x <- c(round_any(do.call("min", tmp), 10, floor), round_any(do.call("max", tmp), 10, ceiling))

tmp <- lapply(c("spots.y", "gbs.y", "main.y"), function(i) range(X[[i]], na.rm = TRUE))

range.y <- c(round_any(do.call("min", tmp), 10, floor), round_any(do.call("max", tmp), 10, ceiling))

## Define observation window

owin <- owin(xrange = range.x, yrange =  range.y, unitname = X$unit)

## Make spots to spatstat ppp objects

spots <- list()
for(i in 1:ncol(X$spots.x)){
x <- na.omit(X$spots.x[,i])
y <- na.omit(X$spots.y[,i])
spots <- c(spots, list(ppp(x, y, marks = as.factor(seq(1, length(x))), window = owin)))}
names(spots) <- colnames(X$spots.x)

# Growth lines

x <- X$gbs.x
y <- X$gbs.y
nlines <- ncol(x)
lines <- colnames(x)

for(j in 1:nlines){
xx <- x[,j][!is.na(x[,j])]
yy <- y[,j][!is.na(y[,j])]
assign(lines[j], psp(x0 = xx[1:(length(xx)-1)], x1 = xx[2:length(xx)], y0 = yy[1:(length(yy)-1)], y1 = yy[2:length(yy)], marks =  rep(lines[j], (length(xx)-1)), window = owin))}

a <- get(lines[1])
for(j in 2:nlines){
a <- superimpose(a, get(lines[j]))}
gbs <- a

# Main axis

x <- X$main.x
x <- x[!is.na(x)]
y <- X$main.y
y <- y[!is.na(y)]

main <- psp(x0 = x[1], x1 = x[2], y0 = y[1], y1 = y[2], marks = c("main"), window = owin)

start <- ppp(x[1], y[1], marks = c("start"), window = owin)
end <- ppp(x[2], y[2], marks = c("end"), window = owin)

df <- list(spots = spots, gbs = gbs, main = main, window = owin, start.main = start, end.main = end, sample.name = X$sample.name, scaling.factor = X$scaling.factor, unit = X$unit)
  
class(df) <- "rawDist"
  
return(df)}

