plot.NightDay <- function(x, maps = "world", add = FALSE, ...) {

   yy <- x$Latitude
   dec <- x$Declination
   GHA <- x$GHA
   Time <- x$Time
   tz <- x$tz

   x <- round(GHA)
   x0 <- 360

   if (x < 180){
     x <- x*(-1)
   }else{
     x <- x0 - x
   }

   if( tz >= 0 ) {
     tz <- paste("+", tz)
   }

   l <- list(dec, GHA)
   names(l) <- c("Dec", "GHA")

   if (!add) {
   map(maps)
   map.axes()
   title(main = paste( Time$year + 1900,"-", Time$mon + 1,"-", Time$mday, "        ",
                Time$hour, ":", Time$min,":", round(Time$sec), "\nGMT", tz, collapse = NULL),
                sub = paste("GHA: ", round(GHA, 2), "   Dec: ", round(dec ,2)))
   abline(h=c(-66.5, -23.5, 23.5, 66.5), col="blue", lty=3)
   abline(h=dec, col="red")
   abline(h=0, v=0, lty=2)
   lines(setdiff((x-length(yy)):(x+length(yy)),0), yy[c(length(yy):1, 1:length(yy))])
   points(x, dec, pch=21, cex=3, bg="yellow")
   }
   polygon(setdiff((x-length(yy)):(x+length(yy)),0), yy[c(length(yy):1, 1:length(yy))], col="grey", density=20, border=NA)
   if(dec < 0){
     rect(xleft=-200, ybottom=max(yy), xright=200, ytop=90, density=20, col="grey", border=NA)
     }
   if(dec > 0){
     rect(xleft=-200, ybottom=-90, xright=200, ytop=min(yy), density=20, col="grey", border=NA)
     }
}
