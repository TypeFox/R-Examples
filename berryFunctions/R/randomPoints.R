#' Distanced random points
#' 
#' Arranges points in square randomly, but with certain minimal distance to each other
#' 
#' @return data.frame with x and y coordinates.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, 2011/2012
#' @seealso \code{\link{distance}}
#' @keywords datagen spatial
#' @export
#' @examples
#' 
#' P <- randomPoints(xmin=200,xmax=700, ymin=300,ymax=680, number=60,mindist=10, asp=1)
#' rect(xleft=200, ybottom=300, xright=700, ytop=680, col=NA, border=1)
#' 
#' format( round(P,4), trim=FALSE)
#' 
#' \dontrun{dev.new(record=TRUE)}
#' ## Rcmd check --as-cran doesn't like to open external devices,
#' ## so this example is excluded from running in the checks.
#' for(i in 1:10)
#' {
#' rp <- randomPoints(xmin=0,xmax=20, ymin=0,ymax=20, number=20, mindist=3, plot=FALSE)
#' plot(rp, las=1, asp=1, pch=16)
#' abline(h=0:30*2, v=0:30*2, col=8); box()
#' for(i in 1:nrow(rp))
#'     circle(rp$x[i],rp$y[i], r=3, col=rgb(1,0,0,alpha=0.2), border=NA)
#' }
#' 
#' @param xmin Minimum x coordinate
#' @param xmax Upper limit x values
#' @param ymin Ditto for y
#' @param ymax And yet again: Ditto.
#' @param number How many points should be randomly + uniformly distributed
#' @param mindist Minimum DIstance each point should have to others
#' @param plot Plot the result? DEFAULT: TRUE
#' @param \dots Further arguments passed to plot
#' 
randomPoints <- function(
xmin,
xmax,
ymin,
ymax,
number,
mindist,
plot=TRUE,
...)
{
# benotigte Abstandsfunktion definieren:
distance <- function(xpt,ypt, xref,yref) sqrt((xref-xpt)^2 + (yref-ypt)^2)
# Zielvektoren fuer zufaellig verteilte Punkte erstellen
x <- y <- rep(NA, number)
# Ersten Wert reinschreiben
x[1] <- runif(1, xmin, xmax) ; y[1] <- runif(1, ymin, ymax)
# number-1  Punkte hinzufuegen
for(i in 2:number)
   {x[i] <- runif(1, xmin, xmax) ; y[i] <- runif(1, ymin, ymax)
   # Wenn minimale Distanz nicht gehalten, Punkt ersetzen
   while(   min(distance(x[i], y[i], x[-i],y[-i]), na.rm=TRUE) < mindist  )
   {x[i] <- runif(1, xmin, xmax) ; y[i] <- runif(1, ymin, ymax) }
   }
# wenn gewollt, Ergebnis plotten
if(plot) plot(x,y,las=1,pch=16, ...)
# Ergebnis ausgeben
return(data.frame(x,y))
}
