## compute the area subtending the trajectory x,y
## the vectors of the coordinates, x,y points

.packageName <- 'mousetrack'

areaunder <- function(x, y){
# require("caTools")

area1 = trapz(abs(x),abs(y));

if (area1 > .5){
    farea = area1 -.5
} else { farea = 0;
     }

    return(farea)

}
