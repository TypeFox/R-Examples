# Calculate path offset and its max

.packageName <- 'mousetrack'

pathoffset <- function(x,y){

    # require("pracma") ## to compute the vector cross-product
    startend = c( (y[length(y)] - y[1]), (x[length(x)] - x[1]))

    startenddistance = sqrt(sum(startend^2))
    perpdistance = vector()

    for (m in 1:length(x)){
        pointstart = c(y[m] - y[1], x[m] - x[1])
        perpdistance = c(perpdistance,
            sqrt( sum( cross( c(startend, 0), c(pointstart,0)
                             )
                      )^2
                 )/ sqrt(sum(startend^2))
            )
    }
    
    pathoffset = perpdistance/startenddistance
    maxpathoffset = max(pathoffset)

    return(maxpathoffset)

}
