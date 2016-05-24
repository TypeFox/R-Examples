#### Function imagevect as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Aug. 18, 2011

imagevect <- function (x, labels, contour = FALSE, gridsize = 20, 
                        axes = TRUE, nlabx = 5, nlaby = 5, ...) 
{
    sort.x <- x[order(formatXY(labels))]
    rrr <- dimension(labels, unique = TRUE, sort = TRUE)
    dims <- c(length(rrr[[2]]), length(rrr[[1]]))
    dim(sort.x) <- dims
    par(xaxs = "i", yaxs = "i")
    image.plot(nnn <- t(sort.x), axes = FALSE, ...)

    if (contour) {
        contour(nnn, add = TRUE, ...)
    }
    
    if (axes) {
        points(0, 0, pch = " ", cex = 3)                      ## invoke the large plot
        
     get.axis.ticks <- 
              function(nlabs = NULL, gridsize = NULL, limit_max = NULL){
     
         ngrid <- (limit_max-0)/gridsize                      ## Obtain number of labels to plot
         per_grid <- 1/(ngrid-1)                              ## Obtain length for each grid 
         start <- 0 - (1/(ngrid-1))/2                         ## starting point for the ticks
         stop <- 1 + (1/(ngrid-1))/2                          ## stopping point for the ticks
         
         lab <-(0:nlabs*(limit_max/nlabs))
         at <- seq(from = start, to = stop, 
                   by = ((1 + per_grid))/((length(lab)-1)))   ## Position of the ticks
         
         return(list(lab, at))
         
     }
       xaxis.position <- get.axis.ticks(nlabs = nlabx, gridsize = gridsize, 
                                limit_max = gridsize * nrow(nnn))
       yaxis.position <- get.axis.ticks(nlabs = nlaby, gridsize = gridsize, 
                                limit_max = gridsize * ncol(nnn))
        
        axis(1, labels = xaxis.position[[1]], at = xaxis.position[[2]])
        axis(2, labels = yaxis.position[[1]], at = yaxis.position[[2]])
    }
    invisible(nnn)
}
