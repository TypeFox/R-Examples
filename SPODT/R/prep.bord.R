prep.bord <-
function(data)
{
    minX <- min(data$x) - (max(data$x) - min(data$x))/50
    maxX <- max(data$x) + (max(data$x) - min(data$x))/50
    minY <- min(data$y) - (max(data$y) - min(data$y))/50
    maxY <- max(data$y) + (max(data$y) - min(data$y))/50

    bord <- data.frame(c(minX, minX, maxX, maxX), c(minY, maxY, maxY, minY))
    colnames(bord) <- c("X", "Y")
    
    return(bord)
}
