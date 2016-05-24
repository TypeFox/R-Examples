displayGraphic <- function (data, colourScale = 1, add = FALSE, windowLimits = NULL, asp = NULL)
{
    dims <- dim(data)
    scale <- getColourScale(colourScale)
    
    if (!is.null(windowLimits))
    {
        data[data < min(windowLimits)] <- min(windowLimits)
        data[data > max(windowLimits)] <- max(windowLimits)
    }
    
    if (add)
    {
        data <- replace(data, which(data==0), NA)
        if (!all(is.na(data)))
        {
            if (is.null(windowLimits))
                image(data, col=scale$colours, add=TRUE, useRaster=TRUE)
            else
                image(data, col=scale$colours, add=TRUE, useRaster=TRUE, zlim=sort(windowLimits))
        }
    }
    else
    {
        if (is.null(asp))
            asp <- dims[2] / dims[1]
        
        oldPars <- par(mai=c(0,0,0,0), bg=scale$background)
        if (is.null(windowLimits))
            image(data, col=scale$colours, useRaster=TRUE, axes=FALSE, asp=asp)
        else
            image(data, col=scale$colours, useRaster=TRUE, axes=FALSE, asp=asp, zlim=sort(windowLimits))
        par(oldPars)
    }
}
