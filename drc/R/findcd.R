"findcd" <- function(x, y, scaleInc = 0.001)
{ 
    yRange <- range(y)
    lenyRange <- scaleInc * diff(yRange)
#    cVal <- yRange[1] - lenyRange  # the c parameter        
#    dVal <- yRange[2] + lenyRange  # the d parameter
    
    c(yRange[1] - lenyRange, yRange[2] + lenyRange)
}