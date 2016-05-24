subtendedAngle <-
function( x1, y1, x2, y2, viewerDistance=58.74, viewerHeight=4.55, resolutionX=1280, resolutionY=1024, screenWidth=33.97, screenHeight=27.31 ) {
    
    d1 <- distance2point(x1, y1, viewerDistance, viewerHeight, resolutionX, resolutionY, screenWidth, screenHeight)
    d2 <- distance2point(x2, y2, viewerDistance, viewerHeight, resolutionX, resolutionY, screenWidth, screenHeight)
    
    dX <- screenWidth * ( ( x2 - x1 ) / resolutionX )
    dY <- screenWidth * ( ( y2 - y1 ) / resolutionY )
    
    screenDistance <- sqrt( dX^2 + dY^2 )
    
    angleRadians <- acos( ( d1^2 + d2^2 - screenDistance^2 ) / ( 2 * d1 * d2 ) )
    
    ( angleRadians / ( 2 * pi ) ) * 360
    
}

