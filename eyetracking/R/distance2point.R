distance2point <-
function( x, y, viewerDistance, viewerHeight, resolutionX, resolutionY, screenWidth, screenHeight ) {
    
    centerX <- screenWidth / 2
    centerY <- screenHeight / 2 - viewerHeight
    
    targetX <- x / resolutionX * screenWidth
    targetY <- y / resolutionY * screenHeight
    
    dX <- targetX - centerX
    dY <- targetY - centerY
    
    screenDistance <- sqrt( dX^2 + dY^2 )
    
    sqrt( ( viewerDistance^2 + screenDistance^2) )
    
}

