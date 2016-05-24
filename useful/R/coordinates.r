#' @title pol2cart
#' @description Converts polar coordinates to caretsian coordinates
#' @details Converts polar coordinates to caretsian coordinates using a simple conversion.  The angle, \code{theta} must be in radians.
#' 
#' Somewhat inspired by http://www.r-bloggers.com/convert-polar-coordinates-to-cartesian/ and https://www.mathsisfun.com/polar-cartesian-coordinates.html
#' @export pol2cart
#' @importFrom dplyr data_frame
#' @aliases pol2cart
#' @author Jared P. Lander
#' @param r The radius of the point
#' @param theta The angle of the point, in radians
#' @param degrees Logical indicating if theta is specified in degrees
#' @return A data.frame holding the (x,y) coordinates and original polar coordinates
#' @examples 
#' 
#' polarRadPosTop <- data.frame(r=c(3, 5, 3, 5, 4, 6, 4, 6, 2), 
#'      theta=c(0, pi/6, pi/4, pi/3, pi/2, 2*pi/3, 3*pi/4, 5*pi/6, pi))
#' polarRadPosBottom <- data.frame(r=c(3, 5, 3, 5, 4, 6, 4, 6, 2), 
#'      theta=c(pi, 7*pi/6, 5*pi/4, 4*pi/3, 3*pi/2, 5*pi/3, 7*pi/4, 9*pi/6, 2*pi))
#' polarRadNegTop <- data.frame(r=c(3, 5, 3, 5, 4, 6, 4, 6, 2), 
#'      theta=-1*c(0, pi/6, pi/4, pi/3, pi/2, 2*pi/3, 3*pi/4, 5*pi/6, pi))
#' polarRadNegBottom <- data.frame(r=c(3, 5, 3, 5, 4, 6, 4, 6, 2), 
#'      theta=-1*c(pi, 7*pi/6, 5*pi/4, 4*pi/3, 3*pi/2, 5*pi/3, 7*pi/4, 9*pi/6, 2*pi))
#' 
#' pol2cart(polarRadPosTop$r, polarRadPosTop$theta)
#' pol2cart(polarRadPosBottom$r, polarRadPosBottom$theta)
#' pol2cart(polarRadNegTop$r, polarRadNegTop$theta)
#' pol2cart(polarRadNegBottom$r, polarRadNegBottom$theta)
#' 
pol2cart <- function(r, theta, degrees=FALSE)
{
    # convert degrees to raidans if so requested
    origTheta <- theta
    if(degrees)
    {
        theta <- theta*pi/180
    }
    
    # compute x
    x <- r*cos(theta)
    # compute y
    y <- r*sin(theta)
    
    data_frame(x=x, y=y, r=r, theta=origTheta)
}


#' @title cart2pol
#' @description Converts polar coordinates to caretsian coordinates
#' @details Converts polar coordinates to caretsian coordinates using a simple conversion.  The angle, \code{theta} must be in radians.
#' 
#' Somewhat inspired by http://www.r-bloggers.com/convert-polar-coordinates-to-cartesian/ and https://www.mathsisfun.com/polar-cartesian-coordinates.html
#' @export cart2pol
#' @importFrom magrittr "%<>%"
#' @importFrom dplyr mutate data_frame
#' @aliases cart2pol
#' @author Jared P. Lander
#' @param x The x-coordinate of the point
#' @param y The y-coordinate of the point
#' @param degrees Logical indicating if theta should be returned in degrees
#' @return A data.frame holding the polar coordinates and the original (x,y) coordinates
#' @examples 
#' 
#' library(dplyr)
#' x1 <- c(1, sqrt(3)/2, sqrt(2)/2, 1/2, 0)
#' y1 <- c(0, 1/2, sqrt(2)/2, sqrt(3)/2, 1)
#' d1 <- data_frame(x=x1, y=y1, Q='I')
#' 
#' x2 <- c(0, -1/2, -sqrt(2)/2, -sqrt(3)/2, -1)
#' y2 <- c(1, sqrt(3)/2, sqrt(2)/2, 1/2, 0)
#' d2 <- data_frame(x=x2, y=y2, Q='II')
#' 
#' x3 <- c(-1, -sqrt(3)/2, -sqrt(2)/2, -1/2, 0)
#' y3 <- c(0, -1/2, -sqrt(2)/2, -sqrt(3)/2, -1)
#' d3 <- data_frame(x=x3, y=y3, Q='III')
#' 
#' x4 <- c(0, 1/2, sqrt(2)/2, sqrt(3)/2, 1)
#' y4 <- c(-1, -sqrt(3)/2, -sqrt(2)/2, -1/2, 0)
#' d4 <- data_frame(x=x4, y=y4, Q='IV')
#' 
#' dAll <- bind_rows(d1, d2, d3, d4)
#' 
#' cart2pol(dAll$x, dAll$y)
#' cart2pol(dAll$x, dAll$y, degrees=TRUE)
#' 
cart2pol <- function(x, y, degrees=FALSE)
{
    # calculate r with sqrt of x and y squared
    r <- sqrt(x^2 + y^2)
    # calculate theta with arctan
    theta <- atan2(y, x)
    
    result <- data_frame(r=r, theta=theta, x=x, y=y)
    
    ## adjust angle for appropriate quadrant
    ## quadrants I and II need no adjustment
    ## quadrant III and IV, add 360
    result %<>% mutate(theta=theta + (y < 0)*2*pi)
    
    # return as degrees if requested
    if(degrees)
    {
        result %<>% mutate(theta=theta*180/pi)
    }
    
    return(result)
}