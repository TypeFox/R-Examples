#' Manually add/detect points to a scatter plot figure.
#'
#' Allows for the user to manually add an unlimited number of points to a 
#' figure image, by left-clicking over a figure's point.  Click on the red  
#' upper-right box called "EXIT" to end recording the position of manually 
#' detected points.
#'
#' @param file The file name and location of a scatter plot figure.  Prompts
#'    for file name if none is explicitly called.  Can also be a binary figure 
#'    image with detected points (an \code{EBImage} object).  
#'    See: \code{\link{figure_detectAllPoints}}
#' @param color The color to paint the manually detected points; default is 
#'    green.
#' @param size The radius of the painted points.
#'
#' @return A data frame with detected points.
#' 
#' @seealso \link{figure_detectAllPoints}
#' 
#' @import grid
#' @importFrom EBImage readImage
#' @importFrom grDevices rgb col2rgb
#' @export

figure_add <- function (file = file.choose(),
                        color = "#009900",
                        size = 0.03) {

  # set up grid window
  grid.newpage()
  vp <- viewport(xscale = c(0, 1), yscale = c(0, 1))
  pushViewport(vp)

  # plot figure image on vp
  grid.raster(ifelse(class(file) == "Image", file, readImage(file)))

  # plot exit box on figure
  grid.rect(x = 0.05, 
            y = 0.95, 
            width = 0.1, 
            height = 0.1, 
            gp = gpar(fill = "red", col = "red", alpha = 0.3))
  grid.text(x = 0.05, 
            y = 0.95, 
            "EXIT", 
            gp = gpar(fontface = "bold"))

  # extract click coordinates and plot click-marker circles
  newPoints <- 0
  
  # capture click events and keep track of [x,y] location
  repeat{
  
    #get click x and y location
    click.loc <- grid.locator()

    if(newPoints == 0) {
      theCoordinates <- matrix(c(unitTrim(click.loc$x), unitTrim(click.loc$y)), 
                               ncol = 2, 
                               byrow = TRUE)
    } else {
      theCoordinates <- rbind(theCoordinates, 
                              matrix(c(unitTrim(click.loc$x), unitTrim(click.loc$y)), 
                                     ncol = 2,
                                     byrow = TRUE))
    }
    # keep track of clicked points
    newPoints <- newPoints + 1  

    # exit manual detection when user clicks within EXIT box region
    if(unitTrim(click.loc$x) < unitTrim(unit(0.1,"native")) & 
         unitTrim(click.loc$y) > unitTrim(unit(0.9,"native"))) {
      # delete last 'exit' click
      theCoordinates <- theCoordinates[-nrow(theCoordinates), ]
      break;
    }
    
    # paint the figure with a marker circle at the clicked location
    theColor <- rgb(t(col2rgb(color)), maxColorValue = 255)
    grid.circle(x = click.loc$x, 
                y = click.loc$y, 
                r = size, 
                gp=gpar(col = theColor,
                        fill = theColor,
                        alpha = 0.3))
    grid.circle(x = click.loc$x, 
                y = click.loc$y, 
                r = 0.003, 
                gp=gpar(col = theColor,
                        fill = theColor))
    
  }

  # organize the coordinates of clicked points in a data frame
  theCoordinates <- data.frame(theCoordinates)
  colnames(theCoordinates) <- c("x","y")
  theCoordinates["cluster"] <- "manual"
                                          
  return(theCoordinates)
}
