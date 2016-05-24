#' Internal function for making molaR plot legends
#' 
#' Function properly scales legends and prints them to the background of rgl
#' devices
#' @param expression it knows what to do...
#'
#' 
#' molaR_bgplot()

molaR_bgplot <- function(expression){
    filename <- tempfile(fileext = ".png")
    png(filename = filename, width = 800, height = 800)
    value <- expression
    dev.off()
    result <- bg3d(texture = filename, col = "white", lit = FALSE)
    invisible(structure(result, value = value))
}