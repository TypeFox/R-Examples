#' Create a matrix summarizing species' attributes within 
#' cells of a PresenceAbsence object
#' 
#' @author Bruno Vilela
#' 
#' @description Summarize species atributes per cell in a presence-absence matrix.
#' 
#' @param x A \code{\link{PresenceAbsence}} object.
#' @param y Species attribute to be considered. It must be a numeric attribute.
#' @param z Species names in the same order as the attributes and exactly the 
#' same as named in the \code{PresenceAbsence} object.
#' @param func A function to summarize the species' atribute in each cell (the function must return only one value).
#' @param ras If \code{TRUE} the raster object will be returned 
#' together with the matrix.
#' 
#' @return The result can be both a \code{matrix} or a \code{list} cointaining 
#' the follow objects:
#' @return \strong{Matrix}: a \code{matrix} object with the cells' geographic 
#' coordinates and the summarized species' attributes within them.
#' @return \strong{Raster}: The summarized species'attributed maped in a 
#' \code{raster} object.
#' 
#' @seealso \code{\link{lets.presab}}
#' @seealso \code{\link{lets.presab.birds}}
#' 
#' 
#' @examples \dontrun{
#' data(PAM)
#' data(IUCN)
#' trait <- IUCN$Description_Year
#' resu <- lets.maplizer(PAM, trait, PAM$S, ras = TRUE)
#' head(resu$Matrix)
#' plot(resu$Raster, xlab = "Longitude", ylab = "Latitude", 
#' main = "Mean description year per site") ; map(add = TRUE)
#' }
#' 
#' @export

lets.maplizer <- function(x, y, z, func = mean, ras = FALSE) {
  
  # Change factor to numbers
  if (is.factor(y)) {
    y <- as.numeric(levels(y))[y]
  }
  
  # To avoid being transformed in NA
  y[y == 0] <- 0.00000000000000000000000000000000000001
  
  # Get the matrix without coordinates
  p <- x[[1]][, -(1:2)]
  
  for(i in 1:ncol(p)) {
    pos <- x[[3]][i] == z
    if (length(pos) > 0) {
      p[, i] <- p[, i] * y[pos]
      pos2 <- p[, i] == 0
      p[pos2, i] <- NA
    } else {
      p[, i] <- NA
    }
  }
  
  func2 <- function(x) {
    pos <- is.na(x) 
    resu <- func(x[!pos])
  }
  
  resum <- apply(p, 1, func2)
  
  # Matrix of result 
  resultado <- cbind(x[[1]][, 1:2], resum)
  resu2 <- na.omit(resultado)
  
  # Name change
  name <- paste("Variable", as.character(substitute(func)),
                sep = "_")
  colnames(resultado)[3] <- name 
  
  # Back to zero
  resultado[resultado <= 0.0000000000000000000000000000000000001] <- 0
  
  # Return result with or without the raster
  if (ras) {
    r <- rasterize(resu2[, 1:2], x[[2]], resu2[, 3])
    return(list(Matrix = resultado, Raster = r))
  } else {
    return(resultado)
  }
}
