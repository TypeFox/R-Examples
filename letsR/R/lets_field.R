#' Create species' values based on the species co-occurrence within focal ranges
#' 
#' @author Bruno Vilela & Fabricio Villalobos
#' 
#' @description Create single species' values based on the attributes of species 
#' co-occurring within individual ranges.
#' 
#' @param x A \code{\link{PresenceAbsence}} object or a presence-absence in \code{matrix}
#' format (see xy argument for matrix use) with the species named in the columns.
#' @param y Species attribute to be considered. It must be a numeric attribute.
#' @param z Species names in the same order as the attributes and exactly the 
#' same as named in the \code{matrix} or in the \code{PresenceAbsence} object.
#' @param weight If \code{TRUE} the value is weighted by species' range size, 
#' if \code{FALSE} the value is the mean of all species that co-occur within the
#'  focal species.
#' @param xy If \code{TRUE} the presence-absence \code{matrix} contains 
#' the cell coordinates in the first two columns.
#' @param count Logical, if \code{TRUE} a counting window will open.
#'  
#' @details If the species do not co-occur with any other species NaN will be 
#' returned. 
#'
#' @references Villalobos, F. and Arita, H.T. 2010. The diversity field of 
#' New World leaf-nosed bats (Phyllostomidae). Global Ecology and Biogeography. 
#' 19, 200-211. 
#' @references Villalobos, F., Rangel, T.F., and Diniz-Filho, J.A.F. 2013. 
#' Phylogenetic fields of species: cross-species patterns of phylogenetic 
#' structure and geographical coexistence. Proceedings of the Royal 
#' Society B. 280, 20122570.
#' 
#' @seealso \code{\link{lets.presab.birds}}
#' @seealso \code{\link{lets.presab}}
#' 
#' @examples \dontrun{
#' data(PAM)
#' range <- lets.rangesize(x = PAM, units = "cell")
#' field <- lets.field(PAM, range, PAM$S, weight = TRUE)
#' }
#' 
#' @export


lets.field <- function(x, y, z, weight = TRUE, 
                       xy = NULL, count = FALSE) {
  
  # Get the matrix without the coordinates
  if (class(x) == "PresenceAbsence") {
    p <- x[[1]][, -(1:2)]
    namesSpe <- x[[3]]
  } 
  
  if (is.matrix(x)) {
    if (is.null(xy)) {
      stop(paste("Please set if your matrix contains ",
                 "coordinates in the first two columns",
                 "(xy argument)."), sep = "")
    }
    if (xy) {
      x <- x[, -(1:2)]
    }
    p <- x
    namesSpe <- colnames(x)
  }
  
  # Change to numeric if factor
  if (is.factor(y)) {
    y <- as.numeric(levels(y))[y]
  }
  
  # Save it in another object
  p2 <- p
  
  
  for(i in 1:ncol(p2)) {
    pos <- z == namesSpe[i]
    if (length(pos) > 0) {
      p2[, i] <- p2[, i] * y[pos]
      pos2 <- p2[, i] == 0
      p2[pos2, i] <- NA
    } else {
      p2[, i] <- NA
    }
  }
  
  media <- numeric(ncol(p))
  n <- length(media)
  
  # With count 
  if (count) {
    if (!"tools:rstudio" %in% search()) {
      dev.new(width = 2, height = 2, pointsize = 12)
      par(mar = c(0, 0, 0, 0))
    }
    for(i in 1:n){
      plot.new()
      text(0.5, 0.5, paste(paste("Total:", n, "\n", "Runs to go: ", (n - i))))
      media[i] <- .InsiLoop(i, p, p2, weight)
    }
    dev.off()
  }
  
  # Without count
  if (!count) {
    for(i in 1:n){          
      media[i] <- .InsiLoop(i, p, p2, weight)
    }
  }
  
  # Return process
  resultado <- cbind(namesSpe, media)
  colnames(resultado) <- c("Species", "Value")
  resultado <- as.data.frame(resultado)
  # Changing the factors
  resultado[, 2] <- as.numeric(levels(resultado[, 2]))[resultado[, 2]]
  return(resultado)
}



# Axuliar function to avoid code repetition inside the loop <<<<<<<<<

.InsiLoop <- function(i, p, p2, weight) {
  pos3 <- p[, i] == 1
  p3 <- p[pos3, -i, drop = FALSE]
  p4 <- p2[pos3, -i, drop = FALSE]
  mult <- p3 * p4
  if (weight) {
    mediai <- mean(mult, na.rm = TRUE)
  }  else {
    mult <- matrix(mult, ncol = (ncol(p) - 1))
    me <- colMeans(mult, na.rm = TRUE)
    mediai <- mean(me, na.rm = TRUE)
  }
  return(mediai)
}