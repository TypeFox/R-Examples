#' Get the local copy of the spnet user manual
#' 
#' This function copies the spnet user manual to a user defined directory.
#' @param where the location where to copy the user manual. Default is the working directory.
#' @param overwrite logical; should existing destination files be overwritten?
#' @export
spnet.get.local.user.manual <- function(where = getwd(), overwrite = FALSE) {
  docpath <- system.file('doc', package = 'spnet')
  file.copy(
    from = file.path(docpath, list.files(docpath)[1]),
    to = where
  )
}


#' Convert colors to contrasted gray level for black and white rendering
#' 
#' This function converts color codes (given in hexadecimal format) to constrasted gray levels. This is useful to enhance readability when printing in black and white. The conversion from color to gray levels is performed using the luminosity method (0.21R + 0.72G + 0.07B), then levels are linearly scaled to [contrast.min;contrast.max].
#' @param x a \code{character}, the vector of color codes given in hexadecimal format (ex "#21AD5C").
#' @param increase.contrast a \code{logical}, if \code{TRUE} the scaling is performed.
#' @param contrast.min the minimal gray value to use in the scaling (0 = white, 1 = black).
#' @param contrast.max the maximal gray value to use in the scaling (0 = white, 1 = black).
#' @importFrom grDevices col2rgb
#' @importFrom grDevices gray
#' @export
#' @examples
#' mycols = c("#BA364E", "#32BAA1", "#12AA91")
#' color2blackwhite(mycols)
#' 
#' barplot(1:3, axes=FALSE, col=mycols)
#' barplot(1:3, axes=FALSE, col=color2blackwhite(mycols, increase.contrast = FALSE))
#' barplot(1:3, axes=FALSE, col=color2blackwhite(mycols))
#' barplot(1:3, axes=FALSE, col=color2blackwhite(mycols, contrast.min = 0, contrast.max = 1))
color2blackwhite <- function(x, increase.contrast = TRUE, contrast.min = 0.02, contrast.max = 0.98) { # hexadecimal
  
  # color encoding: 0 for black, 1 for white, this isn't intuitive for the user so we reverse the order
  contrast.min = 1 - contrast.min
  contrast.max = 1 - contrast.max
  # when reverting color we also need to reverse min and max parameters
  temp = contrast.max
  contrast.max = contrast.min
  contrast.min = temp  
  
  mat <- matrix(c(0.21, 0.72, 0.07), ncol = 1, dimnames=list(c('red', 'green', 'blue'), NULL))
  l <- length(x)
  out <- numeric(l)
  # converting color to gray equivalent using the luminosity method
  for (i in 1:l) {
    out[i] <- sum(col2rgb(x[i])*mat)
  }
  # increasing the contrast between the resulting gray values
  if(increase.contrast) {
    f <- function(x, min, max) {
      return((contrast.max-contrast.min)/(max-min)*(x - min) + contrast.min)
    }
    out <- mapply(f, out, rep(min(out),l),rep(max(out),l))*255
  }
  out <- round(out)
  out <- gray(out/255)
  return(out)
}