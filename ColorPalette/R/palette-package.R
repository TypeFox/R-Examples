#' hsv2rgb convert
#'
#' This function converts the values of a color from hsv color space to rgb.
#' @param h Hue of the color
#' @param s Saturation of the color
#' @param v Value of the color
#' @keywords hsv rgb
#' @export
#' @examples
#' hsv2rgb(150,0.2,0.6)
hsv2rgb <- function(h,s,v) { 
  h <- h/60
  chroma <- v*s
  X <- chroma * (1 - abs((h %% 2) - 1))
  i <- floor(h)

  rgbVector <- switch(as.character(i),  
    "0"=c(chroma,X,0),
    "1"=c(X,chroma,0),
    "2"=c(0,chroma,X),
    "3"=c(0,X,chroma),
    "4"=c(X,0,chroma),
    "5"=c(chroma,0,X)
  )
  m <- v - chroma
  rgbVector <- rgbVector + m

  return(rgb(rgbVector[1],rgbVector[2],rgbVector[3]))
}

#' complement
#'
#' This function returns the complement color of a rgb color
#' @param hex The base color specified as hex
#' @param typeVal Can be specified as split or double. Default is empty.
#' @keywords complement
#' @export
#' @examples
#' complement("#121314")
complement <- function(hex, typeVal="") {
  count <- 1
  rotation <- 180
  scope <- 0

 if(typeVal == "split") {
   count <- 3
   rotation <- 180
   scope <- 180
 } else if(typeVal == "double") {
   count <- 5
   rotation <- 180
   scope <- 180
 }
 
 return(rotational_dispersion(hex,count,"hue",scope,rotation))
}

#' triadic
#'
#' This function returns triadic colors to a given hex color
#' @param hex The base color specified as hex
#' @keywords triadic
#' @export
#' @examples
#' triadic("#121314")
triadic <- function(hex) {
  return(rotational_dispersion(hex,3,"hue",360,0))
}

#' tetradic
#'
#' This function returns tetradic colors to a given hex color
#' @param hex The base color specified as hex
#' @keywords tetradic
#' @export
#' @examples
#' tetradic("#121314")
tetradic<- function(hex) {
  return(rotational_dispersion(hex,4,"hue",360,0))
}

#' pentadic
#'
#' This function returns pentadic colors to a given hex color
#' @param hex The base color specified as hex
#' @keywords pentadic
#' @export
#' @examples
#' pentadic("#121314")
pentadic <- function(hex) {
  return(rotational_dispersion(hex,5,"hue",360,0))
}

degrees <- function(degrees, offset) {
  degrees <- degrees + 360 + offset
  return(degrees %% 360)
}

# pick a point on the wheel, the number of degrees either side to cover and the split
rotational_dispersion <- function(hex,count,typeVal,scope,rotation) {
  hexcol <- col2rgb(hex)
  hsv <- rgb2hsv(
           r=hexcol["red",],
           g=hexcol["green",],
           b=hexcol["blue",]
         )
  h <- hsv["h",]*360
  s <- hsv["s",]*100
  v <- hsv["v",]*100
  palette <- c()
  offset <- 0
  i <- 0

  # if scope is 360, the start and end point are the same color, so should be avoided
  steps <- scope / (count - 1)
  if(typeVal=="hue" && (scope == 360 || scope == 0)){
    steps <- scope / count
  }

  # if scope is 360, start on the current color
  origin <- 0
  if(scope != 360) {
    origin <- degrees(degrees(h, rotation), -1 * scope / 2)
  }

  for(i in 0:(count-1)) {
    offset <- steps * i
    if(typeVal=="hue") {
      palette <- c(palette, c(hsv2rgb(degrees(origin,offset)%%360,s/100,v/100)))
    } else if(typeVal=="saturation") {
      palette <- c(palette, c(hsv2rgb(h,offset/100,v/100)))
    } else if(typeVal=="value" || typeVal == "lightness" || typeVal == "brightness") {
      palette <- c(palette, hsv2rgb(h,s/100,offset/100))
    }
  }
  return(palette)
}

#' Monochromatic
#'
#' This function generates a specified number of monochromatic colors for a given base color
#' @param baseColor The base color specified as hex
#' @param count Number of colors the palette should contain
#' @keywords monochromatic color
#' @export
#' @examples
#' generateMonochromaticColors("#121314", 5)
generateMonochromaticColors <- function(baseColor, count) {
  if(count < 1) {
    return(c())
  } else if (count < 2) {
    return(c(baseColor))
  }

  count <- floor(count)
  base <- baseColor
  colors <- c(baseColor) 
  rgbColor <- col2rgb(baseColor)
  multi <- 2 * 1 / count

  for(i in 1:floor(count/2)){
    colors <- c(
      rgb(((255-rgbColor["red",])*multi*i)/255, 
      ((255-rgbColor["green",])*multi*i)/255, 
      ((255-rgbColor["blue",])*multi*i)/255),
      colors, rgb(
      (rgbColor["red",]*multi*i)/255, 
      (rgbColor["green",]*multi*i)/255, 
      (rgbColor["blue",]*multi*i)/255
      )
    )
  }

  return(tail(colors, n=count))
}

#' complement colors
#'
#' This function generates a color plate with the complement color
#' @param baseColor The base color specified as hex
#' @param count Number of colors the palette should contain
#' @keywords complement color
#' @export
#' @examples
#' complementColors("#121314", 5)
complementColors <- function(baseColor, count) {
  return(c(generateMonochromaticColors(baseColor, ceiling(count/2)), generateMonochromaticColors(complement(baseColor), floor(count/2))))
}

#' triadic colors
#'
#' This function generates a specified number of triadic colors for a given base color
#' @param baseColor The base color specified as hex
#' @param count Number of colors the palette should contain
#' @keywords triadic color
#' @export
#' @examples
#' triadicColors("#121314", 5)
triadicColors <- function(baseColor, count) {
  return(head(unique(c(generateMonochromaticColors(triadic(baseColor)[1], ceiling(count/3)), generateMonochromaticColors(triadic(baseColor)[2], round(count/3)), generateMonochromaticColors(triadic(baseColor)[3], ceiling(count/3)))), n=count))
}

#' tetradic colors
#'
#' This function generates a specified number of tetradic colors for a given base color
#' @param baseColor The base color specified as hex
#' @param count Number of colors the palette should contain
#' @keywords tetradic color
#' @export
#' @examples
#' tetradicColors("#121314", 5)
tetradicColors <- function(baseColor, count) {
  return(head(unique(c(generateMonochromaticColors(tetradic(baseColor)[1], ceiling(count/3)), generateMonochromaticColors(tetradic(baseColor)[2], ceiling(count/3)), generateMonochromaticColors(tetradic(baseColor)[3], ceiling(count/3)), generateMonochromaticColors(tetradic(baseColor)[4], ceiling(count/3)))), n=count))
}

