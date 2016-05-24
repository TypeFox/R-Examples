# Purpose        : Color conversion functions
# Maintainer     : Pierre Roudier (pierre.roudier@landcare.nz);
# Contributions  : Tomislav Hengl (tom.hengl@wur.nl); Dylan Beaudette (debeaudette@ucdavis.edu); 
# Dev Status     : Alpha
# Note           : Converts from HEX to KML and back;


## KML colour follow the scheme #aabbggrr
col2kml <- function(colour){

  # Getting the HEX code
  res <- rgb(t(col2rgb(colour, alpha = TRUE))/255, alpha = TRUE)
  # Converting from HEX to KML
  res <- hex2kml(res)

  res
}

## Convert hex colours to KML specs
hex2kml <- function(hex){

  res <- aaply(hex, 1, function(hex){
    if (str_length(hex) == 9) # if alpha is present
      res <- paste("#", str_sub(tolower(hex), 8, 9), str_sub(tolower(hex), 6, 7), str_sub(tolower(hex), 4, 5), str_sub(tolower(hex), 2, 3), sep = "")
    else
      res <- paste("#ff", str_sub(tolower(hex), 6, 7), str_sub(tolower(hex), 4, 5), str_sub(tolower(hex), 2, 3), sep = "")
    res
    }
  )

  as.vector(res)
}

## KML spec is alpha-BGR
kml2hex <- function(kml) {

  res <- aaply(kml, 1, function(x){
    res <- paste("#", str_sub(toupper(x), 8, 9), str_sub(toupper(x), 6, 7), str_sub(toupper(x), 4, 5), str_sub(toupper(x), 2, 3), sep = "")
    res
    }
  )

  as.vector(res)
}

## From Munsell colour codes to FBGR codes:
munsell2kml <- function(
  the_hue, 
  the_chroma,
  the_value, 
  alpha = 1
  ){

  res <- hex2kml(aqp::munsell2rgb(the_hue, the_value, the_chroma, alpha=alpha))
  
  return(res)
}

# end of script;
