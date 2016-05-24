# Purpose        : List of available aesthetics is given .all_kml_aesthetics along with their default values
# Maintainer     : Pierre Roudier (pierre.roudier@landcare.nz);
# Contributions  : Tomislav Hengl (tom.hengl@wur.nl); Dylan Beaudette (debeaudette@ucdavis.edu); 
# Dev Status     : Pre-Alpha
# Note           : Functionality for constant transparency under development;


.all_kml_aesthetics <- list(
  colour = "black",
  fill = "white",
  shape = paste(get("home_url", envir = plotKML.opts), get("icon", envir = plotKML.opts), sep=""), 
  whitening = "",
  alpha = 1,
  size = get("LabelScale", envir = plotKML.opts),
  width = 1,
  labels = "",
  altitude = 0,
  altitudeMode = "",  
  balloon = FALSE
)

# Parsing a call
#
.parse_call_for_aes <- function(call) {
  called_options <- names(call)
  ind_aes <- charmatch(called_options, names(.all_kml_aesthetics))

  names(.all_kml_aesthetics)[ind_aes[!is.na(ind_aes)]]
}

# Applying aesthetics
#
kml_aes <- function(obj, ...) {

  # Getting parent call
  parent_call <- substitute(list(...))
  parent_call <- as.list(parent_call)[-1]

  # Parse the current call
  called_aes <- .parse_call_for_aes(parent_call)
  aes <- list()

  # Names
  if ("labels" %in% called_aes) {
    # If labels defined using a column of data
    if (is.name(parent_call[['labels']])){
      aes[['labels']] <- as.character(obj[[as.character(parent_call[['labels']])]])
    }
    # if labels given as a vector
    else {
      labels <- eval(parent_call[['labels']])
      if (length(labels) == length(obj))
        aes[['labels']] <- as.character(labels)
      else
        aes[['labels']] <- rep(as.character(labels), length.out = length(obj))
    }
  }
  else {
    if ("data" %in% slotNames(obj)) {

      # If only one data column is represented, we use its values as labels
      if (length(called_aes) == 1) {
        # If its the name of a column
        if (is.name(parent_call[[called_aes]])) {
          # If the column is numeric
          if (is.numeric(obj[[as.character(parent_call[[called_aes]])]]))
            aes[['labels']] <- as.character(round(obj[[as.character(parent_call[[called_aes]])]], digits = 3))
          # Otherwise
          else
            aes[['labels']] <- as.character(obj[[as.character(parent_call[[called_aes]])]])
        }
        # If its a call we have to evaluate it first
        else if (is.call(parent_call[[called_aes]]))
          aes[['labels']] <- as.character(format(eval(parent_call[[called_aes]], obj@data), digits = 3))
        # default behaviour is just numbering using the N first integers
        else
          aes[['labels']] <- as.character(1:length(obj))
      }
      else
        aes[['labels']] <- rownames(obj@data)
    }
    else
      # default behaviour is just numbering using the N first integers
      aes[['labels']] <- as.character(1:length(obj))
  }

  # Colour
  if ("colour" %in% called_aes) {

    # If a column name as been used
    if (is.name(parent_call[['colour']]) | is.call(parent_call[['colour']]) & "data" %in% slotNames(obj)) {

      # Trying to get the colour ramp
      if ("colour_scale" %in% names(parent_call)){
        if("z.lim" %in% names(parent_call)){
          aes[['colour']] <- kml_colour(obj, colour = parent_call[['colour']], colour_scale = parent_call[['colour_scale']], z.lim = eval(parent_call[['z.lim']]))
        } else {
          aes[['colour']] <- kml_colour(obj, colour = parent_call[['colour']], colour_scale = parent_call[['colour_scale']])
        }
      }
      # Otherwise use the default colour ramp
      else {
        if("z.lim" %in% names(parent_call)){
        aes[['colour']] <- kml_colour(obj, colour = parent_call[['colour']], z.lim = eval(parent_call[['z.lim']]))
        } else {
        aes[['colour']] <- kml_colour(obj, colour = parent_call[['colour']])         
        }
    }
    }

    # Otherwise it is interpreted as a colour to use
    else {
      aes[['colour']] <- rep(col2kml(parent_call[['colour']]), length.out = length(obj))
    }
  }
  # using the default value
  else {
    aes[['colour']] <- rep(col2kml(.all_kml_aesthetics[["colour"]]), length.out = length(obj))
  }

  # Alpha
  if ("alpha" %in% called_aes) {
    aes[['colour']] <- kml_alpha(obj, alpha  = eval(parent_call[['alpha']], obj@data), colours = aes[['colour']])
  }
#   else {
#     aes[['alpha']] <- rep(.all_kml_aesthetics[["alpha"]], length.out = length(obj))
#   }

  # Whitening
#     if ("whitening" %in% called_aes) {
#       aes[['colour']] <- kml_whitening(obj, whitening, aes[['colour']], ...)
#     }

  # Shape
  if ("shape" %in% called_aes) {
    aes[["shape"]] <- kml_shape(obj, ...)
  }
  else {
    aes[["shape"]] <- rep(.all_kml_aesthetics[["shape"]], length.out = length(obj))
  }

  # Size
  if ("size" %in% called_aes) {
    # If a column name has been used
    if (is.name(parent_call[['size']])){
      aes[['size']] <- kml_size(obj, size = as.character(parent_call[['size']]))
    }
    # Otherwise it is interpreted as a vector
    else {
      aes[['size']] <- rep(parent_call[['size']], length.out = length(obj))
    }
  }
  else {
    aes[['size']] <- rep(.all_kml_aesthetics[["size"]], length.out = length(obj))
  }

  # Width
  if ("width" %in% called_aes) {
   # If a column name has been used
   if(is.call(parent_call[['width']])){
      aes[['width']] <- eval(parent_call[['width']])    
    } else {
    if (is.name(parent_call[['width']])){
      aes[['width']] <- kml_width(obj, width = as.character(parent_call[['width']]))
    }
    # Otherwise it is interpreted as a vector
    else {
      aes[['width']] <- rep(parent_call[['width']], length.out = length(obj))
    }
  }}
  else {
    aes[['width']] <- rep(.all_kml_aesthetics[["width"]], length.out = length(obj))
  }
  
  # Altitude
  if ("altitude" %in% called_aes) {
    aes[['altitude']] <- kml_altitude(obj, altitude = eval(parent_call[['altitude']], obj@data))
  }
  else {
    aes[['altitude']] <- kml_altitude(obj)
  }                 

  # AltitudeMode
  if ("altitudeMode" %in% called_aes) {
      aes[["altitudeMode"]] <- parent_call[['altitudeMode']]
  }
  else {
      aes[["altitudeMode"]] <- kml_altitude_mode(aes[['altitude']]) 
  }

  # Balloon (pop ups)
  if ("balloon" %in% called_aes) {
    aes[['balloon']] <- eval(parent_call[['balloon']])
  }
  else {
    aes[['balloon']] <- .all_kml_aesthetics[["balloon"]]
  }

  aes
}


.getColourScale <- function(data, z.lim, colour_scale = NULL) {

  ## default colour palettes
  .colour_scale_numeric = get("colour_scale_numeric", envir = plotKML.opts)
  .colour_scale_factor = get("colour_scale_factor", envir = plotKML.opts)

  if (is.null(colour_scale)) {
    # If data is numeric
    if (is.numeric(data)) {
      colour_scale <- .colour_scale_numeric
    }
    # If data is a factor
    else {
      colour_scale <- .colour_scale_factor
    }
  }
  else {
    colour_scale <- eval(colour_scale)
  }

  # creates pal function
  pal <- colorRamp(colour_scale, space = "rgb", interpolate = "linear")

  if (is.numeric(data)) {
    if(missing(z.lim)) { z.lim <- range(data, na.rm = TRUE, finite = TRUE) }
    data <- scales::rescale(data, from=z.lim)
    data <- ifelse(data<0, 0, ifelse(data>1, 1, data))
    cols <- rep("#FFFFFF", length(data))
    cols[!(is.na(data)|is.nan(data))] <- rgb(pal(data[!(is.na(data)|is.nan(data))]) / 255)
  }
  
  # factor variable:
  else {
    # if it is not numeric, it must be a factor:
    data <- as.factor(data)
    values <- 1:length(levels(data))
    cols <- rep("#FFFFFF", length(values))    
    values <- scales::rescale(values) # putting values between 0 and 1
    cols[!(is.na(values)|is.nan(values))] <- rgb(pal(values[!(is.na(values)|is.nan(values))]) / 255)
    # In case of a factor variable, reclassify each level by its corresponding colour:
    levels(data) <- cols
    cols <- as.character(data)
  }

  cols
}


# Colour (points, polygons, lines, raster)
##----------------
kml_colour <- function(obj, colour, colour_scale = NULL, ...){

  # Getting the vector of values to scale
  x <- eval(colour, envir = obj@data) 

  # Retrieving colour scale
  cols <- .getColourScale(data = x, colour_scale = colour_scale, ...)
  cols <- col2kml(cols)

  return(cols)
}


## Shape (points)
##----------------
kml_shape <- function(obj, shape, ...){
  # Simple implementation if a URL is given
  rep(shape, length.out = length(obj))
}

# modification of the colours using that alpha value
.includeAlphaInColourRamp <- function(colours, alpha, RGBA = FALSE) {
  # colours is a vector of characters describing the colour ramp values
  # alpha is a value between 0 and 255
  # RGBA is set to flase if for kml colour, true for RGBA colours

  # Conversion to hex mode
  alpha <- sprintf("%X", alpha)

  colours <- str_sub(colours, 2, str_length(colours))

  if (RGBA)
    b <- cbind("#", colours, alpha)
  else
    b <- cbind('#', alpha, aaply(colours, 1, function(x) str_sub(x, 3, 9)))

  res <-  apply(b, 1, function(x) paste(x, collapse = ""))

  if (RGBA)
    res <- toupper(res)
  else
    res <- tolower(res)

  res
}

# Opacity (points, polygons, lines, raster)
#
# This function modifies the vector of KML colours
kml_alpha <- function(obj, alpha, colours, RGBA = FALSE, ...){

  # Alpha is constant or continuous
  if (is.numeric(alpha)) {
    # rescaling data if continuous case
    if (length(alpha) > 1)
      alpha <- scales::rescale(alpha)

    alpha <- round(255*alpha, digits = 0)
    cols <- .includeAlphaInColourRamp(colours, alpha, RGBA)

#     # Constant transparency
#     if (length(alpha) == 1) {
#       # modification of the KML colours using that alpha value
#       if (RGBA) {
#         # Conversion to hex mode
#         alpha <- sprintf("%X", alpha)
#         # Adding hex A (alpha) to a RGB hex string
#         cols <- aaply(colours, 1, function(x) paste(x, alpha, sep = ""))
#       }
#       else {
#         # Conversion to hex mode
#         alpha <- sprintf("%x", alpha)
#         cols <- aaply(colours, 1, function(x)  paste("#", alpha, str_sub(x, 4, 9), sep = ""))
#       }
#     }
#     # transparency as an aesthetic
#     else {

#     }
  }
  # Categorical data
  else {

    if (is.numeric(obj)) {
      limits <- range(obj, na.rm = TRUE, finite = TRUE)
      brks <- seq(limits[1], limits[2], length.out = length(colours))
      grps <- cut(obj, breaks = brks, include.lowest = TRUE)
    }
    else
      stop("Transparency is only available for numeric data.")
  }

  cols
}

# Whitening of a given coloured layer (update of the colour vector)
# kml_whitening <- function(obj, whitening, col.vect){
# 
# }

# Size (points)
kml_size <- function(obj, size, size.min = get("size_range", envir = plotKML.opts)[1], size.max = get("size_range", envir = plotKML.opts)[2], size.default = 1){

  if (!is.na(size) & "data" %in% slotNames(obj)) {
    ## If data is numeric
    if (is.numeric(obj[[size]])) {
      max.value <- max(obj[[size]], na.rm = TRUE)
      size.values <- size.min + scales::rescale(obj[[size]], from = range(obj[[size]], na.rm = TRUE)) * size.max
      ## fix the missing values:
      size.values[is.na(size.values)] <- size.default
    }
    # Otherwise: factor, character, logical, ...
    else {
      if (!is.factor(obj[[size]])){
        obj[[size]] <- factor(obj[[size]])
      }

      ## compute number of levels
      nl <- nlevels(obj[[size]])
      ## compute the different size values
      sizes.levels <- seq(size.min, size.max, length.out = nl)
      ## assign them to different classes:
      sizes.values <- cut(as.integer(obj[[size]]), breaks = seq(1, nl + 1, by=1), labels = sizes.levels, include.lowest = TRUE)
      size.values <- as.numeric(as.character(sizes.values))
      ## fix the missing values:
      size.values[is.na(size.values)] <- size.default
    }
    
  } else {
    ## If no size aesthetic is asked, or if no data slot
    size.values <- rep(size.default, length.out = length(obj))
  }

  return(size.values)
}

# Width (lines)
kml_width <- function(obj, width, width.min = 1, width.max = 6, width.default = 1){

  if (!is.na(width) & "data" %in% slotNames(obj)) {
    # If data is numeric
    if (is.numeric(obj[[width]])) {
      max.value <- max(obj[[width]], na.rm = TRUE)
      width.values <- round(width.min + scales::rescale(obj[[width]]) * width.max, 0)
    }
    # Otherwise: factor, character, logical, ...
    else {
      if (!is.factor(obj[[width]]))
        obj[[width]] <- factor(obj[[width]])

      # compute number of levels
      nl <- nlevels(obj[[width]])
      # compute the different size values
      width.levels <- seq(width.min, width.max, length.out = nl)
      # affect them to the factor
      width.values <- cut(obj[[width]], breaks = quantile(obj[[width]], probs=seq(0, 1, length.out = nl + 1)), labels = width.levels, include.lowest = TRUE)
      width.values <- as.numeric(as.character(width.values))
    }
  }
  # If no size aesthetic is asked, or if no data slot
  else
    width.values <- rep(width.default, length.out = length(obj))

  width.values  
}

# end of script;