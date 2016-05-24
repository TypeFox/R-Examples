################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG
# 09.01.2016: Added more attributes to result.
# 30.11.2015: More efficient implementation. Added attributes.
# 29.08.2015: Added importFrom.
# 15.12.2014: Changed parameter names to format: lower.case
# 09.12.2014: Moved from PCRsim.

#' @title Height To Peak.
#'
#' @description
#' Internal helper function to convert a peak into a plottable polygon.
#'
#' @details
#' Converts a single height and size value to a plottable 0-height-0 triangle/peak value.
#' Makes 3 data points from each peak size for plotting a polygon representing a peak.
#' Factors in other columns might get converted to factor level.
#' 
#' @param data data frame containing at least columns 'Height' and 'Size'.
#' @param width numeric specifying the width of the peak in bp.
#' @param keep.na logical, TRUE keep empty markers.
#' 
#' @return data.frame with new values.
#' 
#' @export
#' 
#' @importFrom utils str
#' 
#' @keywords internal

heightToPeak <- function(data, width=1, keep.na=TRUE, debug=FALSE){
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("data:")
    print(str(data))
    print(head(data))
    print(tail(data))
    print("width:")
    print(width)
    print("keep.na:")
    print(keep.na)
  }
  
  # CHECK ARGUMENTS -----------------------------------------------------------
  
  if(!is.data.frame(data)){
    stop("'data' must be a data frame.",
         call. = TRUE)
  }
  
  if(!any(grepl("Height",names(data)))){
    stop("Data frame 'data' must contain a column 'Height'.",
         call. = TRUE)
  }
  
  if(!any(grepl("Size",names(data)))){
    stop("Data frame 'data' must contain a column 'Size'.",
         call. = TRUE)
  }
  
  if(!is.numeric(width)){
    stop("Peak width 'width' must be a numeric value.",
         call. = TRUE)
  }
  
  if(!is.logical(keep.na)){
    stop("'keep.na' must be a logical value.",
         call. = TRUE)
  }
  
  if("Size" %in% names(data)){
    if(!is.numeric(data$Size)){
      data$Size <- as.numeric(data$Size)
      message("'Size' must be numeric. Data converted!")
    }
  }
  
  # FUNCTION ------------------------------------------------------------------
  
  if(!keep.na){
    # Remove all rows with no height.
    data <- data[!is.na(data$Height),]
  } else if (keep.na){
    # Replace all NAs with 0s.
    data$Height[is.na(data$Height)] <- 0
  }
  
  # Calculate coordinates for plotting peaks ----------------------------------
  
  # Size should be: [size-x], [size], [size+x].
  nb <- length(data$Size)
  vsize <- vector(mode = "numeric", length = nb * 3)
  vsize[seq(1, nb * 3, 3)] <- data$Size - width/2
  vsize[seq(2, nb * 3, 3)] <- data$Size
  vsize[seq(3, nb * 3, 3)] <- data$Size + width/2
  
  # Height should be: 0, [height], 0.
  vheight <- vector(mode = "numeric", length = nb * 3)
  vheight[seq(1, nb * 3, 3)] <- 0
  vheight[seq(2, nb * 3, 3)] <- data$Height
  vheight[seq(3, nb * 3, 3)] <- 0
  
  # Stretch data frame 3 times (repeat each row).
  data <- data[rep(seq_len(nrow(data)), each=3),]

  # Add new data.  
  data$Size <- vsize
  data$Height <- vheight
  
  # Add attributes to result.
  attr(data, which="heightToPeak, strvalidator") <- as.character(utils::packageVersion("strvalidator"))
  attr(data, which="heightToPeak, call") <- match.call()
  attr(data, which="heightToPeak, date") <- date()
  attr(data, which="heightToPeak, data") <- substitute(data)
  attr(data, which="heightToPeak, width") <- width
  attr(data, which="heightToPeak, keep.na") <- keep.na
  
  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
    print("data:")
    print(str(data))
    print(head(data))
    print(tail(data))
  }

  return(data)
  
}
