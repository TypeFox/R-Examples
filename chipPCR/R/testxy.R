testxy <- function(x, y, txt.x = "Enter abscissa value", 
                   txt.y = "Enter ordinate value", both = TRUE, length = TRUE) {
  
  if (both == TRUE) {
    # Test if x and y exist and have identical lengths.
    if (is.null(x)) 
      stop(txt.x)
    if (is.null(y)) 
      stop(txt.y)
    #   if (is.numeric(x) ) 
    #     stop("Abscissa value must have numeric class")
    #   if (is.numeric(y)) 
    #     stop("Ordinate value must have numeric class")
    if (length)
      if (length(x) != length(y)) 
        stop("Use abscissa and ordinate data with same number of 
  	elements")
  } else {
    # Test if x and y exist and have identical lengths.
    if (is.null(y)) 
      stop(txt.y)
  }
}