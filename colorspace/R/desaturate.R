## desaturate colors (remove chroma in HCL space)
desaturate <- function(col) {
  ## col has to be hex code, otherwise col2rgb is used
  if(is.character(col) &&
    (all(substr(col, 1L, 1L) == "#") & all(nchar(col) %in% c(7L, 9L))))
  {
    ## extract alpha from hex (if any)
    alpha <- substr(col, 8L, 9L)
    ## retain only RGB in hex
    col <- substr(col, 1L, 7L)
    ## convert to colorspace::RGB
    col <- hex2RGB(col)
  } else {
    col <- col2rgb(col, alpha = TRUE)
    ## extract alpha values (if non-FF)
    alpha <- format(as.hexmode(col[4L, ]), width = 2L, upper.case = TRUE)
    alpha[alpha == "FF"] <- ""
    ## retain only RGB
    col <- RGB(t(col[1L:3L, ])/255)
  }
  
  ## convert to HCL and remove chroma
  col <- as(col, "polarLUV")
  col@coords[, 2L] <- 0
  
  ## fix-up extreme luminance cases
  col@coords[col@coords[, 1L] <= 0 | col@coords[, 1L] >= 100, 2L:3L] <- 0
  
  ## convert back to hex and add alpha again (if any)
  col <- hex(col)
  col <- paste(col, alpha, sep = "")
  return(col)
}
