centroidpos <- function(polygon) 
{
  polygon <- na.omit(polygon)
  p <- polygon
  np <- (nrow(p) - 1L)
  if(is.na(p[1L, 1L])) {
    p <- p[2L:(np + 1L),]
    np <- np - 1L
  }
  if((p[1L, 1L] != p[(np + 1L), 1L]) || (p[1L, 2L] != p[(np + 1L), 2L]))
    p[(np + 1L),] <- p[1L,]
	out <- cpos(p, np)

  return(out)
}

