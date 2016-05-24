centroidtext <- function(polygon, poly.name = NULL, counter = "NA", cex = 1, ...) 
{
  pos <- centroidpos(polygon)		
  if(is.null(poly.name))
    txt <- paste(counter)
  else
    txt <- poly.name
  text(pos[1L], pos[2L], txt, cex = cex, ...)

  return(invisible(NULL))
}

