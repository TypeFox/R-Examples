as.CairoPath <- function(path)
{
  x <- lapply(path, function(element) {
    type <- attr(element, "type")
    if (is.null(type))
      stop("Path element must have an integer 'type' attribute")
    as.integer(element)
  })
  class(x) <- "CairoPath"
  x
}

as.CairoGlyph <-
function(x)
{
  x <- as.struct("CairoGlyph", c("index", "x", "y"))
  x[[1]] <- as.numeric(x[[1]])
  x[[2]] <- as.numeric(x[[2]])
  x[[3]] <- as.numeric(x[[3]])
  
  return(x)
}

as.CairoTextCluster <-
  function(x)
{
  x <- as.struct("CairoTextCluster", c("num_bytes", "num_glyphs"))
  x[[1]] <- as.integer(x[[1]])
  x[[2]] <- as.integer(x[[2]])
  
  return(x)
}
