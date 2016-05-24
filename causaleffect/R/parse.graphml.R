parse.graphml <-
function(file, format = c("standard", "internal"), nodes = c(), use.names = TRUE) {
  format <- match.arg(format)
  res <- switch(format, standard = parse.graphml.standard(file, nodes, use.names),
                   internal = parse.graphml.internal(file, nodes, use.names),
                   stop(paste("Unknown file format:", format)))
  return(res)  
}
