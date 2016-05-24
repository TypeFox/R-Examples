# generate a MD5 hash using a line's start / end coordinates 
# after rounding to a fixed precision, related to snapping tolerances
generateLineHash <- function(x, precision=-1, algo='murmur32') {
  the.lines <- slot(x, 'lines')
  # sanity check:
  # (there should only be 1 / Line object)
  lines.per.Line <- sapply(the.lines, function(i) length(i@Lines))
  if(any(lines.per.Line > 1))
    stop('more than 1 line per segment, why?', call. = FALSE)
  
  # unzip lines into single list
  the.coords <- lapply(the.lines, function(i) slot(slot(i, 'Lines')[[1]], 'coords'))
  
  # generate digest from start / end vertices, rounded to defined precision
  res <- sapply(the.coords, function(i) {
    n <- nrow(i)
    start.coord <- i[1, ]
    end.coord <- i[n, ]
    rounded.coords <- round(c(start.coord, end.coord), digits = precision)
    hash <- digest(rounded.coords, algo=algo)
    return(hash)
  })
  
  # check for collisions in hash
  if(any(table(res) > 1))
    stop('collision in hash function, consider increasing precision', call. = FALSE)
  
  return(res)
}
