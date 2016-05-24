SplitSpLines <-
function( sp.lines, into ) {
  # This function divides the sp.lines object into [into] sub sets of Spatial Lines
  # Objects
  #
  # Args:
  #   sp.lines: Object of class SpatialLines calculated by the function Df2SpLines.
  #   into: Number of times that the sp.lines object must be divided
  #   
  # Returns:
  #   A list of size [into] containing Spatial Lines Objects
  
  size <- length(sp.lines)
  
  if (size <= 1)
    stop("The length of the Spatial Lines object must be greater than 1")
  
  if (into >= size)
    stop("Error!")
  
  sp.list <- list()
  
  interval <- length(sp.lines) %/% into
  
  # Serial Execution
  count <- 1
  
  for (i in seq(from=1, to=length(sp.lines), by=interval)) {
    if (count != into){
      sp.list <- c(sp.list, (sp.lines[i:(i + interval - 1)]))
    } else {
      sp.list <- c(sp.list, (sp.lines[i:length(sp.lines)]))
      break
    }
    count <- count + 1
  }
  
  sp.list
}
